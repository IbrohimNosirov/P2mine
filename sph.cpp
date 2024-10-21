#include <cstdio>
#include <cmath>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <omp.h>

#include <fstream>
#include "vec3.hpp"
#include "io.hpp"
#include "params.hpp"
#include "state.hpp"
#include "binhash.hpp"
#include "interact.hpp"
#include "leapfrog.hpp"

#include "timing_harness.hpp"

typedef int (*domain_fun_t)(float, float, float);

int box_indicator(float x, float y, float z)
{
    return (x < 0.5) && (y < 0.75) && (z < 0.5);
}

int circ_indicator(float x, float y, float z)
{
    float dx = (x-0.5);
    float dy = (y-0.5);
    float dz = (z-0.5);
    float r2 = dx*dx + dy*dy + dz*dz;
    return (r2 < 0.25*0.25*0.25);
}

sim_state_t* place_particles(sim_param_t* param, domain_fun_t indicatef)
{
    float h  = param->h;
    float hh = h/1.3;

    // Count mesh points that fall in indicated region.
    int count = 0;
    for (float x = 0; x < 1; x += hh)
        for (float y = 0; y < 1; y += hh)
        	for (float z = 0; z < 1; z += hh)
        		count += indicatef(x,y,z);

    // Populate the particle data structure
    sim_state_t* s = alloc_state(count);
    int p = 0;
    for (float x = 0; x < 1; x += hh) {
        for (float y = 0; y < 1; y += hh) {
            for (float z = 0; z < 1; z += hh) {
                if (indicatef(x,y,z)) {
                    vec3_set(s->part[p].x, x, y, z);
                    vec3_set(s->part[p].v, 0, 0, 0);
                    ++p;
                }
            }
        }
    }
    return s;    
}

void normalize_mass(sim_state_t* s, sim_param_t* param)
{
    s->mass = 1;
    hash_particles(s, param->h);
    compute_density(s, param, param->h);
    float rho0 = param->rho0;
    float rho2s = 0;
    float rhos  = 0;
    for (int i = 0; i < s->n; ++i) {
        rho2s += (s->part[i].rho)*(s->part[i].rho);
        rhos  += s->part[i].rho;
    }
    s->mass *= ( rho0*rhos / rho2s );
}

sim_state_t* init_particles(sim_param_t* param)
{
    sim_state_t* s = place_particles(param, box_indicator);
    normalize_mass(s, param);
    return s;
}

double compute_kinetic_energy(sim_state_t* s)
{
    double total_ke = 0.0;
    int n = s->n;
    particle_t* p = s->part;
    double mass = s->mass;

    #pragma omp parallel for reduction(+:total_ke)
    for (int i = 0; i < n; ++i) {
        float vx = p[i].v[0];
        float vy = p[i].v[1];
        float vz = p[i].v[2];
        double speed_sq = vx*vx + vy*vy + vz*vz;
        double ke = 0.5 * mass * speed_sq;
        total_ke += ke;
    }
    return total_ke;
}

void check_state(sim_state_t* s)
{
    for (int i = 0; i < s->n; ++i) {
        float xi = s->part[i].x[0];
        float yi = s->part[i].x[1];
        float zi = s->part[i].x[2];
        assert( xi >= 0 || xi <= 1 );
        assert( yi >= 0 || yi <= 1 );
        assert( zi >= 0 || zi <= 1 );
    }
}

float compute_max_particle_speed(sim_state_t* s)
{
    float v_max2 = 0.0f; // Keep track of the maximum squared speed
    int n = s->n;
    particle_t* particles = s->part;

    #pragma omp parallel for reduction(max:v_max2)
    for (int i = 0; i < n; ++i) {
        float speed2 = vec3_len2(particles[i].v);
        if (speed2 > v_max2)
            v_max2 = speed2;
    }
    float v_max = sqrtf(v_max2); // Take the square root to get the actual speed
    return v_max;
}

float compute_max_particle_acceleration(sim_state_t* s)
{
    float a_max2 = 0.0f; // Keep track of the maximum squared acceleration
    int n = s->n;
    particle_t* particles = s->part;

    #pragma omp parallel for reduction(max:a_max2)
    for (int i = 0; i < n; ++i) {
        float acceleration2 = vec3_len2(particles[i].a);
        if (acceleration2 > a_max2)
            a_max2 = acceleration2;
    }
    float a_max = sqrtf(a_max2); // Take the square root to get the actual acceleration
    return a_max;
}

int main(int argc, char** argv)
{
    sim_param_t params;
    if (get_params(argc, argv, &params) != 0)
        exit(-1);
    sim_state_t* state = init_particles(&params);
    FILE* fp    = std::fopen(params.fname.c_str(), "w");
    int nframes = params.nframes;
    int npframe = params.npframe;
    float dt    = params.dt;
    int n       = state->n;

    // Initialize variables for predictive rehashing
    int N = 10; // Rehashing interval
    int steps_since_rehash = 0;
    float hash_h = params.h; // Initial hash_h

    double t_start = omp_get_wtime();
    write_header(fp, n, nframes, params.h);
    write_frame_data(fp, n, state, NULL);

    // Initial rehash
    hash_particles(state, hash_h);

    // Use TIMED_CALL to time compute_accel and leapfrog_start
    TIMED_CALL(compute_accel, state, &params, hash_h);
    TIMED_CALL(leapfrog_start, state, dt);

    std::ofstream energy_file("energy_vs_time.out");
    energy_file << "Time KineticEnergy \n";

    check_state(state);
    for (int frame = 1; frame < nframes; ++frame) {
        for (int i = 0; i < npframe; ++i) {
            // Determine whether to rehash
            bool rehash = (steps_since_rehash % N == 0);

            if (rehash) {
                // Predict maximum displacement over N steps
                float v_max = compute_max_particle_speed(state);
                float a_max = compute_max_particle_acceleration(state);

                float d_max = v_max * N * dt +
                    0.5f * a_max * (N * dt) * (N * dt);

                // Update hash_h
                hash_h = params.h + 2.0f * d_max;

                // Rehash particles with updated hash_h
                hash_particles(state, hash_h);
            }

            // Time compute_accel and leapfrog_step
            TIMED_CALL(compute_accel, state, &params, hash_h);
            TIMED_CALL(leapfrog_step, state, dt);
            check_state(state);

            steps_since_rehash++;
        }

        printf("Frame: %d of %d - %2.1f%%\n", frame, nframes,
               100.0f * (float)frame / nframes);

        double ke = compute_kinetic_energy(state);
        double current_time = frame * dt * npframe;

        // Write to file
        energy_file << current_time << " " << ke << "\n";

        // Time write_frame_data if desired
        TIMED_CALL(write_frame_data, fp, n, state, NULL);
    }

    // Close energy file
    energy_file.close();

    double t_end = omp_get_wtime();
    printf("Ran in %g seconds\n", t_end - t_start);

    fclose(fp);
    free_state(state);

    // Write the timing data to a CSV file
    write_timing_data("timing_data.csv");

    return 0;
}
