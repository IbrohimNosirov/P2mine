#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>
#include <omp.h>

#include "vec3.hpp"
#include "zmorton.hpp"
#include "params.hpp"
#include "state.hpp"
#include "interact.hpp"
#include "binhash.hpp"

// Update density function remains largely the same
inline void update_density(particle_t* pi, particle_t* pj, float h2, float C) {
    float r2 = vec3_dist2(pi->x, pj->x);
    float z = h2 - r2;
    if (z > 0) {
        float rho_ij = C * z * z * z;

        // Use atomic operations to prevent race conditions
        #pragma omp atomic
        pi->rho += rho_ij;
        #pragma omp atomic
        pj->rho += rho_ij;
    }
}

// Function to clear densities
void clear_densities(particle_t* particles, int n) {
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        particles[i].rho = 0.0f;
    }
}

// Function to compute self-density contribution
void add_self_density(particle_t* particles, int n, float self_density) {
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        particles[i].rho += self_density;
    }
}

// Compute density with bucketing
void compute_density_with_bucketing(sim_state_t* s, float h, float h2, float C, float
	self_density) {
    int n = s->n;
    particle_t* particles = s->part;
    particle_t** hash = s->hash;

    // Preallocate buckets array
    std::vector<unsigned> buckets(27);

    #pragma omp parallel for private(buckets)
    for (int i = 0; i < n; ++i) {
        particle_t* pi = &particles[i];
        unsigned num_neighbors = particle_neighborhood(buckets.data(), pi, h);

        for (unsigned j = 0; j < num_neighbors; ++j) {
            particle_t* pj = hash[buckets[j]];

            while (pj != nullptr) {
                if (pi < pj) {
                    update_density(pi, pj, h2, C);
                }
                pj = pj->next;
            }
        }
    }
}

// Compute density without bucketing
void compute_density_without_bucketing(sim_state_t* s, float h2, float C, float
	self_density) {
    int n = s->n;
    particle_t* particles = s->part;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        particle_t* pi = &particles[i];
        for (int j = i + 1; j < n; ++j) {
            particle_t* pj = &particles[j];
            update_density(pi, pj, h2, C);
        }
    }
}

// Main compute density function
void compute_density(sim_state_t* s, sim_param_t* params, BucketingOption
	bucketing_option) {
    int n = s->n;
    particle_t* particles = s->part;

    float h = params->h;
    float h2 = h * h;
    float h3 = h2 * h;
    float h9 = h3 * h3 * h3;
    float mass = s->mass;
    float C = (315.0f / (64.0f * M_PI)) * mass / h9;
    float self_density = (315.0f / (64.0f * M_PI)) * mass / h3;

    clear_densities(particles, n);
    add_self_density(particles, n, self_density);

    if (bucketing_option == BucketingOption::ENABLED) {
        compute_density_with_bucketing(s, h, h2, C, self_density);
    } else {
        compute_density_without_bucketing(s, h2, C, self_density);
    }
}

// Update forces function remains largely the same
inline void update_forces(particle_t* pi, particle_t* pj, float h2, float rho0,
	float C0, float Cp, float Cv) {
    float dx[3];
    vec3_diff(dx, pi->x, pj->x);
    float r2 = vec3_len2(dx);
    if (r2 < h2) {
        const float rhoi = pi->rho;
        const float rhoj = pj->rho;
        float q = sqrtf(r2 / h2);
        float u = 1.0f - q;
        float w0 = C0 * u / (rhoi * rhoj);
        float wp = w0 * Cp * (rhoi + rhoj - 2.0f * rho0) * u / q;
        float wv = w0 * Cv;
        float dv[3];
        vec3_diff(dv, pi->v, pj->v);

        // Equal and opposite pressure forces
        vec3_saxpy(pi->a,  wp, dx);
        vec3_saxpy(pj->a, -wp, dx);
        
        // Equal and opposite viscosity forces
        vec3_saxpy(pi->a,  wv, dv);
        vec3_saxpy(pj->a, -wv, dv);
    }
}

// Initialize accelerations with gravity
void initialize_accelerations(particle_t* particles, int n, float g) {
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        vec3_set(particles[i].a, 0.0f, -g, 0.0f);
    }
}

// Compute accelerations with bucketing
void compute_accel_with_bucketing(sim_state_t* s, float h, float h2, float rho0, float
	C0, float Cp, float Cv) {
    int n = s->n;
    particle_t* particles = s->part;
    particle_t** hash = s->hash;

    // Preallocate buckets array
    std::vector<unsigned> buckets(27);

    #pragma omp parallel for private(buckets)
    for (int i = 0; i < n; ++i) {
        particle_t* pi = &particles[i];
        unsigned num_neighbors = particle_neighborhood(buckets.data(), pi, h);

        for (unsigned j = 0; j < num_neighbors; ++j) {
            particle_t* pj = hash[buckets[j]];

            while (pj != nullptr) {
                if (pi < pj) {
                    update_forces(pi, pj, h2, rho0, C0, Cp, Cv);
                }
                pj = pj->next;
            }
        }
    }
}

// Compute accelerations without bucketing
void compute_accel_without_bucketing(sim_state_t* s, float h2, float rho0,
	float C0, float Cp, float Cv) {
    int n = s->n;
    particle_t* particles = s->part;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        particle_t* pi = &particles[i];
        for (int j = i + 1; j < n; ++j) {
            particle_t* pj = &particles[j];
            update_forces(pi, pj, h2, rho0, C0, Cp, Cv);
        }
    }
}

// Main compute acceleration function
void compute_accel(sim_state_t* s, sim_param_t* params, BucketingOption bucketing_option) {
    // Unpack basic parameters
    const float h = params->h;
    const float rho0 = params->rho0;
    const float k = params->k;
    const float mu = params->mu;
    const float g = params->g;
    const float mass = s->mass;
    const float h2 = h * h;

    int n = s->n;
    particle_t* particles = s->part;

    // Rehash the particles
    hash_particles(s, h);

    // Compute density
    compute_density(s, params, bucketing_option);

    // Initialize accelerations with gravity
    initialize_accelerations(particles, n, g);

    // Constants for interaction term
    const float C0 = 45.0f * mass / (M_PI * powf(h, 6));
    const float Cp = k / 2.0f;
    const float Cv = -mu;

    // Compute accelerations
    if (bucketing_option == BucketingOption::ENABLED) {
        compute_accel_with_bucketing(s, h, h2, rho0, C0, Cp, Cv);
    } else {
        compute_accel_without_bucketing(s, h2, rho0, C0, Cp, Cv);
    }
}

