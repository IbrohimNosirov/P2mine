#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <iostream>

#include "vec3.hpp"
#include "zmorton.hpp"

#include "params.hpp"
#include "state.hpp"
#include "interact.hpp"
#include "binhash.hpp"

#define USE_BUCKETING 

inline
void update_density(particle_t* pi, particle_t* pj, float h2, float C)
{
    float r2 = vec3_dist2(pi->x, pj->x);
    float z  = h2-r2;
    if (z > 0) {
        float rho_ij = C*z*z*z;
        #pragma omp atomic 
        pi->rho += rho_ij;
        #pragma omp atomic 
        pj->rho += rho_ij;
    }
}

void compute_density(sim_state_t* s, sim_param_t* params, float hash_h)
{
    int n = s->n;
    particle_t* p = s->part;
    particle_t** hash = s->hash;

    float h  = params->h; // SPH smoothing length remains constant
    float h2 = h * h;
    float h3 = h2 * h;
    float h9 = h3 * h3 * h3;
    float C  = (315.0f / 64.0f / M_PI) * s->mass / h9;

    // Clear densities
    #pragma omp parallel for
    for (int i = 0; i < n; ++i)
        p[i].rho = 0.0f;

#ifdef USE_BUCKETING
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        particle_t* pi = &p[i];
        pi->rho += (315.0f / 64.0f / M_PI) * s->mass / h3;

        unsigned buckets[27];
        unsigned num_neighbors = particle_neighborhood(buckets, pi, hash_h);

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
#else
    for (int i = 0; i < n; ++i) {
        particle_t* pi = s->part+i;
        pi->rho += ( 315.0/64.0/M_PI ) * s->mass / h3;
        for (int j = i+1; j < n; ++j) {
            particle_t* pj = s->part+j;
            update_density(pi, pj, h2, C);
        }
    }
#endif
}

inline
void update_forces(particle_t* pi, particle_t* pj, float h2,
                   float rho0, float C0, float Cp, float Cv)
{
    float dx[3];
    vec3_diff(dx, pi->x, pj->x);
    float r2 = vec3_len2(dx);
    if (r2 < h2) {
        const float rhoi = pi->rho;
        const float rhoj = pj->rho;
        float q = sqrt(r2/h2);
        float u = 1-q;
        float w0 = C0 * u/rhoi/rhoj;
        float wp = w0 * Cp * (rhoi+rhoj-2*rho0) * u/q;
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

void compute_accel(sim_state_t* s, sim_param_t* params, float hash_h)
{
    // Unpack basic parameters
    const float h    = params->h;
    const float rho0 = params->rho0;
    const float k    = params->k;
    const float mu   = params->mu;
    const float g    = params->g;
    const float mass = s->mass;
    const float h2   = h * h;

    // Unpack system state
    particle_t* p = s->part;
    particle_t** hash = s->hash;
    int n = s->n;

    // Compute density and color
    compute_density(s, params, hash_h);

    // Start with gravity and surface forces
    #pragma omp parallel for
    for (int i = 0; i < n; ++i)
        vec3_set(p[i].a,  0.0f, -g, 0.0f);

    // Constants for interaction term
    float C0 = 45.0f * mass / M_PI / (h2 * h2 * h);
    float Cp = k / 2.0f;
    float Cv = -mu;

    // Accumulate forces
#ifdef USE_BUCKETING
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        particle_t* pi = &p[i];
        unsigned buckets[27];
        unsigned num_neighbors = particle_neighborhood(buckets, pi, hash_h);

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
#else
    for (int i = 0; i < n; ++i) {
        particle_t* pi = p+i;
        for (int j = i+1; j < n; ++j) {
            particle_t* pj = p+j;
            update_forces(pi, pj, h2, rho0, C0, Cp, Cv);
        }
    }
#endif
}
