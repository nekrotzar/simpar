//
//  simpar.h
//  simpar
//
//  Created by Phil Marques on 03/04/2019.
//  Copyright © 2019 instituto superior tecnico. All rights reserved.
//

#ifndef simpar_h
#define simpar_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.0005

typedef struct particle_t{
    double x, y;
    double vx, vy;
    double m;
    long cx, cy;
}particle_t;

typedef struct cell_t{
    double x, y;
    double M;
}cell_t;

double t_mass = 0.0;
double t_cx = 0.0, t_cy = 0.0;

void usg_err(void);
long long val_l(const char* arg);
cell_t** init_grid(const long ncside);
void free_grid(cell_t** g, long ncside);
void init_particles(long seed, long ncside, long long n_part, particle_t *par);
void init_env(cell_t** g, long ncside, particle_t* p, long long n_par);
void accellerate_p(double* ax, double* ay, const cell_t* c, double m, double x, double y);
void update_particles(cell_t** grid, long ncside, particle_t* par, long long n_par, int flag);

#endif /* simpar_h */
