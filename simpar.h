#ifndef simpar_h
#define simpar_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.01

#define GRID_MIN_X 0
#define GRID_MAX_X 1
#define GRID_MIN_Y 0
#define GRID_MAX_Y 1

typedef int bool;

#define TRUE 1
#define FALSE 0

struct point2_t
{
    double x;
    double y;
};
typedef struct point2_t point2_t;

struct vector2_t
{
    double x;
    double y;
};
typedef struct vector2_t vector2_t;

struct particle_t 
{
    point2_t position;  /**< position in 2D space. */
    vector2_t velocity; /**< speed in 2D space. */
    vector2_t force;    /**< scalar acceleration. */
    double mass;        /**< mass in kg. */
};
typedef struct particle_t particle_t;

struct cell_t
{
    point2_t center;    /**< center of mass of the cell. */
    double mass;        /**< total mass of the cell. */
};
typedef struct cell_t cell_t;

struct grid_t
{
    long ncside;
    long size;
    cell_t * cells;
};
typedef struct grid_t grid_t;

grid_t * create_grid(long ncside);
void update_grid(grid_t * grid, particle_t * p, long long n_part);
long get_cell_index(grid_t * g, point2_t p);
bool is_cell_empty(cell_t c);

double calculate_force_magnitude(particle_t pA, particle_t pB);
double calculate_distance(point2_t pA, point2_t pB);

void display_particles(long long n_part, particle_t * p);
void init_particles(long seed, long ncside, long long n_part, particle_t * p);
void update_particle(long long index, particle_t * p, double dt);
void update_forces(grid_t * g, particle_t * p, long long n_part);
void update_forces_brute_force(grid_t * g, particle_t * p, long long n_part);

point2_t center_of_mass(long long n_part, particle_t * p);

#endif /* simpar_h */