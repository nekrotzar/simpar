/**
 * @file simpar.h
 * @author Luís Fonseca
 * @date 9 Apr 2019
 * @brief File containing the particle simulation data structures and functions.
 *
 * Both the structs and relevant functions are here defined for usage in
 * the n-body gravitational simulation. Some physical constants are defined
 * as variables. For a detailed overview of the project visit the link below.
 * @see https://github.com/nekrotzar/simpar
 */

#ifndef simpar_h
#define simpar_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND0_1 ((double)random() / ((long long)1 << 31))
#define G 6.67408e-11
#define EPSLON 0.01

#define GRID_MIN_X 0
#define GRID_MAX_X 1
#define GRID_MIN_Y 0
#define GRID_MAX_Y 1

typedef int bool;

#define TRUE 1
#define FALSE 0

/**
 * @brief Point in 2D space.
 * 
 * A geometry point with two coordinates to represent particle positions.
 **/
struct point2_t
{
    double x; /**< x component of the point. */
    double y; /**< y component of the point. */
};
typedef struct point2_t point2_t;

/**
 * @brief Vector in 2D space.
 * 
 * A geometry vector with two coordinates to represent phyisics vectors.
 **/
struct vector2_t
{
    double x; /**< x component of the point. */
    double y; /**< y component of the point. */
};
typedef struct vector2_t vector2_t;

/**
 * @brief Particle.
 * 
 * A particle containing information about its position,
 * velocity and force at a given time-step and its constant mass.
 **/
struct particle_t
{
    point2_t position;  /**< position in 2D space. */
    vector2_t velocity; /**< speed in 2D space. */
    vector2_t force;    /**< scalar acceleration. */
    double mass;        /**< mass in kg. */
};
typedef struct particle_t particle_t;

/**
 * @brief Grid cell in 2D space.
 * 
 * A grid cell containing a point representing the center of mass and the total
 * mass of all the particles inside this cell..
 **/
struct cell_t
{
    point2_t center; /**< center of mass of the cell. */
    double mass;     /**< total mass of the cell. */
};
typedef struct cell_t cell_t;

/**
 * @brief Grid in 2D space.
 * 
 * A wrap around grid limited by the points (0,0) and (1,1) containing
 * ncside * ncside cells.
 **/
struct grid_t
{
    long ncside;
    long size;
    cell_t *cells;
};
typedef struct grid_t grid_t;

/**
 * @brief Initialize particle system with random physical properties.
 *
 * Initialization of an array of type particle_t, using a 
 * given random seed to allow for result replication. Each particle
 * position, velocity and mass are initialized randomly. 
 * 
 * @param seed Seed for the random number generator.
 * @param ncside Size of the grid (number of cells on the side).
 * @param n_part Number of particles.
 * @param p Pointer to an array of particles.
 */
void init_particles(long seed, long ncside, long long n_part, particle_t *p);

/**
 * @brief Initialize 2D grid where particles are located.
 *
 * Initialization of an array of type cell_t, using the 
 * default values for cell position and total mass. Assignment
 * of the grid number of sides and total size.
 * 
 * @param ncside Size of the grid (number of cells on the side).
 * @param g Pointer to a struct of type grid_t.
 **/
void init_grid(long ncside, grid_t *g);

/**
 * @brief Determine the center of mass of each cell.
 *
 * Calculate the center of mass of each cell in the grid by checking
 * each particle position in order to update the center and mass of the
 * cell which the particle belongs to. The function time complexity is
 * O(n + c) where n is the number of particles and c the number of cells.
 * 
 * @param grid Pointer to a struct of type grid_t.
 * @param p Pointer to an array of particles.
 * @param n_part Number of particles.
 */
void update_grid(grid_t *grid, particle_t *p, long long n_part);

/**
 * @brief Get the cell index object.
 * 
 * @param g Pointer to a struct of type grid_t.
 * @param p Particle position to check.
 * @return long The index of the cell where the particle lies.
 */
long get_cell_index(grid_t *g, point2_t p);

/**
 * @brief Check if a cell is empty.
 * 
 * Check if the cell is empty by looking at the total mass
 * of the cell.
 * 
 * @param c A struct of type cell_t.
 * @return true If mass is equal to zero.
 * @return false Mass is not zero.
 */
bool is_cell_empty(cell_t c);

/**
 * @brief Compute the gravitational force applied to each particle.
 * 
 * The force on a particle is calculated only from the centers of masses of
 * its current and all adjacent cells. The function time complexity is
 * O(n * c) where n is the number of particles and c the number of cells.
 * Given this, this is the main target for parallellization.
 * 
 * @param g Pointer to a struct of type grid_t.
 * @param p Pointer to an array of particles.
 * @param n_part Number of particles.
 */
void update_forces(grid_t *g, particle_t *p, long long n_part);

/**
 * @brief The magnitude of the force between particles A and B.
 * 
 * The function calculates the gravitacional force between two particles given
 * the formula: Fa,b = Fb,a = G mA * mB / dA,B ^2
 * 
 * @param pA Particle A.
 * @param pB Partible B.
 * @return double Force magnitude.
 */
double calculate_force_magnitude(particle_t pA, particle_t pB);

/**
 * @brief Calculate distance between two points in 2D space.
 * 
 * The function calculates the distance between two particles given
 * the formula: d = √((pB.x - pA.x)^2 + (pB.y - pA.y)^2)
 * 
 * @param pA Position of A.
 * @param pB Position of B.
 * @return double Distance between two points.
 */
double calculate_distance(point2_t pA, point2_t pB);

/**
 * @brief Calculate the new velocity and then the new position of each particle.
 * 
 * Update each particle in the system state at a given time-step of the simulation.
 * The function time complexity is O(n) where n is the number of particles.
 * 
 * @param p Pointer to an array of particles.
 * @param n_part Number of particles.
 */
void update_particles(particle_t * p, long long n_part);

/**
 * @brief Calculate the particle new values for velocity and position.
 * 
 * The simulation time-step is discretized so each time-step corresponds to 1s.
 * After setting each new position, we check if the new position goes outside the grid
 * and if so determine its new position inside the grid accordingly.
 * 
 * @param index Particle index. 
 * @param p Pointer to an array of particles.
 * @param dt Simulation time-step.
 */
void update_particle(long long index, particle_t *p, double dt);

/**
 * @brief Display the particle system information.
 * 
 * @param p Pointer to an array of particles.
 * @param n_part Number of particles. 
 */
void display_particles(particle_t *p, long long n_part);

/**
 * @brief Determine the center of mass of the particle system.
 * 
 * @param p Pointer to an array of particles.
 * @param n_part Number of particles. 
 * @return point2_t The center of mass 2D point.
 */
point2_t get_center_of_mass(particle_t *p, long long n_part);

#endif /* simpar_h */