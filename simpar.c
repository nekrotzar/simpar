#include "simpar.h"

int main(int argc, char const *argv[])
{
    if (argc != 5)
    {
        printf("ERROR: Invalid number of arguments (%d).\n\n", (argc - 1));
        printf("Usage - Required exactly 4 arguments i.e ./simpar seed gridsize particles timestep\n");
        printf("\tseed\t\t - (int) : seed for the random number generator\n");
        printf("\tgridsize\t - (int) : size of the grid (number of cells on the side\n");
        printf("\tparticles\t - (int) : number of particles\n");
        printf("\ttimestep\t - (int) : number of time-steps\n");
        exit(1);
    }

    long seed = atoi(argv[1]);
    long ncside = atoi(argv[2]);
    long long n_part = atoi(argv[3]);
    long n_step = atoi(argv[4]);

    particle_t * p = malloc(sizeof(particle_t) * n_part); 
    init_particles(seed, ncside, n_part, p);

    grid_t * g = malloc(sizeof(grid_t));
    init_grid(ncside, g);
    
    for(long step = 0; step < n_step; step++)
    {
        update_grid(g, p, n_part);
        update_forces(g, p, n_part);
        update_particles(p, n_part);
    }

    point2_t c = get_center_of_mass(p, n_part);

    printf("%.2f %.2f\n", p[0].position.x, p[0].position.y);
    printf("%.2f %.2f", c.x, c.y);

    free(p);
    free(g->cells);
    free(g);
    
    return 0;
}

void init_particles(long seed, long ncside, long long n_part, particle_t * p)
{
    long long i;
    
    srandom(seed);

    for(i=0; i < n_part; i++)
    {
        p[i].position.x = RND0_1;
        p[i].position.y = RND0_1;
        p[i].velocity.x = RND0_1 / ncside / 10.0;
        p[i].velocity.y = RND0_1 / ncside / 10.0;
        p[i].force.x = 0.0;
        p[i].force.y = 0.0;
        p[i].mass = RND0_1 * ncside / (G * 1e6 * n_part);
    }
}

void init_grid(long ncside, grid_t * g)
{
    g->ncside = ncside;
    g->size = ncside * ncside;
    g->cells = malloc(sizeof(cell_t) * g->size);

    for(int i = 0; i < g->size; i++)
    {
        g->cells[i].center.x = 0.0;
        g->cells[i].center.y = 0.0;
        g->cells[i].mass = 0.0;
    }
}

void update_grid(grid_t * g, particle_t * p, long long n_part)
{
    for(long long i = 0; i < n_part; i++)
    {
        long c = get_cell_index(g, p[i].position);

        g->cells[c].center.x += p[i].position.x * p[i].mass;
        g->cells[c].center.y += p[i].position.y * p[i].mass;
        g->cells[c].mass += p[i].mass;
    }
    
    for(long j = 0; j < g->size; j++)
    {
        if(is_cell_empty(g->cells[j]))
        {
            g->cells[j].center.x /= g->cells[j].mass;
            g->cells[j].center.x /= g->cells[j].mass;
        }
    }
}

long get_cell_index(grid_t * g, point2_t p)
{
    double xstep = GRID_MAX_X / g->ncside;
    double ystep = GRID_MAX_Y / g->ncside;

    long x = (long) (p.x / xstep);
    long y = (long) (p.y / ystep);

    return x + y * g->ncside;
}

bool is_cell_empty(cell_t c)
{
    return c.mass == 0;
}

void update_forces(grid_t * g, particle_t * p, long long n_part)
{
    for(long long i = 0; i < n_part; i++)
    {
        for(long j = 0; j < g->size; j++)
        {
            if(!is_cell_empty(g->cells[j]))
            {
                vector2_t direction;

                particle_t b;
                b.position.x = g->cells[j].center.x;
                b.position.y = g->cells[j].center.y;
                b.mass = g->cells[j].mass;

                direction.x = b.position.x - p[i].position.x;
                direction.y = b.position.y - p[i].position.y;

                double magnitude = calculate_force_magnitude(p[i], b);
                p[i].force.x += direction.x * magnitude;
                p[i].force.y += direction.y * magnitude;
            }
        }
    }
}

double calculate_force_magnitude(particle_t pA, particle_t pB)
{
    double d = calculate_distance(pA.position, pB.position);
    double f = G * ((pA.mass * pB.mass) / (d * d));

    if(d < EPSLON)
        return 0.0;

    return f;
}

double calculate_distance(point2_t pA, point2_t pB)
{
    double dx = (pB.x - pA.x);
    double dy = (pB.y - pA.y);
    return sqrt((dx * dx) + (dy * dy));
}

void update_particles(particle_t * p, long long n_part)
{
    for(long long i = 0; i < n_part; i++)
        update_particle(i, p, 1);
}

void update_particle(long long i, particle_t * p, double dt)
{
    // Get acceleration from force applied to particle
    vector2_t acceleration;
    acceleration.x = (p[i].force.x / p[i].mass);
    acceleration.y = (p[i].force.y / p[i].mass);

    // Update particle velocity
    p[i].velocity.x = p[i].velocity.x + acceleration.x * dt;
    p[i].velocity.y = p[i].velocity.y + acceleration.y * dt;

    // Update particle position
    p[i].position.x = p[i].position.x + p[i].velocity.x + 0.5 * acceleration.x * (dt * dt);
    p[i].position.y = p[i].position.y + p[i].velocity.y + 0.5 * acceleration.y * (dt * dt);

    // Update position according to grid limitations.
    p[i].position.x = (p[i].position.x >= GRID_MAX_X) ? (p[i].position.x - GRID_MAX_X) : (p[i].position.x <= GRID_MIN_X) ? (p[i].position.x + GRID_MAX_X) : p[i].position.x;
    p[i].position.y = (p[i].position.y >= GRID_MAX_Y) ? (p[i].position.y - GRID_MAX_Y) : (p[i].position.y <= GRID_MIN_Y) ? (p[i].position.y + GRID_MAX_Y) : p[i].position.y;
}

void display_particles(particle_t * p, long long n_part)
{
    for(int i=0; i < n_part; i++)
        printf("Particle #%d - p = (%f %f), v = (%f %f), f = (%f %f) | \tm = %f\n", i, p[i].position.x, p[i].position.y, p[i].velocity.x, p[i].velocity.y, p[i].force.x, p[i].force.y, p[i].mass);
}

point2_t get_center_of_mass(particle_t * p, long long n_part)
{
    double tmass = 0.0;
    point2_t center;
    center.x = 0.0;
    center.y = 0.0;

    for(int i = 0; i < n_part; i++)
    {
        center.x += p[i].mass * p[i].position.x;
        center.y += p[i].mass * p[i].position.y;
        tmass += p[i].mass;
    }

    center.x /= tmass;
    center.y /= tmass;

    return center;
}