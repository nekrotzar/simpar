#include "simpar.h"

void usg_err()
{
    printf("\t[-] usage : ./simpar <seed> <ncside> <n_par> <n_step>");
    printf("\t\t[-] int <seed> : seed for random number generation.\n");
    printf("\t\t[-] int <ncside> :  size of the grid (number of cells on the side.\n");
    printf("\t\t[-] int <n_par> :  number of particles\n");
    printf("\t\t[-] int <n_par> :  number of time-steps\n");
    exit(1);
}

long long val_l(const char* arg)
{
    char *endptr;
    long long x = strtol(arg, &endptr, 10);
    if (endptr == arg) {
        printf("[-] ERROR: Invalid number: %s\n", arg);
        return 0;
    } else if (*endptr) {
        printf("[-] ERROR: Trailing characters after number: %s\n", arg);
    } else if (x <= 0) {
        printf("[-] ERROR: Number must be positive: %llu\n", x);
        return 0;
    }
    return x;
}

cell_t** init_grid(const long ncside)
{
    cell_t** grid = (cell_t**) calloc(ncside, sizeof(cell_t*));
    for(long c=0; c<ncside; c++)
    {
        grid[c] = (cell_t*) calloc(ncside, sizeof(cell_t));
        
        if(grid[c]==NULL) 
            exit(0);
    }
    return grid;
}

void free_grid(cell_t** g, long ncside)
{
    for(long c=0; c<ncside; c++)
    {
        free(g[c]);
    }
    free(g);
}

void init_particles(long seed, long ncside, long long n_part, particle_t *par)
{
    long long i;
    srandom(seed);

    for(i=0; i < n_part; i++)
    {
        par[i].x = RND0_1;
        par[i].y = RND0_1;
        par[i].vx = RND0_1 / ncside / 10.0;
        par[i].vy = RND0_1 / ncside / 10.0;
        par[i].m = RND0_1 * ncside / (G * 1e6 * n_part);
    }
}

void init_env(cell_t** g, long ncside, particle_t* p, long long n_par)
{
    for(long long i=0; i<n_par; i++)
    {
        p[i].cx = (long) p[i].x * ncside;
        p[i].cy = (long) p[i].y * ncside;
        g[p[i].cx][p[i].cy].M += p[i].m;
        g[p[i].cx][p[i].cy].x += p[i].m * p[i].x;
        g[p[i].cx][p[i].cy].y += p[i].m * p[i].y;
    }
}

void accellerate_p(double* ax, double* ay, const cell_t* c, double m, double x, double y)
{
    if((c->M) == 0.0)
        return;

    //double dirx = 1.0, diry = 1.0,
    double magnitude;
    double dx = ((c->x)/(c->M)) - x;
    double dy = ((c->y)/(c->M)) - y;

    double d_2 = (dx*dx)+(dy*dy);
    
    if(sqrt(d_2) < EPSLON)
        return;

    //if(dx<0.0){ dirx = -1.0; }else if(dx == 0.0){ dirx = 0.0; }
    //if(dy<0.0){ diry = -1.0; }else if(dy == 0.0){ diry = 0.0; }
    magnitude = (((c->M)*G)/d_2);
    *ax += dx * magnitude;
    *ay += dy * magnitude;
}

void update_particles(cell_t** x, long ncside, particle_t* par, long long n_par, long n_step, long step)
{
    for(long long i=0; i<n_par; i++)
    {
        double m = par[i].m;
        double px = par[i].x;
        double py = par[i].y;
        long cx = (long) px * ncside, nx;
        long cy = (long) py * ncside, ny;
        long ux = cx+1, uy = cy+1, lx = cx-1, ly = cy-1;
        
        if(ux >= ncside)
            ux = 0;
        else if(lx < 0)
            lx = ncside-1;
        
        if(uy >= ncside)
            uy = 0;
        else if(ly < 0)
            ly = ncside-1;

        double ax = 0.0;
        double ay = 0.0;

        accellerate_p(&ax, &ay, &(x[cx][cy]), m, px, py); // current cell
        accellerate_p(&ax, &ay, &(x[ux][cy]), m, px, py); // right cell
        accellerate_p(&ax, &ay, &(x[lx][cy]), m, px, py); // left cell
        //upper adjacents
        accellerate_p(&ax, &ay, &(x[cx][uy]), m, px, py); // upper cell
        accellerate_p(&ax, &ay, &(x[lx][uy]), m, px, py); // upper left cell
        accellerate_p(&ax, &ay, &(x[ux][uy]), m, px, py); // upper right cell
        //lower adjacents
        accellerate_p(&ax, &ay, &(x[cx][ly]), m, px, py); // lower cell
        accellerate_p(&ax, &ay, &(x[lx][ly]), m, px, py); // lower left cell
        accellerate_p(&ax, &ay, &(x[ux][ly]), m, px, py); // lower right cell

        //update velocity
        par[i].vx += ax;
        par[i].vy += ay;

        //update position
        par[i].x += par[i].vx + ax*0.5;
        
        while(par[i].x >= 1.0)
            par[i].x -= 1.0;
        while(par[i].x < 0.0)
            par[i].x += 1.0;

        par[i].y += par[i].vy + ay*0.5;
        
        while(par[i].y >= 1.0)
            par[i].y -= 1.0;
        while(par[i].y < 0.0)
            par[i].y += 1.0;

        //update cells if cell changed maybe outside loop?
        nx = (long) par[i].x*ncside;
        ny = (long) par[i].y*ncside;
        if(cx-nx || cy-ny)
        {
            if(cx-nx) par[i].cx = nx;
            if(cy-ny) par[i].cy = ny;

            x[cx][cy].M -= m;
            x[cx][cy].x -= m * px;
            x[cx][cy].y -= m * py;

            x[nx][ny].M += m;
            x[nx][ny].x += m * par[i].x;
            x[nx][ny].y += m * par[i].y;
        }

        if(n_step-1-step == 0)
        {
            t_mass += par[i].m;
            t_cx += par[i].m * par[i].x;
            t_cy += par[i].m * par[i].y;
        }
    }
}

int main(int argc, const char * argv[]) {
    if(argc != 5)
    {
        printf("[-] ERROR: Invalid number of arguments... Expected 4 but got %d\n", argc-1);
        usg_err();
    }

    const long seed = (long) val_l(argv[1]);
    const long ncside = (long) val_l(argv[2]);
    const long long n_par = val_l(argv[3]);
    const long n_step = (long) val_l(argv[4]);
    
    if(!(seed*ncside*n_par*n_step))
        usg_err();

    clock_t start_t, end_t;
    double elapsed_t;
    start_t = clock();

    particle_t* par = (particle_t*)calloc(n_par, sizeof(particle_t));
    init_particles(seed, ncside, n_par, par);

    cell_t** grid = init_grid(ncside);
    if(grid==NULL || par == NULL) 
        exit(0);

    init_env(grid, ncside, par, n_par);

    for(long step = 0; step < n_step; step++)
    {
        update_particles(grid, ncside, par, n_par, n_step, step);
    }

    t_cx /= t_mass;
    t_cy /= t_mass;
    printf("%.2f %.2f\n", par[0].x, par[0].y);
    printf("%.2f %.2f\n", t_cx, t_cy);

    end_t = clock();
    elapsed_t = ((double) (end_t - start_t)) / CLOCKS_PER_SEC;
    //printf("%f (s)\n", elapsed_t);

    free(par);
    free_grid(grid, ncside);
    return 0;
}
