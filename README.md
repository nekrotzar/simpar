# Particle Simulation
The project implements three different methods of a gravitational n-body particle simulation .  
The solvers are implemented in C language where one version is serial and the other two are parallel versions, using OpenMP and MPI.

## Usage

All the commands should be executed on the root directory of the project.  
The **-lm** flag is required when compiling the parallel code to link the math library.  

### Compile the source code
> Serial
```shell
$ gcc -lm -fopenmp -o simpar simpar.c
```
> Parallel
```shell
$ gcc -lm -fopenmp -o simpar-omp simpar-omp.c
```

### Execute the source code
All the commands can receive the following arguments:  
`seed` **required** Seed for the random number generator  
`ncside` **required** Size of the grid (number of cells on the side)  
`npart` **required** Number of particles  
`nstep` **required** number of time-steps

* **On Linux/macOS:**
```
$ ./simpar seed ncside npart nstep

# Example
$ ./simpar 1 3 10 1
```

### Set the number of threads
Use the command line / terminal to set the number of threads when running the OpenMP implementation.  

**On the command line**  
`export OMP_NUM_THREADS=2`


### Input/Ouput format
All the implementations follow the same input and output format.
	
**Output**  
The program output consists of two lines where the first line is the final position of particle #0 and the second line is center of mass of the overall space at the end of simulation. 

Output File Example: 
```
0.87 0.42
0.55 0.59
```
## Testing

Run all the tests in the [input](https://github.com/nekrotzar/simpar/test/input) directory by executing the following commands. The parallel can receive, as an optional argument, the number of desired threads *nthreads*.
> Serial
```shell
$ ./run-tests.sh simpar
```
> Parallel
```shell
$ ./run-tests.sh simpar-omp [nthreads]
```
## Contributors
* [Filipe Marques](https://github.com/Akorra)
* [Luís Fonseca](https://github.com/nekrotzar)

## License  
Licensed under MIT. See [LICENSE](LICENSE) for more information. 
