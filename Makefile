OS := $(shell uname)
ifeq ($(OS), Darwin)
	CC=gcc-8
else
	CC=gcc
endif

all: simpar simpar-omp

simpar: simpar.c
		$(CC) -fopenmp -o simpar simpar.c

simpar-omp: simpar-omp.c
		$(CC) -fopenmp -o simpar-omp simpar-omp.c

clean:
		-rm -f test/output/*.out
		-rm -f *.o
		-rm -f simpar
		-rm -f simpar-omp