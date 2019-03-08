#!/bin/bash
PROGRAM=$1
NUM_THREADS=${2:-4}

if [[ -n "$PROGRAM" ]]; then
	#Set number of threads
	
	printf "\n*******************************************\n"
	printf "\Run all particle simulations.\n\n"

	if [[ "$PROGRAM" == "simpar" ]]; then
		echo "MODE: serial"
	elif [[ "$PROGRAM" == "simpar-omp" ]]; then
		export OMP_NUM_THREADS="$NUM_THREADS"

		echo "MODE: openMP parallel"
		echo "NUMBER OF THREADS: $NUM_THREADS"
	fi
	printf "*******************************************\n\n"

    echo "TEST: seed = 1, ncside = 3, n_particle = 10, nstep = 1"
    time ./"$PROGRAM" 1 3 10 1
	echo "=========================================="

    echo "TEST: seed = 1, ncside = 3, n_particle = 1000000, nstep = 20"
    time ./"$PROGRAM" 1 3 1000000 20
	echo "=========================================="

    echo "TEST: seed = 1, ncside = 10, n_particle = 2000000, nstep = 10"
    time ./"$PROGRAM" 1 10 2000000 10
	echo "=========================================="

    echo "TEST: seed = 1, ncside = 30, n_particle = 20000000, nstep = 10"
    time ./"$PROGRAM" 1 30 20000000 10
	echo "=========================================="


else
	echo "Error passing arguments: must pass program name as argument, e.g sudoku-serial ."
fi

