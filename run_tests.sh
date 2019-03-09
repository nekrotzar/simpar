#!/bin/bash
PROGRAM=$1
NUM_THREADS=${2:-4}

if [[ -n "$PROGRAM" ]]; then
	#Set number of threads
	
	if [[ "$PROGRAM" == "simpar" ]]; then
		echo "MODE: serial"
	elif [[ "$PROGRAM" == "simpar-omp" ]]; then
		export OMP_NUM_THREADS="$NUM_THREADS"
		echo "MODE: OpenMP parallel"
		echo "NUMBER OF THREADS: $NUM_THREADS"
	fi
	
	COUNT=0
	SUCCESS=0
	FAILURE=0

	for filename in test/input/*.in
	do
		echo "__________________________________________"
		COUNT=$((COUNT+1))
				
		# Get arguments from file
		while read name; do
			echo "Name read from file - $name"
		done < "$filename"

		echo "TEST #$COUNT : ./"$PROGRAM" "$name""
		
		time ./"$PROGRAM" $name > test/output/test$COUNT.out

		file1="test/validation/test$COUNT.txt"
		file2="test/output/test$COUNT.out"

		# Run the program
		# Check if the output is correct
		STATUS="$(diff "$file1" "$file2")"

		if [ -n "$STATUS" ] ; then
			echo "FAILURE :"
			echo "$STATUS"
			FAILURE=$((FAILURE+1))
		else
			echo
			echo "SUCCESS"
			SUCCESS=$((SUCCESS+1))
		fi
	done
	echo "__________________________________________"
	echo "FINISHED: $SUCCESS out of $COUNT successfull tests"
	echo "------------------------------------------"
else
	echo "Error passing arguments: must pass program name as argument, e.g ./run.test.sh simpar"
fi