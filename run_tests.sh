#!/bin/bash
OS=`uname`
PROGRAM=$1
NUM_THREADS=${2:-4}
EXE="$1"

if [ ! -f "$EXE" ]; then
    echo "$0: Executable '${EXE}' not found."
	exit 1
fi

echo "===================================="
echo 'System Info'
echo "------------------------------------"

if [[ "$OS" == "Darwin" ]]; then
	echo 'OS: macOS'
	echo 'Processor Vendor: ' `sysctl -n machdep.cpu.vendor`
	echo 'Processor Name/Speed: ' `sysctl -n machdep.cpu.brand_string`
	echo 'Total Number of Cores: ' `sysctl -n machdep.cpu.core_count`
	echo 'Total Number of Threads: ' `sysctl -n machdep.cpu.thread_count`
	echo 'Memory Size (MB):' "$((`sysctl -n hw.memsize` / (1024*1024)))"
else
	echo 'OS: Linux'
	echo 'Processor Vendor: ' `cat /proc/cpuinfo | grep 'vendor' | uniq`
	echo 'Processor Name/Speed: ' `cat /proc/cpuinfo | grep 'model name' | uniq`
	echo 'Total Number of Cores: ' `cat /proc/cpuinfo | grep 'processor' | wc -l`
	echo 'Memory Size (MB)' "$((`grep MemTotal /proc/meminfo | awk '{print $2}'` / 1024))"
fi
echo "===================================="
echo

if [[ -n "$PROGRAM" ]]; then
	
	# Create output directory if not exists
	mkdir -p test/output
	
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
		echo "____________________________________"
		COUNT=$((COUNT+1))
				
		# Get arguments from file
		while read name; do
			echo "Name read from file - $name"
		done < "$filename"

		echo "TEST #$COUNT : ./"$PROGRAM" "$name""
		echo "Execution times : " `{ time ./"$PROGRAM" $name > test/output/test$COUNT.out 2>&1 ; } 2>&1`

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
			echo "SUCCESS"
			SUCCESS=$((SUCCESS+1))
		fi
	done
	echo "____________________________________"
	echo "FINISHED: $SUCCESS out of $COUNT successfull tests"
	echo "____________________________________"
else
	echo "Error passing arguments: must pass program name as argument, e.g ./run.test.sh simpar"
fi