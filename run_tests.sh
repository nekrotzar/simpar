#!/bin/bash
OS=`uname`
PROGRAM=$1
NUM_THREADS=${2:-4}
EXE="$1"

if [ ! -f "$EXE" ]; then
    echo "$0: Executable '${EXE}' not found."
	exit 1
fi

mkdir -p data
if [ -f data/system_info.txt ]; then
	rm data/system_info.txt
fi

echo "===================================="
echo 'System Info'
echo "------------------------------------"

if [[ "$OS" == "Darwin" ]]; then
	echo 'OS: macOS' | tee -a data/system_info.txt
	echo 'Processor Vendor: ' `sysctl -n machdep.cpu.vendor` | tee -a data/system_info.txt
	echo 'Processor Name/Speed: ' `sysctl -n machdep.cpu.brand_string` | tee -a data/system_info.txt
	echo 'Total Number of Cores: ' `sysctl -n machdep.cpu.core_count` | tee -a data/system_info.txt
	echo 'Memory Size (MB):' "$((`sysctl -n hw.memsize` / (1024*1024)))"  | tee -a data/system_info.txt
else
	echo 'OS: Linux' | tee -a data/system_info.txt
	echo 'Processor Vendor: ' `cat /proc/cpuinfo | grep 'vendor' | uniq` | tee -a  data/system_info.txt
	echo 'Processor Name/Speed: ' `cat /proc/cpuinfo | grep 'model name' | uniq` | tee -a data/system_info.txt
	echo 'Total Number of Cores: ' `cat /proc/cpuinfo | grep 'processor' | wc -l` | tee -a data/system_info.txt
	echo 'Memory Size (MB)' "$((`grep MemTotal /proc/meminfo | awk '{print $2}'` / 1024))" | tee -a data/system_info.txt
fi
echo "===================================="
echo

# Save system info in a file


if [[ -n "$PROGRAM" ]]; then
		
	if [[ "$PROGRAM" == "simpar" ]]; then
		label='serial'
		NUM_THREADS=1
		echo "MODE: $label"
	elif [[ "$PROGRAM" == "simpar-omp" ]]; then
		export OMP_NUM_THREADS="$NUM_THREADS"
		label='parallel'
		echo "MODE: OpenMP $label"
		echo "NUMBER OF THREADS: $NUM_THREADS"
	fi
	
	COUNT=0
	SUCCESS=0
	FAILURE=0

	# Remove file if it exists
	if [ -f data/"$label"_cpu"$NUM_THREADS"_results.csv ]; then
		rm data/"$label"_cpu"$NUM_THREADS"_results.csv
	fi
	printf "Seed,Grid size,Particles,Steps,Processors,Real time, User time,System time,CPU percentage,Label,Correctness\n" >> data/"$label"_cpu"$NUM_THREADS"_results.csv
	
	# Create output directory if not exists
	mkdir -p test/output
	for filename in test/input/*.in
	do
		echo "____________________________________"
		COUNT=$((COUNT+1))
				
		# Get arguments from file
		while read name; do
			echo Name read from file - "$name"
		done < "$filename"

		seed="$(cut -d' ' -f1 <<< $name)"
		ncside="$(cut -d' ' -f2 <<< $name)"
		npart="$(cut -d' ' -f3 <<< $name)"
		step="$(cut -d' ' -f4 <<< $name)"

		echo "TEST #$COUNT : ./"$PROGRAM" "$name""

		results=`{ (TIMEFORMAT="%R %U %S %P"; time ./"$PROGRAM" $name > test/output/test$COUNT.out 2>&1) ; } 2>&1`
		
		realt="$(cut -d' ' -f1 <<< $results)"
		usert="$(cut -d' ' -f2 <<< $results)"
		systt="$(cut -d' ' -f3 <<< $results)"
		cpupc="$(cut -d' ' -f4 <<< $results)"

		echo "Execution results"
		echo "	Elapsed time: $realt s"
		echo "	Number of CPU seconds spent in user mode: $usert s"
		echo "	Number of CPU seconds spent in system mode: $systt s"
		echo "	CPU percentage: $cpupc %"

		file1="test/validation/test$COUNT.txt"
		file2="test/output/test$COUNT.out"

		# Run the program
		# Check if the output is correct
		STATUS="$(diff "$file1" "$file2")"

		if [ -n "$STATUS" ] ; then
			test_status='FAILURE'
			echo "$test_status:"
			echo "$STATUS"
			FAILURE=$((FAILURE+1))
		else
			test_status='SUCCESS'
			echo "$test_status"
			SUCCESS=$((SUCCESS+1))
		fi

		# Store values in a csv file
		printf "$seed,$ncside,$npart,$step,$NUM_THREADS,$realt,$usert,$systt,$cpupc,$label,$test_status\n" >> data/"$label"_cpu"$NUM_THREADS"_results.csv

	done
	echo "____________________________________"
	echo "FINISHED: $SUCCESS out of $COUNT successfull tests"
	echo "____________________________________"
else
	echo "Error passing arguments: must pass program name as argument, e.g ./run.test.sh simpar"
fi