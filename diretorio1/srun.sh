#!/bin/bash

TIMEFORMAT=%R
dir=$PWD
EXECUTIONS=11
files=("mpi" "sequential")

mkdir output

path="progs/sequential"
cd $path
for ((i=1;i<=${EXECUTIONS};i++)); do
	echo -e "\n--> sequential $i/${EXECUTIONS}\n\n"
	(time ./bin/is.intel --class huge) 2>> ../../output/tempo_sequential.txt
done
cd -

for element in ${files[@]}; do
	path="progs"/$element
	cd $path
	for t in 2 4 8 16; do
		for ((i=1;i<=${EXECUTIONS};i++)); do
           	echo -e "\n--> $element ($t flows) $i/${EXECUTIONS}\n"
        	(time mpirun -n $t ./bin/is.intel --class huge) 2>> ../../output/tempo_${element}_${t}_threads.txt
		done
    done
	cd -
done
