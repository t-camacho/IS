#!/bin/bash

dir=$PWD
files=("mpi" "sequential" "openmp-CAP")

for element in ${files[@]}
do
    path="progs"/$element
    cd $path
    make clean -f makefile
	cd -
done
