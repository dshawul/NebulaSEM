#!/bin/bash

set -e

#display help
usage() {
    echo "Usage: $0 [options]" >&2
    echo
    echo "   -n,--np       Number of processors to use."
    echo "   -c,--case     Path to grid name under a test case directory."
    echo "   -b,--bin-path Path to biary files mesh, prepare and solvers."
    echo "   -h,--help     Display this help message."
    echo
    exit 0
}

#default test case
case=examples/cavity/cavity
nprocs=1
bin_path=../bin/

#process command line args
while :; do
    case $1 in
        --help|-h) usage; exit 0;;
        --np|-n) shift; nprocs=$1;; 
        --case|-c) shift; case=$1;; 
        --bin-path|-b) shift; bin_path=$1;; 
        *) break
    esac
    shift
done

#get base and dir names
grid=$(basename $case)
case=$(dirname $case)

#run solver
run() {
	echo "Starting job with " $2 " processors"

	#generate grid
	${bin_path}mesh $1 -o grid_0.bin

	#solve
	solver=$(grep -v "^#" ./controls | grep solver | awk '{print $2}')
	mpirun -n $2  ${bin_path}${solver} ./controls
}

#prepare directory
rundir=run$nprocs
rm -rf $rundir
cp -r $case $rundir
cd $rundir

#run
run $grid $nprocs

#find the control field
fields=$(cat controls | grep fields)
fields="${fields//\{/ }"
fields="${fields//\}/ }"
field=$(echo $fields | awk '{print $3}')

#merge
if [ $nprocs -gt 1 ]; then
    t=1
    file=(./grid0/$field$t.*)
    while [ -f "$file" ]; do
        mpirun -n $nprocs ${bin_path}prepare ./controls -merge -start $t
        t=$((t+1))
        file=(./grid0/$field$t.*)
    done
fi

#vtk
t=0
file=(./$field$t.*)
while [ -f "$file" ]; do
    ${bin_path}prepare ./controls -vtk -start $t
    t=$((t+1))
    file=(./$field$t.*)
done

cd ..
