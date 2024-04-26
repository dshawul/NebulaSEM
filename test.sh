#!/bin/bash

platform=$(uname)

set -e

#display help
usage() {
    echo "Usage: $0 [options]" >&2
    echo
    echo "   -n,--np       Number of MPI processes to launch."
    echo "   -c,--case     Path to grid file name that is under a test case directory."
    echo "   -b,--bin-path Path to binaries: mesh, prepare and solvers."
    echo "   -s,--steps    Number of time steps, which overwrites the one in control file."
    echo "   -h,--help     Display this help message."
    echo
    exit 0
}

#default test case
case=examples/cavity/cavity
nprocs=1
bin_path=../bin/
steps=0

#process command line args
while :; do
    case $1 in
        --help|-h) usage; exit 0;;
        --np|-n) shift; nprocs=$1;; 
        --case|-c) shift; case=$1;; 
        --bin-path|-b) shift; bin_path=$1;; 
        --steps|-s) shift; steps=$1;;
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
    format=$(grep -v "^#" ./controls | grep write_format | awk '{print $2}')
    if [ "$format" == "TEXT" ]; then
        ${bin_path}/mesh $1 -o grid_0.txt
    else
        ${bin_path}/mesh $1 -o grid_0.bin
    fi

    bind_numa="-bind-to numa"
    if [ "$platform" == "Darwin" ]; then
       bind_numa=
    fi

    #solve
    solver=$(grep -v "^#" ./controls | grep solver | awk '{print $2}')
    mpirun -n $2 ${bind_numa} ${bin_path}/${solver} ./controls | tee log.txt
}

#check case/grid exists
if ! [ -f ${case}/${grid} ]; then
  echo "Grid file does not exist: "${case}/${grid}
  exit 0
fi

#prepare directory
rundir="run$nproces-${case//\//-}"
rm -rf $rundir
cp -r $case $rundir
cd $rundir

if [ $steps -ne 0 ]; then
    sed -i "/^[[:space:]]*end_step/c\    end_step      $steps" ./controls
    sed -i "/print_time/d" ./controls
fi

#run
run $grid $nprocs

#find the control field
fields=$(grep fields ./controls)
fields="${fields//\{/ }"
fields="${fields//\}/ }"
field=$(echo $fields | awk '{print $3}')

#merge
if [ $nprocs -gt 1 ]; then
    t=1
    file=(./grid0/$field$t.*)
    while [ -f "$file" ]; do
        mpirun -n $nprocs ${bin_path}/prepare ./controls -merge -start $t
        t=$((t+1))
        file=(./grid0/$field$t.*)
    done
fi

#vtk
t=0
file=(./$field$t.*)
while [ -f "$file" ]; do
    mpirun -n 1 ${bin_path}/prepare ./controls -vtk -start $t
    t=$((t+1))
    file=(./$field$t.*)
done

cd ..
