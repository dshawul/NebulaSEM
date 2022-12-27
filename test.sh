#!/bin/sh

set -ex

#set case file
case=examples/cavity/
grid=cavity

#run solver
function run() {
	echo "Starting job with " $2 " processors"

	#generate grid
	../bin/mesh $1 >grid_0

	#solve
	mpirun -n $2  ../bin/solver ./controls
}

#prepare directory
rundir=run$1
rm -rf $rundir
cp -r $case $rundir
cd $rundir

#run
run $grid $1

#merge
if [[ $1 -gt 1 ]]; then
    t=1
    while [ -f "./grid0/U$t" ]; do
	    mpirun -n $1 ../bin/prepare ./controls -merge -start $t
        t=$((t+1))
    done
fi

#vtk
../bin/prepare ./controls -vtk

cd ..
