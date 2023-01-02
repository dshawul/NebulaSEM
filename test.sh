#!/bin/bash

set -e

#set case file
case=examples/cavity/
grid=cavity

#run solver
run() {
	echo "Starting job with " $2 " processors"

	#generate grid
	../bin/mesh $1 -o grid_0.bin

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
if [ $1 -gt 1 ]; then
    t=1
    file=(./grid0/U$t.*)
    while [ -f "$file" ]; do
        mpirun -n $1 ../bin/prepare ./controls -merge -start $t
        t=$((t+1))
        file=(./grid0/U$t.*)
    done
fi

#vtk
t=0
file=(./U$t.*)
while [ -f "$file" ]; do
    ../bin/prepare ./controls -vtk -start $t
    t=$((t+1))
    file=(./U$t.*)
done

cd ..
