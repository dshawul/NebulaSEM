#!/bin/sh

#set case file
case=examples/backface/ 
grid=backface

#run solver
function run() {
	echo "Starting job with " $2 " processors"

	#generate grid
	../bin/mesh $1 >grid

	#solve
	mpirun --np $2  ../bin/solver controls
}

#prepare directory
rundir=run$1
rm -rf $rundir
cp -r $case $rundir
cd $rundir

#modify controls
if [ ! -z $2 ]; then
    sed -i "s/3 {2 1 1}/3 {"$2" "$3" 1}/g" ./controls
fi

#run
run $grid $1

#vtk
../bin/prepare ./controls -vtk

cd ..
