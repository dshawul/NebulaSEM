#!/bin/sh

#set case file
case=examples/cavity/ 
grid=cavity

#run solver
function run() {
	echo "Starting job with " $2 " processors"

	#generate grid
	../bin/mesh $1 >grid_0

	#solve
	mpirun --np $2  ../bin/solver controls
}

#prepare directory
rm -rf run$1
cp -r $case run$1

#change controls
cd run$1
sed -i '' -e "s/3 {2 1 1}/3 {"$2" "$3" 1}/g" controls
run $grid $1
cd ..
