echo "Starting job with " $2 " processors"

#generate grid
../bin/mesh $1 >grid

#decompose domain
../bin/prepare controls

#solve
mpirun --np $2  ../bin/solver controls
