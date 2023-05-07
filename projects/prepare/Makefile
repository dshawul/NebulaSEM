############################
# Target executable and files
############################
EXE = prepare
OBJ = mesh.o tensor.o util.o vtk.o field.o dg.o mp.o prepare.o prepareApp.o

#############################
# paths
############################
ALLDIR   = field mesh tensor util vtk mp prepare solvers/solver
METISDIR = /usr/local
INC      = -I$(METISDIR)
LINC     = -lmetis -L$(METISDIR)/lib

#############################
# include
############################
include ../../Make.inc
