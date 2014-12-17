############################
# Target executable and files
############################
EXE = prepare
OBJ = mesh.o tensor.o util.o vtk.o field.o dg.o mp.o prepare.o prepareApp.o

#############################
# paths
############################
METISDIR = /usr/local
ALLDIR   = field mesh tensor util vtk mp prepare solvers/piso $(METISDIR)/include
LINC     = -lmetis -L$(METISDIR)/lib

#############################
# include
############################
include ../../Make.inc