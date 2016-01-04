############################
# Target executable and files
############################
EXE = solver
OBJ = solve.o mesh.o tensor.o util.o solver.o mp.o ke.o kw.o les.o realizableke.o rngke.o mixing_length.o field.o dg.o turbulence.o

#############################
# paths
############################
ALLDIR   = field mesh tensor util turbulence turbulence/ke turbulence/kw turbulence/rngke turbulence/realizableke turbulence/mixing_length turbulence/les mp decompose solvers solvers/solver
METISDIR = /usr/local
INC      = -I$(METISDIR)
LINC     = -lmetis -L$(METISDIR)/lib

#############################
# include
############################
include ../../Make.inc
