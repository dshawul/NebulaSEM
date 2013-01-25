############################
# Target executable and files
############################
EXE = solver
OBJ = solve.o mesh.o tensor.o util.o solver.o mp.o ke.o kw.o les.o realizableke.o rngke.o mixing_length.o field.o 
ICPC_WARN = -wd280 -wd593 -wd869 -wd981 -wd383  -wd1419 -wd2259 -wd1572 -wd1599 -wd444 -wd1418

#############################
# paths
############################
BDIR     = ../..
ALLDIR   = field mesh tensor util turbulence turbulence/ke turbulence/kw turbulence/rngke turbulence/realizableke turbulence/mixing_length turbulence/les mp decompose solvers solvers/solver
EXEDIR   = $(BDIR)/Projects/Release
OBJDIR   = ./Release
OBJLIST  = $(addprefix $(OBJDIR)/, $(OBJ))
VPATH    = $(addprefix   $(BDIR)/, $(ALLDIR))
INC      = $(addprefix -I$(BDIR)/, $(ALLDIR))

######################
# Rules
######################
RM = rm -rf
DEFINES = 

default:
	$(MAKE) gcc
clean:
	$(RM) $(OBJ) $(EXE)

strip:
	strip $(EXE)

gcc:   
	$(MAKE) \
	CXX="mpiCC $(INC)" \
	CXXFLAGS="-O3 -msse -Wall -fomit-frame-pointer -fstrict-aliasing -fno-exceptions -fno-rtti" \
	LXXFLAGS="-lm" \
	all

icpc:   
	$(MAKE) \
	CXX="mpiCC $(INC)" \
	CXXFLAGS="-O3 -Wall -fomit-frame-pointer -fstrict-aliasing -fno-exceptions -fno-rtti $(ICPC_WARN)" \
	LXXFLAGS="-lm" \
	all

all: $(EXEDIR)/$(EXE)
	
##############
# Dependencies
############## 

$(EXEDIR)/$(EXE): $(OBJ)
	$(CXX) $(LXXFLAGS) $(DEFINES) -o $@ $(OBJLIST)


%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $< -o $(OBJDIR)/$@
