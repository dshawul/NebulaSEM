############################
# Target executable and files
############################
EXE = solver
OBJ = solve.o mesh.o tensor.o util.o solver.o mp.o ke.o kw.o les.o realizableke.o rngke.o mixing_length.o field.o dg.o turbulence.o
ICPC_WARN = -wd280 -wd593 -wd869 -wd981 -wd383  -wd1419 -wd2259 -wd1572 -wd1599 -wd444 -wd1418

#############################
# paths
############################
BDIR     = ../..
SRCDIR   = $(BDIR)/src
ALLDIR   = field mesh tensor util turbulence turbulence/ke turbulence/kw turbulence/rngke turbulence/realizableke turbulence/mixing_length turbulence/les mp decompose solvers solvers/solver
EXEDIR   = $(BDIR)/bin
OBJDIR   = ../objs
OBJLIST  = $(addprefix $(OBJDIR)/, $(OBJ))
VPATH    = $(addprefix   $(SRCDIR)/, $(ALLDIR))
INC      = $(addprefix -I$(SRCDIR)/, $(ALLDIR))

######################
# Rules
######################
RM = rm -rf
DEFINES = 

default:
	$(MAKE) gcc strip

clean:
	$(RM) $(OBJLIST) $(EXEDIR)/$(EXE)

strip:
	strip $(EXEDIR)/$(EXE)

gcc:   
	$(MAKE) \
	CXX="mpic++ $(INC)" \
	CXXFLAGS="-O3 -msse -Wall -fomit-frame-pointer -fstrict-aliasing -fno-rtti" \
	LXXFLAGS="-lm" \
	all

icpc:   
	$(MAKE) \
	CXX="mpiCC $(INC)" \
	CXXFLAGS="-O3 -Wall -fomit-frame-pointer -fstrict-aliasing  -fno-rtti $(ICPC_WARN)" \
	LXXFLAGS="-lm" \
	all

all: $(EXEDIR)/$(EXE)
	
##############
# Dependencies
############## 

$(EXEDIR)/$(EXE): $(OBJLIST)
	$(CXX) $(DEFINES) -o $@ $(OBJLIST) $(LXXFLAGS)


$(OBJDIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $< -o $@
