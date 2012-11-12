############################
# Target executable and files
############################
EXE = decompose
OBJ = mesh.o tensor.o util.o decompose.o decomposeApp.o
ICPC_WARN = -wd280 -wd593 -wd869 -wd981 -wd383  -wd1419 -wd2259 -wd1572 -wd1599 -wd444 -wd1418

#############################
# paths
############################
BDIR     = ../..
ALLDIR   = field mesh tensor util vtk mp decompose solvers/piso
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
	CXXFLAGS="-O3 -Wall -fomit-frame-pointer -fstrict-aliasing -fno-exceptions -fno-rtti $(ICPC_WARN) " \
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
