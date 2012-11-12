############################
# Target executable and files
############################
EXE = mesh
OBJ = mesh.o tensor.o util.o hexMesh.o mshMesh.o meshApp.o

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
	CXX="g++ $(INC)" \
	CXXFLAGS="-O3 -msse -Wall -fomit-frame-pointer -fstrict-aliasing -fno-exceptions -fno-rtti" \
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
