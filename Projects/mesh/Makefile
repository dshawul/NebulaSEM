############################
# Target executable and files
############################
EXE = mesh
OBJ = mesh.o tensor.o util.o hexMesh.o mshMesh.o meshApp.o

#############################
# paths
############################
BDIR     = ../..
SRCDIR   = $(BDIR)/src
ALLDIR   = field mesh tensor util vtk mp decompose solvers/piso
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
	$(RM) $(OBJ) $(EXE)

strip:
	strip $(EXEDIR)/$(EXE)

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
