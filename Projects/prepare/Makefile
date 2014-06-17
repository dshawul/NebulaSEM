############################
# Target executable and files
############################
EXE = prepare
OBJ = mesh.o tensor.o util.o vtk.o field.o mp.o prepare.o prepareApp.o
ICPC_WARN = -wd280 -wd593 -wd869 -wd981 -wd383  -wd1419 -wd2259 -wd1572 -wd1599 -wd444 -wd1418

#############################
# paths
############################
BDIR     = ../..
SRCDIR   = $(BDIR)/src
ALLDIR   = field mesh tensor util vtk mp prepare solvers/piso
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

$(EXEDIR)/$(EXE): $(OBJLIST)
	$(CXX) $(LXXFLAGS) $(DEFINES) -o $@ $(OBJLIST)


$(OBJDIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $< -o $@
