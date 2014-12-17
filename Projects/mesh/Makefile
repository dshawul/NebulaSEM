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
ALLDIR   = mesh tensor util
EXEDIR   = $(BDIR)/bin
OBJDIR   = $(BDIR)/build
OBJLIST  = $(addprefix $(OBJDIR)/, $(OBJ))
VPATH    = $(addprefix   $(SRCDIR)/, $(ALLDIR))
INC      = $(addprefix -I$(SRCDIR)/, $(ALLDIR))
LINC     =
######################
# Compiler choice
######################
DEBUG = 0
COMP = gcc
STRIP = strip
RM = rm -rf
DEFINES =

ifeq ($(COMP),gcc)
	CXX=g++
else ifeq ($(COMP),icpc)
	CXX=icpc
endif

###########################
#  Compiler flags
###########################
CXXFLAGS = -Wall -fomit-frame-pointer -fstrict-aliasing -fno-rtti $(INC)
LXXFLAGS = -lm  $(LINC)

ifeq ($(COMP),gcc)
	CXXFLAGS += -msse
else ifeq ($(COMP),icpc)
	CXXFLAGS += -wd280 -wd593 -wd869 -wd981 -wd383  -wd1419 -wd2259 -wd1572 -wd1599 -wd444 -wd1418
endif

ifeq ($(DEBUG),3)
	ifeq ($(COMP),icpc)
		CXXFLAGS += -O2 -prof-gen
	else
		CXXFLAGS += -O2 -fprofile-generate
		LXXFLAGS += -lgcov
	endif
	STRIP=
else ifeq ($(DEBUG),2)
	ifeq ($(COMP),icpc)
		CXXFLAGS += -g -pg 
	else
		CXXFLAGS += -g -pg 
		LXXFLAGS += -g -pg
	endif
	STRIP=
else ifeq ($(DEBUG),1)
   	ifeq ($(COMP),icpc)
		CXXFLAGS += -prof-use -fast -fomit-frame-pointer
	else
		CXXFLAGS += -fprofile-use -Ofast -fomit-frame-pointer -flto
		LXXFLAGS += -lgcov
	endif
else
	ifeq ($(COMP),icpc)
		CXXFLAGS += -fast -fomit-frame-pointer
	else
		CXXFLAGS += -Ofast -fomit-frame-pointer -flto
	endif
endif


######################
# Rules
######################

default:
	$(MAKE) all strip

clean:
	$(RM) $(OBJLIST) $(EXEDIR)/$(EXE)

strip:
	$(STRIP) $(EXEDIR)/$(EXE)

all: $(EXEDIR)/$(EXE)
	
help:
	@echo ""
	@echo "1. make [DEBUG=n] [COMP=c]"
	@echo ""
	@echo "  n ="
	@echo "	0: Compile optimized binary (-03)"
	@echo "	1: Compile with profile guided optimization (PGO)"
	@echo "	2: Compile for deugging (default)"
	@echo "	3: Prepare for PGO"
	@echo ""
	@echo "  c ="
	@echo "	gcc    :  g++ compiler"
	@echo "	icpc   :  intel compiler"
	@echo ""
	@echo "2. make clean - removes all files but source code"
	@echo "3. make strip - strips executable of debugging/profiling data"
	@echo ""

##############
# Dependencies
############## 

$(EXEDIR)/$(EXE): $(OBJLIST)
	$(CXX) $(DEFINES) -o $@ $(OBJLIST) $(LXXFLAGS)


$(OBJDIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $< -o $@
