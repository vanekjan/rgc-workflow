##
## RGC skim Makefile
##
## Build from the parent project folder:
##
##   cd rgcskim
##   source setup_env.sh
##   make clean
##   make
##
## Source code:
##   src/rgcskim.cc
##
## Executable:
##   ./rgcskim
##
## The executable is kept in the parent folder because scripts/run_all_hipo.sh
## expects:
##
##   EXE=./rgcskim
##

# ------------------------------------------------------------
# HIPO libraries
# ------------------------------------------------------------

HIPOCFLAGS  := -I$(HIPO)/hipo4
HIPOLIBS    := -L$(HIPO)/lib -lhipo4

LZ4LIBS     := -L$(HIPO)/lz4/lib -llz4
LZ4INCLUDES := -I$(HIPO)/lz4/lib

# ------------------------------------------------------------
# ROOT libraries
# ------------------------------------------------------------

ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLIBS    := $(shell root-config --libs) $(shell root-config --glibs)

# ------------------------------------------------------------
# Compiler settings
# ------------------------------------------------------------

CXX      := g++
CXXFLAGS := -O2 -Wall -fPIC -std=c++17

# ------------------------------------------------------------
# Project files
# ------------------------------------------------------------

TARGET := rgcskim
SRC    := src/rgcskim.cc
OBJ    := build/rgcskim.o

# ------------------------------------------------------------
# Build rules
# ------------------------------------------------------------

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

$(OBJ): $(SRC) | build
	$(CXX) $(CXXFLAGS) $(HIPOCFLAGS) $(LZ4INCLUDES) $(ROOTCFLAGS) -c $< -o $@

build:
	mkdir -p build

clean:
	@echo 'Removing all build files'
	@rm -rf build
	@rm -f $(TARGET)
	@rm -f *~