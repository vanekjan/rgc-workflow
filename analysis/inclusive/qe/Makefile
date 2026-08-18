#============================================================
# Makefile for RGC QE tensor analysis and pair diagnostics
# Location:
#   ../rgc/analysis/qe/
# Build:
#   make
# Build only one executable:
#   make qe_analysis
#   make qe_pair_scan
# Clean:
#   make clean
#============================================================

CXX = g++

QE_ANALYSIS_TARGET = qe_analysis
QE_ANALYSIS_SRC    = qe_analysis.cc

QE_PAIR_SCAN_TARGET = qe_pair_scan
QE_PAIR_SCAN_SRC    = qe_pair_scan.cc

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)

CXXFLAGS = -O2 -Wall -std=c++17

all: $(QE_ANALYSIS_TARGET) $(QE_PAIR_SCAN_TARGET)

$(QE_ANALYSIS_TARGET): $(QE_ANALYSIS_SRC)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) -o $(QE_ANALYSIS_TARGET) $(QE_ANALYSIS_SRC) $(ROOTLIBS)

$(QE_PAIR_SCAN_TARGET): $(QE_PAIR_SCAN_SRC)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) -o $(QE_PAIR_SCAN_TARGET) $(QE_PAIR_SCAN_SRC) $(ROOTLIBS)

clean:
	rm -f $(QE_ANALYSIS_TARGET) $(QE_PAIR_SCAN_TARGET)

clean_outputs:
	rm -f outputs/*.root outputs/*.png outputs/*.pdf outputs/*.csv

.PHONY: all clean clean_outputs