#Make file for fci code
# JDWhitfield 2017
# Dartmouth College

CXX= g++
#PROFILER_FLAG=-pg 
WARNING_FLAG=-Wall -Wno-deprecated -Wfatal-errors
DEBUG_FLAG=-g  -gdwarf-2 
OPT_FLAG=  -O3 
CXXFLAGS= -std=c++11 -w


INCLUDESLibint	=-I/usr/local/include -I/usr/local/include/eigen -I/usr/local/include/libint2 -I/opt/local/include
LIBSLibint      =-L/usr/local/lib -lint2
LIBSLapack	=-L/opt/local/lib/lapack -llapack  -lblas -llapacke  
LIBS	        =-lm -lgmp -lgmpxx $(LIBSLapack) $(LIBSLibint)
INCLUDES        =$(INCLUDESLibint) 

all: weyl
	
weyl: weyl.cpp
	g++ weyl.cpp -std=c++11 -o weyl; ./weyl > output

minimal.o: minimal.cpp ci_matrix.cpp
	$(CXX) minimal.cpp $(CXXFLAGS) $(INCLUDESLibint) $(INCLUDESEig) -c -o minimal.o

minimal:minimal.o 
	$(CXX) $(CXXFLAGS) minimal.o -o minimal

fci:fci.o ci_matrix.o parser.o libint_interface.o
	$(CXX) $(CXXFLAGS) fci.o libint_interface.o ci_matrix.o parser.o $(LIBS) -o fci

parser.o: parser.cpp
	$(CXX) parser.cpp $(CXXFLAGS) $(INCLUDESLibint) $(INCLUDESEig) -c -o parser.o

libint_interface.o: libint_interface.cpp
	$(CXX) libint_interface.cpp $(INCLUDES) -std=c++11 -c -o libint_interface.o

ci_matrix.o: ci_matrix.cpp 
	$(CXX) ci_matrix.cpp $(CXXFLAGS) $(INCLUDESLibint) -c -o ci_matrix.o

fci.o: fci.cpp ci_matrix.cpp parser.cpp
	$(CXX) fci.cpp $(CXXFLAGS) $(INCLUDESLibint) $(INCLUDESEig) -c -o fci.o


