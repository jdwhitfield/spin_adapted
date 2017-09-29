#Make file for fci code
# JDWhitfield 2017
# Dartmouth College

CXX= g++ -I/usr/local/include
#PROFILER_FLAG=-pg 
WARNING_FLAG=-Wall -Wno-deprecated -Wfatal-errors
DEBUG_FLAG=-g  -gdwarf-2 
OPT_FLAG=  -O3 
CXXFLAGS= -std=c++11 -w


INCLUDESLibint2	=-I/usr/local/include -I/usr/include/eigen3/Eigen -I/usr/include/eigen3/ -I/usr/local/libint/2.4.0-beta.3/include -I/opt/local/include
INCLUDESEig2	=-I/usr/include/eigen3/Eigen -I/usr/include/eigen3/ 

INCLUDESLibint	=-I/usr/local/include -I/usr/local/include/libint2 -I/opt/local/include -I/usr/local/libint/2.4.0-beta.3/include
INCLUDESEig     =-I/usr/local/include/Eigen -I/usr/local/include/ -I/usr/include/eigen/Eigen -I/usr/include/eigen3/Eigen -I/usr/include/eigen/Eigen
LIBSLibint      =-L/usr/local/lib  -L/usr/local/libint/2.4.0-beta.3/lib -lint2
LIBSLapack	=-L/opt/local/lib/lapack -llapack  -lblas -llapacke  
LIBS	        =-lm -lgmp -lgmpxx $(LIBSLapack) $(LIBSLibint)
INCLUDES        =$(INCLUDESLibint) 
EigenMac        =-I/usr/local/include/Eigen -framework accelerate -llapacke
LibintMac       =-I/usr/local/include/libint2



all: fci
	
weyl: weyl.cpp
	g++ weyl.cpp -std=c++11 -o weyl; ./weyl > output

minimal: minimal.cpp parser.o 
	$(CXX) minimal.cpp $(LibintMac)                  $(EigenMac) parser.o  -lint2 -o minimal

#	$(CXX) minimal.cpp $(CXXFLAGS) $(INCLUDESLibint) $(EigenMac) parser.o  -lint2 -o minimal

fci:fci.cpp ci_matrix.o parser.o libint_interface.o
	g++ -g -I/usr/local/include -std=c++11 fci.cpp $(LibintMac) $(EigenMac) libint_interface.o ci_matrix.o parser.o  -lint2 -o fci

test_eigen: test_eigen.cpp
	g++ -g -std=c++11 test_eigen.cpp $(EigenMac) -o test_eigen

libint_interface.o: libint_interface.cpp
	$(CXX) libint_interface.cpp $(INCLUDESLibint) -std=c++11 -c -o libint_interface.o

#fci:fci.cpp ci_matrix.o parser.o libint_interface.o
#	g++ -std=c++11 fci.cpp libint_interface.o ci_matrix.o parser.o $(EigenMac) $(INCLUDESLibint) -o fci

fci.o: fci.cpp ci_matrix.cpp parser.cpp
	$(CXX) fci.cpp $(CXXFLAGS) $(INCLUDESLibint) $(INCLUDESEig) -c -o fci.o

parser.o: parser.cpp
	$(CXX) parser.cpp $(CXXFLAGS) $(INCLUDESLibint) $(INCLUDESEig) -c -o parser.o


ci_matrix.o: ci_matrix.cpp 
	$(CXX) ci_matrix.cpp $(CXXFLAGS) $(INCLUDESLibint) -c -o ci_matrix.o



