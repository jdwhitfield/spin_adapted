#Make file for fci code
# JDWhitfield 2017
# Dartmouth College

CXX= g++ 
#PROFILER_FLAG=-pg 
WARNING_FLAG=-Wall -Wno-deprecated -Wfatal-errors
DEBUG_FLAG=-g  -gdwarf-2 
OPT_FLAG=  -O3 
CXXFLAGS= -std=c++11 -w



INCLUDESLibint	=-I/usr/local/include -I/usr/local/libint/2.4.0-beta.3/include -I/usr/local/include/libint2
LIBSLibint      =-L/usr/local/lib  -L/usr/local/libint/2.4.0-beta.3/lib -lint2 
INCLUDESEig     =-I/usr/local/include/eigen/Eigen -I/usr/local/include/eigen -I/usr/local/include/ -I/usr/include/eigen/Eigen -I/usr/include/eigen3/Eigen -I/usr/include/Eigen/ -I/usr/include/Eigen/Eigen
LIBSLapack	=-llapack   -llapacke  
LIBSMac         =-framework accelerate
LIBS	        =-lm -lgmp -lgmpxx $(LIBSLapack) $(LIBSLibint) 
INCLUDES        =$(INCLUDESLibint)  $(INCLUDESEig)

all: minimal	

test_wigner_projectors: test_wigner_projectors.cpp weyl.o young.o
	$(CXX) test_wigner_projectors.cpp $(CXXFLAGS) $(INCLUDESLibint) $(INCLUDESEig) weyl.o young.o $(LIBS) -o test_wigner_projectors

test_irrep_orthogonality: test_irrep_orthogonality.cpp weyl.o young.o
	$(CXX) test_irrep_orthogonality.cpp $(CXXFLAGS) $(INCLUDESLibint) $(INCLUDESEig) weyl.o young.o $(LIBS) -o test_irrep_orthogonality

young.o: young.cpp 
	$(CXX) $(CXXFLAGS) $(DEBUG_FLAG) $(INCLUDESEig) $(INCLUDESLibint) young.cpp -c -o young.o

young: young.cpp weyl.o
	$(CXX) $(CXXFLAGS) $(DEBUG_FLAG) $(INCLUDESEig) $(INCLUDESLibint) young.cpp weyl.o -o young; ./young > out; beep

weyl: weyl.cpp
	g++ weyl.cpp -std=c++11 -o weyl; ./weyl > output

weyl.o: weyl.cpp
	g++ weyl.cpp -std=c++11 -c -o weyl.o

minimal: minimal.cpp weyl.o young.o
	$(CXX) minimal.cpp $(CXXFLAGS) $(INCLUDESLibint) $(INCLUDESEig) weyl.o young.o $(LIBS) -o minimal; ./minimal > out

#minimal.o: minimal.cpp ci_matrix.cpp
#	$(CXX) minimal.cpp $(CXXFLAGS) $(INCLUDESLibint) $(INCLUDESEig) -c -o minimal.o

test_eigen: test_eigen.cpp
	g++ -std=c++11 test_eigen.cpp -I/usr/local/include/eigen/Eigen -llapack -llapacke -o test_eigen

fci:fci.o ci_matrix.o parser.o libint_interface.o weyl.o
	$(CXX) $(CXXFLAGS) fci.o libint_interface.o weyl.o ci_matrix.o parser.o $(LIBS) -o fci

libint_interface.o: libint_interface.cpp
	$(CXX) -g libint_interface.cpp $(INCLUDES)  -std=c++11 -c -o libint_interface.o

ci_matrix.o: ci_matrix.cpp 
	$(CXX) ci_matrix.cpp $(CXXFLAGS) $(INCLUDESLibint) $(INCLUDESEig) -c -o ci_matrix.o

fci.o: fci.cpp ci_matrix.cpp parser.cpp
	$(CXX) -g fci.cpp $(CXXFLAGS) $(INCLUDESLibint) $(INCLUDESEig) -c -o fci.o

parser.o: parser.cpp
	$(CXX) parser.cpp $(CXXFLAGS) $(INCLUDESLibint) $(INCLUDESEig) -c -o parser.o



