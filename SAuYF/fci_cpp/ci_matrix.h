/******************************************************************************
 * 
 * ci_matrix.h
 *
 * This file has functions to form a Hamiltonian and compute its eigenvalues.  
 * 
 * JDWhitfield 
 * Dartmouth, 2016-2017
 * 
 *
 * Parts of the code are borrowed from the open-source projects 
 * Libints and PyQuante under the GNU General Public License.
 *
 ***
 *
 *      This program is free software: you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation, either version 2 of the License, or
 *      (at your option) any later version.
 * 
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 * 
 *      You should have received a copy of the GNU General Public License
 *      along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 *****************************************************************************/


// //start of header guard
#ifndef CI_LIB
#define CI_LIB

#include<cmath>
#include<stdio.h>
#include<iostream>
#include<string>
#include<random>
#include<vector>
// Eigen matrix algebra library
//#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#include<Eigen> 
//<eigen3/Eigen/Dense>
//#include<Eigen/Dense> //<eigen3/Eigen/Dense>
// Libint Gaussian integrals library
#include <libint2.hpp>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        Matrix;  // import dense, dynamically sized Matrix type from Eigen;
                 // this is a matrix with row-major storage (http://en.wikipedia.org/wiki/Row-major_order)
                 // to meet the layout of the integrals returned by the Libint integral library
                 
using libint2::Shell;
using std::vector;
// **********************************/
//prototypes

double*
ci_matrix_packed_upper(int M, int N, std::vector<double> MOInts, Matrix h2);

Matrix
ci_matrix(int M, int N, std::vector<double> MOInts, Matrix h2);

double*
ci_eigenvals(int M, int N, std::vector<double> MOInts, Matrix h2, int DEBUG=0);

int 
spin2spatial(int);

//Hall-Robinson formula
int
num_weyl(int M, int N, int ms);


//based on my Matlab implementations
int 
nchoosek(double n, double k);

std::vector<int>
CboI(int idx, int M, int N);


int
idx(int n,int m, int M);

std::vector<size_t> map_shell_to_basis_function(const std::vector<Shell>& );


int 
term4(int, int, int, int);

std::vector<double> 
compute_2body_ints_old(const std::vector<libint2::Shell>&);


Matrix 
compute_1body_ints_old(const std::vector<libint2::Shell>& shells,
                   libint2::Operator obtype,
                   const std::vector<libint2::Atom>& atoms = std::vector<libint2::Atom>());



//makes a Harr random unitary based on the shit Eigen SVD decomp
Matrix
randU(int M, std::default_random_engine& generator);

//based on PyQuante implementation 

std::vector<double>
transform4(const int M, const Matrix& C, const vector<double>& AOInts, int debug=0);

/* taken from Libint implementation of Hartree-Fock
 * Copyright (C) 2004-2014 Edward F. Valeev
 *
 * Used with rights derived from the GNU General Public License
 */
void
get_libints(const std::vector<libint2::Shell>& shells, const std::vector<libint2::Atom>& atoms,
            Matrix &S, Matrix &hV, Matrix &hT, vector<double>& AOInts);


double
zero_excitation(const vector<int> detL, const vector<int> detR,
                const vector<int> virt,  const vector<int> occ,
                std::vector<double> MOInts, Matrix h2,
                bool spatial_orbs=true);

double
single_excitation(const vector<int> detL,  const vector<int> detR,
                  const vector<int> virt,  const vector<int> occ,
                  std::vector<double> MOInts, Matrix h2,
                  bool spatial_orbs=true);

double
double_excitation(const vector<int> detL,  const vector<int> detR,
                  const vector<int> virt,  const vector<int> occ,
                  std::vector<double> MOInts, Matrix h2,
                  bool spatial_orbs=true);

int
SlaterDiffs(const vector<int>&, const vector<int>&, vector<int>& , vector<int>& , int&);


#endif
