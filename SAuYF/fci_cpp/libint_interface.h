/******************************************************************************
 * 
 * libint_interface.h
 *
 * This file has functions to interface with libints
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


//
// //start of header guard
#ifndef LIBINT_INTERFACE
#define LIBINT_INTERFACE

#include <libint2.hpp>
#include "parser.h"
#include "ci_matrix.h"

int
get_integrals(const char* basis_fname, const char* nuc_fname, int& M, int& N, Matrix& h2, std::vector<double>& h4);

void
get_libints(const libint2::BasisSet& basisset, const std::vector<libint2::Atom>& atoms,
            Matrix &S, Matrix &hV, Matrix &hT, vector<double>& AOInts);

Matrix
compute_1body_ints(libint2::Operator, 
		   const libint2::BasisSet, const std::vector<libint2::Atom>& atoms=std::vector<libint2::Atom>());

std::vector<double>
compute_2body_ints(const libint2::BasisSet bfs);

#endif
// //end of header guard
