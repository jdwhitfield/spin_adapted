/******************************************************************************
 * 
 * fci.cpp
 *
 * This program loads basis functions and nuclear field information to 
 * form a Hamiltonian and compute its eigenvalues.  
 * 
 * JDWhitfield 
 * Dartmouth, 2016-2017
 *
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

#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<assert.h>
#include<iomanip>
#include<random>
#include<chrono>
#include<algorithm>
#include"parser.h"
#include"libint_interface.h"
#include"ci_matrix.h"
#include"weyl.h"


int 
main(int argc, char *argv[])
{
    bool debug=true;
    if(debug)
    {
	    std::cout.precision(10);
	    std::cout << std::scientific;
    }

    int M,N;
    Matrix h2;
    std::vector<double>  h4;

    // ****************************
    // * Parse commandline inputs *
    // ****************************
    if (argc > 1)
    {	
	    if(!strcmp(argv[1],"-h"))
	    {
		    std::cout << "fci [basis_func] [nuc_field]";
		    std::cout << "\n\tDefault files are used when called with no parameters.\n\n";
		    return 0;
	    }
    }

    const char* basis_fname=(argc>1) ? argv[1] : "basis_funcs";
    const char* nuc_fname  =(argc>2) ? argv[2] : "nuc_field";


    if(get_integrals(basis_fname, nuc_fname, M, N,  h2, h4)!=0)
    {
	    std::cout << "Problem getting integrals\n";
    }

    M=4; N=2; 
    std::cout << "M:" << M <<", N:" << N <<"\n";

    int multiplicity;

    
    if( N % 2 ) multiplicity=2;//doublet
    else        multiplicity=1;//triplet

    int nweyl=num_weyl(M/2,N,multiplicity);
    int ndets=nchoosek(M,N);
    int D=nchoosek(M,N);

    if(debug)
    {
	    std::cout << "Number of Weyl tableau: " << nweyl << std::endl;
	    std::cout << "Number of determinants: " << D << std::endl;
    }


    /*
    std::vector<int> frame;
    frame.clear();
    frame.push_back(2);
    //frame.push_back(1);

    std::vector<std::vector<int>> weyl_list;
    std::vector<int> T;
    get_init_tableau(M/2,frame,T);
    weyl_list.push_back(T);

    for(auto t : T)
	   std::cout << t << " ";
    std::cout << "\n";

    for(int i=0; i<nweyl-1; i++)
    {

	get_next_tableau(M/2,frame,T);
<<<<<<< HEAD
   	for(auto t : T)
=======
    	weyl_list.push_back(T);

   	for(t : T)
>>>>>>> cf18a1813d9393140c3115edf8290ec17de17bca
	   std::cout << t << " ";
        std::cout << "\n";
	   
    }

    std::cout << "print back Weyl list\n";
    int ctr=0;
    for(auto T:weyl_list)
    {
	    std::cout << ctr++ << ":" ;

	    for(auto t : T )
	   	std::cout << t << " ";

	    std::cout << "\n";

    }
    */

    //Eij(0,0,weyl_list[0],weyl_list[0])


    auto start_time=std::chrono::high_resolution_clock::now();
    //here is the basic CI algorithm
    double* w=ci_eigenvals(M,N,h4,h2,false);
    auto time_elapsed = std::chrono::duration_cast<std::chrono::duration<double>>
	                 (std::chrono::high_resolution_clock::now() - start_time);
    std::cout << "Computing integrals, making matrix and getting eigenvalues took " << time_elapsed.count() 
              << "s to get eigenvalues.\n";
    std::cout << "\nEigenvalues( "<< D << " ) :\n --------- \n";
    for(int i=0; i<D; i++)
        std::cout << w[i]  << "\n";


    return 0;
}


