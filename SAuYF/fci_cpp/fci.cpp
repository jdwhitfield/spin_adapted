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
//#include"ci_matrix.h"
#include"weyl.h"


int 
main(int argc, char *argv[])
{

    int M,N;
    Matrix h2;
    std::vector<double>  h4;
    bool debug=true;

    if(debug)
    {
	    std::cout.precision(10);
	    std::cout << std::scientific;
    }

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

    std::cout << "M:" << M <<", N:" << N <<"\n";

    int ndets=nchoosek(M,N);
    int multiplicity;

    
    if( N % 2 ) multiplicity=2;//doublet
    else        multiplicity=1;//triplet

    int nweyl=num_weyl(M/2,N,multiplicity);
    int D=nchoosek(M,N);

    if(debug)
    {
	    std::cout << "Number of determinants: " << D << std::endl;
	    std::cout << "Number of Weyl tableau: " << nweyl << std::endl;
    }

    std::vector<int> frame;
    frame.clear();
    frame.push_back(2);
    //frame.push_back(1);

    std::vector<int> T;
    get_init_tableau(M/2,frame,T);

    for(t : T)
	   std::cout << t << " ";
    std::cout << "\n";

    for(int i=0; i<nweyl-1; i++)
    {
	get_next_tableau(M/2,frame,T);
   	for(t : T)
	   std::cout << t << " ";
        std::cout << "\n";
	   
    }


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

int
main2(int argc, char *argv[])
{ 
/*
    using std::chrono::high_resolution_clock;

    std::cout.precision(10);
    std::cout << std::scientific;

    int debug=2;

    //timing variables
    //std::chrono::system_clock::time_point start_time;
    std::chrono::duration<double> time_elapsed;

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

    int M,N;
    Matrix h2;
    std::vector<double>  h4;


    if(get_integrals(basis_fname, nuc_fname, M, N,  h2, h4)!=0)
    {
	    std::cout << "Problem getting integrals\n";
    }

    // ****************************
    // * Output files             *
    // ****************************
    std::fstream fonebody;
    system("touch ham_ov.dat");
    system("rm ham_ov.dat");
    system("touch ham_ov.dat");
    fonebody.open("ham_ov.dat");

    // Format ham_ov.dat-
    // number of basis
    // Ov matrix
    // Core Hamiltonian matrix
    // Nuclear Repulsion

    std::fstream ftwobody;
    system("touch 2_ele.dat");
    system("rm 2_ele.dat");
    system("touch 2_ele.dat");
    ftwobody.open("2_ele.dat");
    // Format 2_ele.dat-
    // number of reduced integral
    //  i,j,k,l, ee(i,j,k,l)

    //increase the printout precision
    fonebody.precision(10);
    ftwobody.precision(10);
    fonebody << std::scientific << std::showpos;
    ftwobody << std::scientific;

    fonebody << M << "\n";
    fonebody << h2 << "\n";

    //temporary header to be overwritten by the number of integrals
    ftwobody  << "h4       \n";
    int nints=0;
    for(auto p=0; p!=M/2; p++) //unique integral labels, looping scheme from libint
            for(auto q=0; q<=p; q++)
                for(auto r=0; r<=p; r++)
                    for(auto s=0; s<= (p==r ? q : r) ; s++)
                        if(std::abs(h4[term4(p,q,r,s)]) > 1e-5)
			{
				//std::cout << term4(p,q,r,s) << ": ["
				//	  << p << q << "|" << r << s << "] = " 
				//	  << h4[term4(p,q,r,s)] << std::endl;

 			        ftwobody << p << " " << q << " " << r << " " << s << " " 
				         << h4[term4(p,q,r,s)] << std::endl;
				nints++;
			}
    //go back to write the number of integrals at the top
    ftwobody.seekp(0);
    //header
    ftwobody  << nints;

    //number of 2e- integrals
    int m=(M/2) - 1;
    int nints=term4(m,m,m,m)+1;

    fonebody.close();
    ftwobody.close();



    //FCI matrix dimension
    int D=nchoosek(M,N);

    //LAPACK
    auto start_time=std::chrono::high_resolution_clock::now();


    //here is the basic CI algorithm
    double* w=ci_eigenvals(M,N,h4,h2,debug);

    time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
    std::cout << "Making matrix and getting eigenvalues took " << time_elapsed.count() 
              << "s to get eigenvalues.\n";

    std::cout << "\nEigenvalues( "<< D << " ) :\n --------- \n";

    std::cout.precision(10);
    std::cout << std::scientific;
    for(int i=0; i<D; i++)
        std::cout << w[i]  << "\n";

    */
    return 0;
}
