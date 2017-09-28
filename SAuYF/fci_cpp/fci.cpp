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
#define EIGEN_USE_LAPACKE
#include"parser.h"
#include"ci_matrix.h"
#include"weyl.h"
#include"libint_interface.h"

int
get_integrals(const char* basis_fname, const char* nuc_fname, int& M, int& N, Matrix& h2, std::vector<double>& h4)
{

    bool debug=true; 
    bool CANONICAL=true;// if false use symmetric orthogonalization

    //parse files
    auto atoms = parse_nucfile(nuc_fname,N);
    libint2::BasisSet shells = parse_basisfile(basis_fname);

    //set M
    M=2*shells.nbf();

    //internal matricies
    Matrix hT,hV;
    Matrix S, h;
    std::vector<double> AOInts;

    //return variables
    //    Matrix h2;
    //    std::vector<double> MOInts;

    auto start_time = std::chrono::high_resolution_clock::now();

    //libints
    get_libints(shells,atoms,S,hV,hT,AOInts);
    h=hT+hV;

    auto end_time=std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_elapsed= end_time-start_time;

    if(debug)
	    std::cout << "Took " << time_elapsed.count() << " seconds to get AO Ints\n";

    if(Matrix(S).isIdentity(1e-9))
    {
        //nothing changes if S is the identity matrix
        h2=h;
        h4=AOInts;
	return 0;
        //std::cout<<"We don't need to orthogonalize\n";
    }
    else
    {// *     ORTHOGONALIZATION    *


        /***********************************************************/ 
        // Transformation matricies
        //
        //Basis transform to an orthogonal space
	//X=Us^{-1/2} when S=U^+ s U
        //C = C(W) = XW for any WW^\dag =\id
        //

	if(debug)
		std::cout << S << std::endl;

	if(debug)
		std::cout << S.eigenvalues() << std::endl;

        Eigen::SelfAdjointEigenSolver<Matrix> Eigensystem(S);

	if(debug)
		std::cout << "here\n";

        Matrix s=Matrix(Eigensystem.eigenvalues());
        //check which eigenvalues are non-trivial to avoid linear dependence
        //then compute s^{-1/2}
        Matrix shalf;
        shalf.resize(M/2,M/2);    // set size
        shalf.fill(0);            // make sure all elements are zero
        //make sure its a column vector
        if(s.rows()<s.cols())
            s.transposeInPlace();

        for(int j=0; j<s.rows(); j++)
        {
            if(std::abs(s(j,0))>1e-6)
                shalf(j,j)=1/sqrt(s(j,0));
            else
                shalf(j,j)=0;
        }

      
        Matrix C;
        if(CANONICAL)  //if canonical W=\id
            C=Matrix(Eigensystem.eigenvectors()*shalf);
        else  	       //if symmetric W=U^+
            C=Matrix(Eigensystem.eigenvectors()*shalf*Eigensystem.eigenvectors().adjoint());

	if(debug) //time the transformations
		start_time = std::chrono::high_resolution_clock::now();

        h2       =Matrix(C.adjoint())*h*C;   //One-body transform
        h4       =transform4(M,C,AOInts);//two-body transform

	if(debug)
	{
		time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
		std::cout << "Transformation to (canonical) orthogonal basis took " << time_elapsed.count() << "ms \n";
	}

    }


    return 0;


}

int
main(int argc, char *argv[])
{ 
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
    /*
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
    */


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

    return 0;
}
