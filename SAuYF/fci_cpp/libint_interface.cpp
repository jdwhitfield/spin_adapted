/******************************************************************************
 * 
 * libint_interface.cpp
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


#include "libint_interface.h"

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

        Eigen::SelfAdjointEigenSolver<Matrix> es(S);
        Matrix s=Matrix(es.eigenvalues());
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
            C=Matrix(es.eigenvectors()*shalf);
        else  	       //if symmetric W=U^+
            C=Matrix(es.eigenvectors()*shalf*es.eigenvectors().adjoint());

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
void
get_libints(const libint2::BasisSet& basisset, const std::vector<libint2::Atom>& atoms,
            Matrix &S, Matrix &V, Matrix &T, vector<double>& AOInts)
{

	using namespace libint2;
	S=compute_1body_ints(Operator::overlap,basisset);
	T=compute_1body_ints(Operator::kinetic,basisset);
	V=compute_1body_ints(Operator::nuclear,basisset,atoms);
	AOInts=compute_2body_ints(basisset);

	return;
}



Matrix
compute_1body_ints(libint2::Operator type, 
		   const libint2::BasisSet bfs, const std::vector<libint2::Atom>& atoms)
{

	using namespace libint2;
	using namespace std;

	// maps shell index to basis function index
	auto shell2bf=bfs.shell2bf(); 

	//place to store output
  	const auto n = bfs.nbf();
  	Matrix result(n,n);

	libint2::initialize();  // safe to use libint now


	//ENGINE SETUP
	Engine engine(type,bfs.max_nprim(),bfs.max_l());
	if(type== Operator::nuclear) //put in point charges when needed
	{
		engine.set_params(make_point_charges(atoms));
	}
	//this is a pointer to results of computation updated after compute() calls
	const auto& result_buf = engine.results();

	for(auto s1=0; s1!=bfs.size(); ++s1) {
		//first bf in each shell
		auto bf1=shell2bf[s1];
		//number of bf in each shell
		auto n1= bfs[s1].size();
		   for(auto s2=0; s2!=bfs.size(); ++s2) {
			engine.compute(bfs[s1], bfs[s2]);

			auto ints_shellset = result_buf[0];  // location of the computed integrals
			if (ints_shellset == nullptr)
				continue;  // nullptr returned if the entire shell-set was screened out

			auto bf2=shell2bf[s2];
			auto n2= bfs[s2].size();

			// integrals are packed into ints_shellset in row-major (C) form
			// this iterates over integrals in this order
			//
			//For example, if the number of functions in each shell is na, nb, nc, and nd,respectively,
			//then the integral (ab|cd) is found at position abcd =((anb+b)nc+c)nd+d.
		
			// "map" buffer to a const Eigen Matrix, and copy it to the corresponding blocks of the result
			// 
			// JDW: I really don't know what this code does and I got it from hartree-fock.cc in libint test suite
			Eigen::Map<const Matrix> buf_mat(result_buf[0], n1, n2);
			result.block(bf1, bf2, n1, n2) = buf_mat;
			if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!
			result.block(bf2, bf1, n2, n1) = buf_mat.transpose();

			//for(auto f1=0; f1!=n1; ++f1)
			//	for(auto f2=0; f2!=n2; ++f2)
			//		cout <<"<"<< bf1+f1 << "|" << bf2+f2 << ">=" << ints_shellset[f1*n2+f2] << endl;
	   }
	}

	libint2::finalize();  // do not use libint after this

	return result;
}

std::vector<double>
compute_2body_ints(const libint2::BasisSet bfs)
{
	using namespace libint2;
	using namespace std;

	// maps shell index to basis function index
	auto shell2bf=bfs.shell2bf(); 

	//number of spatial orbitals (including px, py, pz separately)
	const auto M_spatial=bfs.nbf();

	//number of 2e- integrals
	int nints =term4(M_spatial-1,M_spatial-1,M_spatial-1,M_spatial-1)+1;

	//set up array of ints
	std::vector<double> AOInts(nints, 0.0);

	initialize();
 
	// this engine will compute electronic repulsion integrals
	Engine eri_engine(Operator::coulomb, bfs.max_nprim(), bfs.max_l());

	//pointer to the results, gets updated later
	const auto& results = eri_engine.results();

	//this is a pointer to results of computation updated after compute() calls
	const auto& eri_buf_vec = eri_engine.results();


	for(auto s1=0; s1!=bfs.size(); ++s1) {
		for(auto s2=0; s2!=bfs.size(); ++s2) {
			for(auto s3=0; s3!=bfs.size(); ++s3) {
			   for(auto s4=0; s4!=bfs.size(); ++s4) {

				eri_engine.compute(bfs[s1],bfs[s2],bfs[s3],bfs[s4]);
				auto ints_shellset = eri_buf_vec[0];  // location of the computed integrals
				if (ints_shellset == nullptr)
					continue;  // nullptr returned if the entire shell-set was screened out

				auto bf1=shell2bf[s1];
				auto bf2=shell2bf[s2];
				auto bf3=shell2bf[s3];
				auto bf4=shell2bf[s4];
				//number of bf in each shell
				auto n1= bfs[s1].size();
				auto n2= bfs[s2].size();
				auto n3= bfs[s3].size();
				auto n4= bfs[s4].size();

				// integrals are packed into ints_shellset in row-major (C) form
				// this iterates over integrals in this order
				//
				//For example, if the number of functions in each shell is na, nb, nc, and nd,respectively,
				//then the integral (ab|cd) is found at position abcd =((a nb+b) nc+c) nd+d.
				for(auto f1=0; f1!=n1; f1++)
				{
					for(auto f2=0; f2!=n2; f2++)
					{
						for(auto f3=0; f3!=n3; f3++)
						{
							for(auto f4=0; f4!=n4; f4++)
							{
                    		//unpack and save the data the easy way but with a lot of overwriting
                            AOInts[term4(f1+bf1,f2+bf2,f3+bf3,f4+bf4)]=
                              ints_shellset[ f4+ n4 * (f3 + n3*(f2 + n2 *f1)) ];

							//chemist notation: (ab|cd) = \int\int \phi_a(1) \phi_b(1) r^{-1}_{12} \phi_c(2) \phi_d(2) d1 d2
							//if(ints_shellset[((f1*n2+f2)*n3+f3)*n4+f4])
							//cout <<"("<< bf1+f1 << " " << bf2+f2 
							//     <<"|"<< bf3+f3 << " " << bf4+f4 << ")" 
				     		//	   <<"="<<ints_shellset[((f1*n2+f2)*n3+f3)*n4+f4] << endl;

							}
						}
					}
				}

			   }
			}
		}
	}



	finalize();  // do not use libint after this
	



	return AOInts;

}







