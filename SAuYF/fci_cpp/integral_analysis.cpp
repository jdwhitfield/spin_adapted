#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<assert.h>
#include<string>
#include<iomanip>
#include<random>
#include<chrono>
#include<algorithm>
#include"parser.h"
#include"libint_interface.h"
#include"ci_matrix.h"

using std::string; 
using std::chrono::high_resolution_clock;
using std::cout;
using std::endl;
using std::ostream;
using std::vector;

void
output_to_matlab(ostream& out, Matrix h2, std::vector<double> h4, string name);


int
main(int argc, char *argv[])
{ 

    int debug=1;


    //timing variables
    std::chrono::system_clock::time_point start_time;
    std::chrono::duration<double> time_elapsed;

    // ****************************
    // * Parse commandline inputs *
    // ****************************
    if (argc > 1)
    {	
	    if(!strcmp(argv[1],"-h"))
	    {
            std::cout << "integral_analysis [basis_func] [nuc_field]";
            std::cout << "\n\tDefault files are used when called with no parameters.\n\n";
		    return 0;
	    }
    }
 

    const auto basis_fname=(argc>1) ? argv[1] : "basis_funcs";
    const auto nuc_fname  =(argc>2) ? argv[2] : "nuc_field";

    // ****************************
    // * Parse data files         *
    // ****************************
    
    auto shells = parse_basisfile(basis_fname);


    if(debug)
    {
        printf("Basis set\n");
        for( auto s : shells)
            std::cout << s << "\n";
    }

    //number of spin orbitals in the basis set
    int M=2*shells.nbf();
    int M_spatial=shells.nbf();
    
    //number of electrons
    int N=-99;

    //std::cout << "\nNints="<<term4(M_spatial-1,M_spatial-1,M_spatial-1,M_spatial-1)+1 << "\n";

    /*
    //number of 2e- integrals
    int m=(M/2) - 1;
    int nints=term4(m,m,m,m)+1;
    //FCI matrix dimension
    int D=nchoosek(M,N);
    */

    auto atoms = parse_nucfile(nuc_fname,N);

    if(N==-99)
    {
        std::cout << "Nelec not set\n";
        throw;
    }
    
    std::cout << "\nNelec=" << N << std::endl;
    //output header
    /*
    cout << "\%p(h_{pp}) &=& u_{[-10,0]}(h_{pp}) \n"
         << "\%p(h_{pq}) &=& u_{[-1,1]}(h_{pq}) \n"
         << "\%p(h_{pqqp}) &=& u_{[-0.5,0.5]}(h_{pqqp}) \n"
         << "\%p(h_{pqqr}) &=& \\frac{1}{2 \\cdot 0.2} e^{-|h_{pqqr}|/0.2} \n"
         << "\%p(h_{pqrs}) &=& \\frac{1}{2 \\cdot 0.1} e^{-|h_{pqrs}|/0.1}\n";
    */

    Matrix S, hT,hV;
    std::vector<double> AOInts;

    //libints
    start_time = high_resolution_clock::now();
    get_libints(shells,atoms,S,hV,hT,AOInts);
    time_elapsed=(high_resolution_clock::now()-start_time);
    if(debug)
	    std::cout << "Took " << time_elapsed.count() << " seconds to get AO Ints\n";

    std::cout <<"\n with M_spatial="    << M_spatial     << " should have " 
              << term4(M_spatial-1,M_spatial-1,M_spatial-1,M_spatial-1)+1 
              << " integrals and we have " << AOInts.size() << " ints\n";

    

    /*
    std::ofstream file("out2"); 
    output_to_matlab(file,hT,AOInts,"_pos");

    return 0;
    */

    if(debug)
    {
        std::cout << "S\n" << S << "\n";
        std::cout << "hV\n" << hV << "\n";
        std::cout << "hT\n" << hT << "\n";
        std::cout << "AO 2body integrals\n";
        for(auto p=0; p!=M/2; p++) //unique integral labels, looping scheme from libint
            for(auto q=0; q<=p; q++)
                for(auto r=0; r<=p; r++)
                    for(auto s=0; s<= (p==r ? q : r) ; s++)
                        if(std::abs(AOInts[term4(p,q,r,s)]) > 1e-5)
                            printf("%i: [%i%i|%i%i] = %lf\n",
                                    term4(p,q,r,s),p,q,r,s,AOInts[term4(p,q,r,s)]);
    }

    Matrix h=hT+hV;
    Matrix h2;
    std::vector<double> MOInts;

    //now we'll bin the one-body integrals
    double bin_hpp[M_spatial];
    double bin_spq[M_spatial*(M_spatial-1)/2];
    double bin_hpq[M_spatial*(M_spatial-1)/2];


    if(debug) cout << "Starting 1-body sorting..." << endl;

    int ctr=0;
    for(auto p=0; p<M_spatial; p++) 
    {
	    for(auto q=p; q<M_spatial; q++) 
	    {
		    if(p==q)
			    bin_hpp[p]=h(p,p);
		    else
		    {
			    bin_spq[ctr]=S(p,q);
			    bin_hpq[ctr]=h(p,q);
			    ctr++;
		    }
	    }
    }

    if(debug)
    {
	    cout << "...Done 1-body sorting" << endl;
	    cout << "bin_hpp: ";
	    for( auto i : bin_hpp) cout << i << " ";
	    cout << endl;


	    cout << "bin_hpq: ";
	    for( auto i : bin_hpq) cout << i << " ";
	    cout << endl;

	    cout << "bin_spq: ";
	    for( auto i : bin_spq) cout << i << " ";
	    cout << endl;
    }



    if(debug)
    {
	    cout << "Now for the two-body integrals:\n";
	    cout << "We will classify them as two integral, three integral or four integral\n"
		 << "For each of these classes one might want to add or subtract the integrals\n" 
		 << "appropriately. For example, the integrals of the form h_{ijji} and h_{ijij}\n"
		 << "often appear as h_{ijji}-h_{ijij}.\n"
		 << "\nI won't accommodate this right now.  We'll just group them as:\n"
		 << "{h_ijji,h_ijij}, {h_iiij} {h_ijjk,h_ijkj}, {h_pqrs,h_pqsr}\n [ ] print integral indices \n [ ] coverage of integrals\n";
    }


    cout << "two index integrals\n";

    double bin_hpppp[M_spatial*(M_spatial+1)/2];
    double bin_hpqqq[M_spatial*(M_spatial+1)/2];
    double bin_hppqq[M_spatial*(M_spatial+1)/2];
    int ctr2=0;
    for(auto p=0; p<M_spatial; p++) //unique integral labels, looping scheme from libint
    {
         for(auto q=0; q<=p; q++)
	 {
		if(p==q)
		{
		    cout << term4(p,p,p,p)  << " ";
	 	    ctr2++; 
		}
	 }
    }

    cout << "\nCounted " << ctr2 << " two-body integrals with one index." 
	 << endl;

    ctr2=0;
    for(auto p=0; p<M_spatial; p++) //unique integral label loop from libint
    {
         for(auto q=0; q<=p; q++)
	 {
	        if(p!=q)
		{
			cout << term4(p,q,p,q)  << " " << term4(p,p,q,q)  << " ";
			cout << term4(p,q,q,q)  << " " << term4(q,p,p,p) << " ";
			ctr2++;
			ctr2++;
			ctr2++;
			ctr2++;
		}
	 }
    }
    cout << "\nCounted " << ctr2 << " two-body integrals with two indices." 
	 << endl;

    ctr2=0;
    for(auto p=0; p<M_spatial; p++) //unique integral label loop from libint
    {
         for(auto q=0; q<=p; q++)
	 {

	     if(p!=q)
	     {
                for(auto r=0; r<=p; r++)
	        {
	           if(r!=p && r<q)
		   {
			cout << "*p=" << p ;
			cout << ",q=" << q ;
			cout << ",r=" << r <<"* ";
			cout << term4(p,q,q,r)  << " " << term4(p,r,q,q)  << " ";
			cout << term4(p,p,q,r)  << " " << term4(p,r,p,q)  << " ";
			cout << term4(p,r,r,q)  << " " << term4(p,q,r,r)  << " ";
			ctr2++;
			ctr2++;
			ctr2++;
			ctr2++;
			ctr2++;
			ctr2++;
		   }
		}
	     }
	 }
    }

    cout << "\nCounted " << ctr2 << " two-body integrals with three indices." 
	 << endl;

    double bin_hpqqr[M_spatial*(M_spatial-1)/2];
    double bin_hpqrs[M_spatial*(M_spatial-1)*(M_spatial-2)*(M_spatial-3)/192];

    //number of 2e- integrals
    int nints =term4(M_spatial-1,M_spatial-1,M_spatial-1,M_spatial-1)+1;

//------------------------------------------------------------------------------
    cout << nints << endl;
    cout << "N_hpqrs : " 
	 << M_spatial*(M_spatial-1)*(M_spatial-2)*(M_spatial-3)/192 
	 << endl;

    int ctr_pq=0;
    int ctr_pqqp=0;
    int ctr_pqqr=0;
    int ctr_pqrs=0;
    for(auto p=0; p<M_spatial; p++) //unique integral labels, looping scheme from libint
    {
         for(auto q=0; q<=p; q++)
	 {

	     //cout << term4(p,p,q,q) << " ";
             ctr_pqqp++;
		
             if(p!=q)
	     {
                 ctr_pq++;
	     }

             for(auto r=0; r<=p; r++)
	     {
		 if(r!=p && r!=q)
		 {
		     cout << term4(p,r,r,q)  << " " << term4(p,q,r,r) << " " ;
		     ctr_pqqr++;
		 }

                 for(auto s=0; s<= (p==r ? q : r) ; s++)
			ctr_pqrs++;
	     }
	 }
    }
    cout << endl;
 
    /*
    for(int p=0; p<M_spatial; p++)
    {
    	for(int q=p+1; q<M_spatial; q++)
	{

		bin_spq[ctr_pq]=S(p,q);
		bin_hpq[ctr_pq]=h(p,q);
		cout << term4(p,p,q,q) << " ";

		//cout << ctr_pq << " ";

		for(int r=0; r<M_spatial; r++)
		{
			//bin_hpqqr[ctr_pqqr]=
			if(r==p || r==q)
				continue;
			ctr_pqqr++;
			for(int s=0; s<M_spatial; s++)
			{
				if(r==s || s == q || s ==p)
					continue;
				ctr_pqrs++;

			}
		}
	}

    }
    */
    cout << ctr_pq << " ?= " << M_spatial*(M_spatial-1)/2 << endl;
    cout << ctr_pqqp << " ?= " << M_spatial*M_spatial << endl;
    cout << ctr_pqqr << " ?= " << M_spatial*(M_spatial-1)*(M_spatial-2)/2 << endl;
    cout << ctr_pqrs << " ?= " << M_spatial*(M_spatial-1)*(M_spatial-2)*(M_spatial-3)/192 << endl;
    /*
    cout << --ctr_pq << " ?= " << M_spatial*(M_spatial-1)/2 << endl;

    */

    // ***********************************************************/ 
    // *     ORTHOGONALIZATION    *
    // ***********************************************************/ 
    // Transformation matrices

    /*
    //Basis transform to an orthogonal space
    //C = C(W) = XW for any WW^\dag =\id
    Eigen::SelfAdjointEigenSolver<Matrix> EigensystemS(S);
    Matrix s=Matrix(EigensystemS.eigenvalues());
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

    Matrix X=Matrix(EigensystemS.eigenvectors()*shalf);
    //Basis transform to an orthogonal space with X=U_S s^{-1/2}
    //C = C(W) = XW for any WW^\dag =\id


    //CANONICAL ORTHOGONALIZATION basis using W=\id
    Matrix C;
    bool CANONICAL=true;
    if(!CANONICAL)
        X=Matrix(EigensystemS.eigenvectors()*shalf*EigensystemS.eigenvectors().adjoint());

    start_time = high_resolution_clock::now();
    h2       =Matrix(X.adjoint())*h*X;   //One-body transform
    Matrix h2_V      =Matrix(X.adjoint())*hV*X;   //One-body transform of potential
    MOInts   = transform4(M,X,AOInts);   //two-body transform
    time_elapsed = high_resolution_clock::now() - start_time;

    string name("_canonical");
    std::ofstream file1("ints_canonical.m"); 
    output_to_matlab(file1,h2,MOInts,name);


   
    /*
    //SYMMETRICAL ORTHOGONALIZATION with W=U_S
    Eigen::SelfAdjointEigenSolver<Matrix> pos_basis(h2_V);
    Matrix W = EigensystemS.eigenvectors();

    C=X*W;
    
    start_time      = high_resolution_clock::now();
    h2              = Matrix(C.adjoint())*h *C;   //One-body transform
    Matrix h2_V_pos = Matrix(C.adjoint())*hV*C;   //One-body transform
    Matrix h2_T_pos = Matrix(C.adjoint())*hV*C;   //One-body transform
    MOInts          = transform4(M,C,AOInts);     //Two-body transform
    time_elapsed    = high_resolution_clock::now() - start_time;

    std::ofstream file2("ints_posbasis.m"); 
    output_to_matlab(file2,h2,MOInts,"_pos");



    
    //POSITION BASIS with W=eigvect(V) 
    Eigen::SelfAdjointEigenSolver<Matrix> pos_basis(h2_V);
    Matrix W = pos_basis.eigenvectors();

    C=X*W;
    
    start_time      = high_resolution_clock::now();
    h2              = Matrix(C.adjoint())*h *C;   //One-body transform
    Matrix h2_V_pos = Matrix(C.adjoint())*hV*C;   //One-body transform
    Matrix h2_T_pos = Matrix(C.adjoint())*hV*C;   //One-body transform
    MOInts          = transform4(M,C,AOInts);     //Two-body transform
    time_elapsed    = high_resolution_clock::now() - start_time;

    std::ofstream file2("ints_posbasis.m"); 
    output_to_matlab(file2,h2,MOInts,"_pos");

    //RANDOM BASIS 

    //random engine initialization
    std::default_random_engine generator{};
    unsigned rand_seed=0;
    if(rand_seed==0)
        rand_seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed(rand_seed);
 
    //get a Haar random unitary
    W=randU(h2.cols(),generator);
    C=X*W;
    
    start_time      = high_resolution_clock::now();
    h2              = Matrix(C.adjoint())*h *C;   //One-body transform
    h2_V_pos = Matrix(C.adjoint())*hV*C;   //One-body transform
    h2_T_pos = Matrix(C.adjoint())*hV*C;   //One-body transform
    MOInts          = transform4(M,C,AOInts);     //Two-body transform
    time_elapsed    = high_resolution_clock::now() - start_time;

    std::ofstream file3("ints_randbasis.m"); 
    output_to_matlab(file3,h2,MOInts,"_rand");
    */

    return 0;
}



void
output_to_matlab(ostream& out, Matrix h2, std::vector<double> h4, string name)
{
    out.precision(10);
    out << std::scientific;


    int M=h2.cols();

    out << "\nhpp" << name << "=[";
    for(int p=0; p<M; p++)
        out << h2(p,p) << " ";
    out <<"];\n";

    out << "\nhpq" << name << "=[";
    for(int p=0;p<M;p++)
        for(int q=p+1;q<M;q++)
            out << h2(p,q) << " ";
    out<<"];\n";

    vector<double> hpqqp;
    vector<double> hpqqr;
    vector<double> hpqrs;
    int ctr=0;
    for(int p=0; p<M; p++)
        for(int q=0;q<=p; q++)
            for(int r=0;r<=p; r++)
                for(int s=0;s<=(p==r?q:r); s++)
                {
                    auto flag=1;

                    if((p==q && r!=s) || (p!=q && r==s))
                    {
                        hpqqr.push_back(h4[term4(p,q,r,s)]);
                       
                        //hpqqr.push_back(h4[term4(p,s,q,r)]);
                      
                        //elem=distro_hpqqr(generator);
                        flag=0;
                    }

                    if(p==q && r==s)
                    {
                        ctr++;
                        //elem=  distro_hpqqp(generator);
                        //elem*= 2*(randsign(generator)-.5);
                        hpqqp.push_back(h4[term4(p,q,r,s)]);
                        flag=0;
                    }

                    if(flag)
                    {
                        hpqrs.push_back(h4[term4(p,q,r,s)]);
                    }
                    
                }

    std::cout << "\n Printing "<< ctr << " hpqqp integrals\n";
    out << "hpqqp"<< name << "=[";
    for( auto h : hpqqp )
        out << h << " ";
    out << "];\n";

    out << "hpqqr"<< name << "=[";
    for( auto h : hpqqr )
        out << h << " ";
    out << "];\n";

    out << "hpqrs"<< name << "=[";
    for( auto h : hpqrs )
        out << h << " ";
    out << "];\n";


    return;
}
