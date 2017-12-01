#include<stdio.h>
#include<iostream>
#include<algorithm>
#include"weyl.h"
#include"young.h"

int
Eij(int i, int j, std::vector<int> Tl, std::vector<int> Tr)
{

	if( Tl == Tr)
	{
		if(i==j)
		{
			//count number of i in T
			int ctr=0;
			for(int ii=0; ii<Tl.size(); ii++)
				if(Tl[ii]== i) ctr++;

			std::cout << "E_{"<<i <<","<< j << "}(T,T)=" << ctr;
		}
	}
	return 0;
}

void
gelfand(int M, std::vector<int> tableau, const std::vector<int> frame_rows)
{
	std::vector<int> gm;
	std::vector<int> f=frame_rows;
	auto T=tableau;
	std::vector<int> f_reduced;
	std::vector<int> T_reduced;
	
	std::cout << "In gelfand\n";
	print_tableau(f,T);

	std::cout << "\n";
	for(int rowg=0; rowg < M ; rowg++)
	{

		int M_max=M-rowg;


		//subduct
		int removed=0;
		T_reduced.clear();
		f_reduced=f;
		for( int i=0; i<T.size(); i++)
		{
			int entry=T[i];
			if(entry==M_max)
			{
				int row;
				int col;
				
				pos_to_row_col(i,f,row,col);

				//remove box from that row of the frame
				f_reduced[row]--;

			}
			else
				T_reduced.push_back(entry);

		}
		//assign for next loop
		T=T_reduced;
		f=f_reduced;

		//std::cout << rowg << " : ";
		//std::cout << "frame: ";
		int ctr=0;
		for( auto row :f)
		{
			if(ctr< M-rowg)
				std::cout << row << " ";
			ctr++;
		}
		if(M-f.size()-rowg>0)
			for(int j=0; j< (M-(int)f.size() -rowg); j++)
				std::cout << "0 ";
		std::cout << "\n";
  		//std::cout << "Tableau:";
		//print_tableau(f,T);


	

	}
}

//Gelfand tableau entries
int 
mij(int i, int j, int M, std::vector<int> T,std::vector<int> frame_rows)
{
	//m_ij is the ith row after the j-1 largest values have been removed 
	//                    after removing M-j values

	int m=0;

	if(i > frame_rows.size() - 1)
	{
		return 0;
	}

	m=frame_rows[i];

	//for col in   > M-1 ) m--; 

	
	for(int col=0; col<frame_rows[i]; col++) //j=0, deduct nothing
	                                         //j=1, deduct values of M-1	    i.e. discount T[i,j]>M-2
	                                         //j=2, deduct values of M-1, M-2 i.e. discount T[i,j]>M-3
	                                         //j=s, deduct values of M-1, M-2, ..., M-s i.e. discount T[i,j]>M-(s+1)
		if( T[tableau_pos(i,col,frame_rows)] > M-(j+1)) m--;


			
	return m;
}

std::vector<basis_func>
excitation_op(int occ, int vir, const std::vector<basis_func> F)
{
    bool debug=false;
    std::vector<int> P;
    std::vector<basis_func> G=F;
    bool found_occ=false;
    
    int n=G.size();
    //sum over basis functions
    for(int i=0; i<n; i++)
    {
    	//copy in
    	P=G[i].orbs;
    
    	//check for occupied, replace with virtual
    	//loop over orbitals in N-body basis functions
    	for(int j=0;j<P.size(); j++)
    	{
    	    if(P[j]==occ)
    	    {
    	    	if(debug) std::cout 
    	    		  << "found occupied orbital at index " 
    	    		  << j << "\n";
    	    	P[j]=vir;
    	    	found_occ=true;
    	    	if(j<P.size()-1) 
    	    	    for(int k=j+1;k<P.size();k++)
    	    	    {
    	    	    	if(P[k]==occ)
    	    	    	{
    	    	    	    if(debug) std::cout 
    	    	    		      << "found occupied orbital at index " 
    	    	    		      << j+1 << "\n";

    	    	    	    basis_func bf=G[i];
    	    	    	    bf.coeff=G[i].coeff;
    	    	    	    //revert previous change
    	    	    	    bf.orbs[j]=occ;
    	    	    	    //update changed orbital
    	    	    	    bf.orbs[j+1]=vir;
    	    	    	    //additional term
    	    	    	    G.push_back(bf);

    	    	    	    if(debug) print_bf(G[G.size()-1]);

			    break;
    	    	    	}
    	    	    }
    	    		
    	    	break;
    	    }
    	}
    
    	//save orbital listing
    	G[i].orbs=P;
    	if(debug) print_bf(G[i]);
    }
    return G;
}

int 
main(int argc, char *argv[])
{
    using std::cout;

    bool debug=true;
    if(debug)
    {
	    std::cout.precision(10);
	    std::cout << std::scientific;
    }

    std::vector<int> F;
    F.clear();
    F.push_back(2);
    F.push_back(1);

    std::vector<basis_func> wf_init= initial_state(3,{0,0,2});
    cout << "wf_init= \n"; 
    print_wf(wf_init);
    cout << "\n";

    auto wf=wf_init;



    std::vector<std::vector<basis_func>> C = get_irrep_basis(F,wf_init);

    cout << "Obtained states that span the irreducible space\n";
    cout << "C = {";
    for(auto c : C)
    {   print_wf(c); cout << "\n"; }
    cout << "}\n";

    cout << "Orthogonalization of this basis is done using the Wigner ops.\n";
    
    int N=3; //number of electrons
    int gS=2;//gS = 2S+1 -> S=1/2
    std::vector<perm_a> W00=wigner_op(0,0,N,gS,C);
    std::vector<perm_a> W10=wigner_op(1,0,N,gS,C);

    cout << "W00= {\n";
    for( auto w : W00)
	    print_perm(w);
    cout << "}\nand\nW10 ={\n";
    for( auto w : W10)
	    print_perm(w);
    cout << "}.\nAs a sanity check consider: W00 W10 and W10 W00 with\nW10 W00={\n";
    for( auto w : multiply(W10,W00))
	    print_perm(w);
    cout <<"}\nand W00 W10\n";
    for( auto w : multiply(W00,W10))
	    print_perm(w);
    cout << "\n";



    cout << "Acting these on the initial wave function gives:\n";

    std::vector<basis_func> projected_0=multiply(W00,wf_init);
    std::vector<basis_func> projected_1=multiply(W10,wf_init);

    cout << "wf_p0 = W00 wf_init :\n";
    print_wf(projected_0);
    cout << "\n and wf_p1 = W10 wf_init :\n";
    print_wf(projected_1);
    cout <<"\n";

    cout << "\nWe can check that they are orthogonal:\n dot(wf_p0,wf_p1) = " 
	 << dot(projected_0,projected_1) << "\n and that they span the same"
	 << " space as C.\n";

    cout << "\nLet us check the other possibilities:\n"
	 << "dot(wf_p0,wf_p1) = " << dot(projected_0,projected_1) << "\n"
	 << "dot(wf_p1,wf_p0) = " << dot(projected_1,projected_0) << "\n"
	 << "dot(wf_p1,wf_p1) = " << dot(projected_1,projected_1) << "\n"
	 << "dot(wf_p0,wf_p0) = " << dot(projected_0,projected_0) << "\n";

    cout << "\nNone are zero but if we consider the adjoint of W00\n"
	 << "dot(wf, W00 wf) = " << dot(wf, multiply(W00,wf))     << "\n"
	 << "dot(W00 wf, wf) = " << dot(multiply(W00,wf), wf)     << "\n";

    cout << "\n\n-------------------------------------\n\n";
    /*************************************************************************/
    //The wave functions are not orthogonal. Why not?

    std::vector<perm_a> A;
    std::vector<perm_a> invA= A;
    A.clear();
    invA.clear();

    perm_a P123; P123.coeff=1; P123.perm={2,0,1};
    perm_a P132; P132.coeff=1; P123.perm={1,2,0};

    A.push_back(P123);
    invA.push_back(P132);

    double lambda=1;
    A[0].coeff=lambda;

    invA[0].coeff=1.0/lambda;

    auto Awf=multiply(A,wf_init);
    auto invAwf=multiply(invA,wf_init);

    cout << "As a test of the dot product, using matrix A and its inverse A^-1.\n"
         << "Apply A = a I to the left and A^-1=(1/a) I to the right and take\n"
	 << "the inner product of the wave function is: <wf|wf> =\n\t";
    cout << dot(wf_init,wf_init) << "\n\n"; 
    cout << "Compare with <A wf | invA wf > =\n\t";
    cout << dot(multiply(A,wf_init),multiply(invA,wf_init)) << "\n\n"; 
    cout << "And compare\n<A wf | A wf > =\t";
    cout << dot(multiply(A,wf_init),multiply(A,wf_init)) << "\n"; 
    cout << "with\n<wf | A(A wf) > =? <wf | (A^2) wf >\t";
    cout << dot(wf_init,multiply(A,multiply(A,wf_init))) << " =? " 
	 << dot(wf_init,multiply(multiply(A,A),wf_init)) << "\n\n";
    cout << "and\n<wf | A(A wf) > =? <wf | (A^2) wf >\t";
    cout << dot(wf_init,multiply(A,multiply(A,wf_init))) << " =? " 
	 << dot(wf_init,multiply(multiply(A,A),wf_init)) << "\n\n";



    basis_func f;
    f.coeff= 1;
    f.orbs.clear();
    f.orbs.push_back(1);f.orbs.push_back(0);f.orbs.push_back(0);

    wf.clear();
    wf.push_back(f);

    cout << "Let us consider a simpler wave function.\nwf:";
    print_wf(wf);
    cout << "The projection via W00.\nW00 wf:";
    print_wf(multiply(W00,wf));
    cout << "\nThe projection via W10.\nW10 wf:";
    print_wf(multiply(W10,wf));
    cout << "\n\n-------------------------------------\n\n";



    Awf=multiply(A,wf);
    invAwf=multiply(invA,wf);


    cout << "As a test of the dot product, using matrix A and its inverse A^-1.\n"
         << "Apply A = a I to the left and A^-1=(1/a) I to the right and take\n"
	 << "the inner product of the wave function is: <wf|wf> =\n\t";
    cout << dot(wf,wf) << "\n\n"; 
    cout << "Compare with <A wf | invA wf > =\n\t";
    cout << dot(multiply(A,wf),multiply(invA,wf)) << "\n\n"; 
    cout << "And compare\n<A wf | A wf > =\t";
    cout << dot(multiply(A,wf),multiply(A,wf)) << "\n"; 
    cout << "with\n<wf | A(A wf) > =? <wf | (A^2) wf >\t";
    cout << dot(wf,multiply(A,multiply(A,wf))) << " =? " 
	 << dot(wf,multiply(multiply(A,A),wf)) << "\n\n";

    cout << "printing A and A^2\nA={\n";
    for(auto p : A) print_perm(p);
    cout << "}\nA^2={";
    for(auto p : multiply(A,A)) print_perm(p);
    cout << "}\n\n";

    cout << "printing wf, A wf, and A^2 wf.\n";
    cout << "wf: \n";       print_wf(wf); cout <<"\n";
    cout << "A wf: \n";     print_wf(multiply(A,wf)); cout <<"\n";
    cout << "A^2 wf: \n";   print_wf(multiply(multiply(A,A),wf)); cout <<"\n";



    cout << "\n\n-------------------------------------\n\n";




    /*************************************************************************/



    cout << 
"Consider P ket{C[0]} = dot(p_0,C[0]) ket{p_0} + dot(p_1,C[0]) ket{p_1}. If \n"
	 << "this is fully in space P then, dot(p_1,C[0])^2 + dot(p_0,C[0])^2\n"
	 << "should be equal to  dot(C[0],C[0]). Ip so facto for C[1].\n";

    cout << "\ndot(p_1,C[0])^2 + dot(p_0,C[0])^2 =\t" 
         << dot(projected_0,C[0])*dot(projected_0,C[0]) 
         + dot(projected_1,C[0])*dot(projected_1,C[0]);
    cout << "\ndot(C[0],C[0])=\t" << dot(C[0],C[0]);
    cout << "\n";

    return 0;

    /*
    //                                i j  
    std::vector<basis_func> wfL0=multiply(W00,C[0]);
    std::vector<basis_func> wfL1=multiply(W10,C[0]);

    std::vector<basis_func> wfL= initial_state(3,{0,0,1});
    std::vector<basis_func> HwfL= excitation_op(1,2,wfL);
    
    cout << "HwfL= ";
    cout << "\n";
    for(auto b : HwfL)
	    print_bf(b);
    cout << "\n";
    

    //consider the inner product of the excited and initial wave functions
    cout << "wfR= ";
    print_wf(wf_init);
    cout << "\n";
    cout << "wfL= "; 
    print_wf(wfL);
    cout << "\n";
    cout << "HwfL= ";
    print_wf(HwfL);
    cout << "\n";
    cout << "\n";

    cout << "dot(wf_init,a_i^\dag a_j wfL)=" << dot(wf_init,HwfL) << "\n";

    return 0;


    int M,N;
    M=5; N=3; 
    std::cout << "M:" << M <<", N:" << N <<"\n";

    int multiplicity;

    
    if( N % 2 ) multiplicity=2;//doublet
    else        multiplicity=1;//triplet

    int nweyl=num_weyl(M/2,N,multiplicity);
    //

    return 0;

    std::vector<std::vector<basis_func>> C = get_irrep_basis(F,wfL);

    int ndets=nchoosek(M,N);
    //int D=nchoosek(M,N);

    if(debug)
    {
        std::cout << "Number of Weyl tableau: " << nweyl << std::endl;
        //std::cout << "Number of determinants: " << D << std::endl;
    }


    std::vector<int> frame_rows;
    frame_rows.clear();
    frame_rows.push_back(2);
    //frame_rows.push_back(1);

    std::vector<std::vector<int>> weyl_list;
    std::vector<int> T;
    get_init_tableau(M/2,frame_rows,T);
    weyl_list.push_back(T);

    for(t : T)
	   std::cout << t << " ";
    std::cout << "\n";

    for(int i=0; i<nweyl-1; i++)
    {

	get_next_tableau(M/2,frame_rows,T);
    	weyl_list.push_back(T);

   	for(t : T)
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

    M=8;
    T.clear();
    T.push_back(0);
    T.push_back(2);
    T.push_back(1);
    T.push_back(5);
    T.push_back(7);
    frame_rows.clear();
    frame_rows.push_back(2);
    frame_rows.push_back(2);
    frame_rows.push_back(1);
    print_tableau(frame_rows,T);

    / *
    M=3;
    T.clear();
    frame_rows.clear();
    frame_rows.push_back(2);
    T.push_back(1);
    T.push_back(2);

    for(int j=M-1; j>-1;j--)
    {
    	for(int i=0; i<j+1;i++)
	{
		    std::cout << "m_"<< i << j << " ";
		    std::cout << "=" << mij(i,j,M,T,frame_rows);
		    std::cout << "\t";
	}

	std::cout << "\n";
    }

    * /
    std::cout << "\n";
    gelfand(M, T, frame_rows);
    */
    

    return 0;
}
