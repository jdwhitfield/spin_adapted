// Functions for the symmetric group approach
#include<iostream>
#include<vector>
#include"weyl.h"

class perm_a
{
	public:
		std::vector<int> perm;
		double coeff=1;
		int    n=0;
	// TODO: overload the multiply operator
};

void
print_perm(perm_a P)
{
	std::cout << P.coeff << " * ";
	for(auto p : P.perm)
		std::cout << p <<" ";
	std::cout << "\n";
	return;
}

perm_a
eye(int n)
{
	//unit permutation
	perm_a P;
	P.perm.clear();
	for(int k=0; k<n; k++)
		P.perm.push_back(k);
	P.coeff=1;
	return P;
}

//overload the * operator
perm_a 
perm_multiply(const perm_a& P1, const perm_a& P2)          
{

	perm_a result;
	std::vector<int> p_result;

	p_result.clear();
	if(P1.n != P2.n)
	{
		std::cout << "Error in perm_multiply, perms must be same size" 
			  << std::endl;
		return result;
	}

	p_result.resize(P1.n);

	for(int j=0; j<P1.n; j++)
		p_result[j]=P1.perm[P2.perm[j]];

	result.coeff=P1.coeff*P2.coeff;
	result.perm=p_result;
	
	return result;
}

//overload the * operator
std::vector<perm_a> 
multiply(std::vector<perm_a> A, std::vector<perm_a> B)
{
	std::vector<perm_a>  product;
	for( auto p : A )
		for( auto q : B)
			product.push_back(perm_multiply(p,q));

	//simplify using std::algorithm
	//i.e. combine like terms 
	
	return product;
}

std::vector<perm_a>
Srow(std::vector<int> frame,std::vector<int> tableau)
{
	bool debug=true;
	if(debug) std::cout << "here (young::Srow)\n";
	int N=0; //number of electrons
	for(int r=0; r<frame.size(); r++)
		N=N+frame[r];
	
	std::vector<perm_a> S;
	std::vector<perm_a> Sr;
	for(int row=0; row<frame.size(); row++)
	{
		Sr.clear();

		if(frame[row]==1)
			continue;
		perm_a P;
		P.perm.clear();

		int col0=tableau_pos(row,0,frame);
		int col1=tableau_pos(row,1,frame);
		

		P.perm.clear();
		for(int k=0; k< N; k++)
			P.perm.push_back(k);
		P.perm[tableau[col0]]=tableau[col1];
		P.perm[tableau[col1]]=tableau[col0];

		P.coeff=1;

		Sr.push_back(eye(N));
		Sr.push_back(P);
		
		if(debug) std::cout << "S_before (young::Srow)\n";


		S=multiply(S,Sr);
		if(debug) std::cout << "S_before (young::Srow)\n";
	}
	return S;
}

int Ey(std::vector<int> frame_rows,std::vector<int> yT)
{
	std::cout << "here (young::Ey)\n";
	std::vector<perm_a> S=Srow(frame_rows,yT);
	for(int k=0; k<S.size(); k++)
	{
		std::cout << "k : " << k << "\n";
		print_perm(S[k]);
	}
	return 0;
}

int
main()
{
	int n=3;

	perm_a P12;
	P12.n=n;
	P12.perm={1,0,2};
	P12.coeff=.5;

	//std::cout << ".5 * P12\n";
	//print_perm(P12);

	perm_a id;
	id.n=n;
	id.perm={0,1,2};
	id.coeff=-1;

	//std::cout << "Id\n";
	//print_perm(id);

	std::vector<perm_a> pvec;
	pvec.clear();
	pvec.push_back(id);


	//print_perm(perm_multiply(P12,P12));

	std::vector<int> frame_rows;
	frame_rows.clear();
	frame_rows.push_back(2);
	frame_rows.push_back(2);
	frame_rows.push_back(1);

	std::vector<int> young_tableau;
	young_tableau.clear();
	young_tableau.push_back(0);
	young_tableau.push_back(1);
	young_tableau.push_back(2);
	young_tableau.push_back(3);
	young_tableau.push_back(4);
		

	Ey(frame_rows,young_tableau);
	
	return 0;
}

