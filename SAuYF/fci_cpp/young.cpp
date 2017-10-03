// Functions for the symmetric group approach
#include<iostream>
#include<vector>
#include<algorithm>
#include<Eigen>
#include"weyl.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        Matrix;  // import dense, dynamically sized Matrix type from Eigen;

class perm_a
{
	public:
		std::vector<int> perm;
		double coeff=1;
		int    n=0;
	// TODO: overload the multiply operator
};

//test functions
int                   test_symmetrizer();
//symmetrizer
std::vector<perm_a>   Srow(std::vector<int> frame,std::vector<int> tableau);


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
perm_I(int n)
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

	bool debug=false;
	perm_a result;
	std::vector<int> p_result;

	p_result.clear();
	if(P1.perm.size() != P2.perm.size())
	{
		std::cout << "Error in perm_multiply, perms must be same size" 
			  << std::endl;
		return result;
	}

	p_result.resize(P1.perm.size());

	for(int j=0; j<P1.perm.size(); j++)
	{
		p_result[j]=P1.perm[P2.perm[j]];
		if(debug)
		{
			std::cout << "p_result element " << j << " "
				  << p_result[j] ;
			std::cout << " (young::perm_multiply)\n";
		}
	}



	result.coeff=P1.coeff*P2.coeff;
	result.perm=p_result;
	
	return result;
}

//overload the * operator
std::vector<perm_a> 
multiply(std::vector<perm_a> A, std::vector<perm_a> B)
{
	bool debug=false;
	std::vector<perm_a>  product;
	for( auto p : A )
		for( auto q : B)
		{
			perm_a pq=perm_multiply(p,q);
			if(debug)
			{
				std::cout << "pq (young::multiply)\n";
				print_perm(pq);
			}
			product.push_back(pq);
		}

	//simplify using std::algorithm
	//i.e. combine like terms 
	
	if(debug)
	{
		std::cout << "product (young::multiply)\n";
		for( perm_a p : product)
			print_perm(p);
	}
	
	return product;
}

int
pos_ymax(const int pos, const int M, const std::vector<int> frame, const std::vector<int> tableau)
{
	int max_val=M-1;
	return max_val;
}

int
pos_ymin(const int pos,const std::vector<int> frame, std::vector<int> tableau)
{

	int above_pos;
	int min_idx=0;
	int row,col;
	bool debug=true;

	//get row from pos
	pos_to_row_col(pos,frame,row,col);

	// strictly monotonic down columns
	if(row_col_to_pos(above_pos,frame,row-1,col))
	{
		min_idx=tableau[above_pos]+1;

		if(debug)
		{
			std::cout << "pos " << pos 
				  << " -> row,col ";
		        pos_to_row_col(pos,frame,row,col);
			std::cout << row << "," << col << "\n"
				  <<"above_pos "<< above_pos
				  << "-> row,col ";
			pos_to_row_col(above_pos,frame,row,col);
			std::cout << row << "," << col << "\n"
			     "tableau[above_pos]+1 " 
			     << tableau[above_pos]+1<< "\n";
		}
	}
	
	// weakly monotonic across rows
	if(col)
	{
		min_idx=std::max(min_idx,tableau[pos-1]);
		if(debug)
		{
			std::cout << "col: "<< col << ", tableau[pos-1] "
				  << tableau[pos-1] << "\npos " << pos 
				  << " -> row,col ";
			pos_to_row_col(pos,frame,row,col);
			std::cout << row << "," << col 
				  << "\n";

			std::cout << "pos-1 -> row,col ";
			pos_to_row_col(pos-1,frame,row,col);
			std::cout << row << "," << col 
				  << "\n";
		}
	}


	if(debug)
	{
		if(min_idx==3)
			std::cout << "\n3 found\n";


		std::cout << "returning "<< min_idx <<"\n\n";
	}

	return min_idx;

}

int
reset_prev_young(const int pos0, const int M, const std::vector<int> frame, std::vector<int>& tableau)
{

	bool debug=false;


	if(pos0==tableau.size()) //already at the end, do nothing, no error
		return 0;


	for(int p=pos0; p<tableau.size(); p++)
	{
		//get the min value allowed here based on tableau
		int min_val=pos_ymin(p,frame,tableau);
		
		//check the value 
		if(min_val>M-1)			//if out of bounds
			return 1;
		else				//if reasonable
			tableau[p]=min_val;
	}
	return 0;
}

int
get_next_ytableau(const std::vector<int> frame, std::vector<int>& tableau)
{
	int ctr=0;

	//number of boxes
	int N=0; 
	for(int r=0; r<frame.size(); r++) 
		N=N+frame[r];

	for(int pos=tableau.size()-1; pos>0-1; pos--)
	{
		if(tableau[pos]==N-1)
		{//at the max allowed, go to next pos
			continue;
		}
		else
		{//not at the max allowed
			//advance tableau position
			tableau[pos]++;

			//try to reset lower positions
			int first_pos_to_be_reset=pos+1;
			auto out=reset_prev_young(first_pos_to_be_reset,N,frame,tableau);
			if(0==out) //everything worked
			       	return 0;
			
			//if it didn't work, just go to next position
			
		}

		//breaker 
		if(ctr++>1000)
			break;
	}

	return 1;
}

int 
get_init_ytableau(const int M, const std::vector<int> frame, std::vector<int>&tableau)
{
	//empty tableau
	tableau.clear();
	
	int n=0;
	for(int row=0; row<frame.size(); row++)
		n=n+frame[row];

	int ctr=0;
	for(int k=0;k<n; k++)
	{
		tableau.push_back(ctr);
		ctr++;
	}

	return 0;
}

Matrix
perm_matrix(std::vector<int> perm)
{
	int n = perm.size();
	Matrix P;
	P.resize(n,n);
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<n; j++)
		{
			if(j==perm[i])
				P(i,j)=1;
			else
				P(i,j)=0;
		}
	}
	return P;
}

int 
sgn(std::vector<int> perm)
{
	return perm_matrix(perm).determinant();
}

int Ey(std::vector<int> frame_rows,std::vector<int> yT)
{
	
	return 0;
}

int
num_young(int N, int ms)
{

	bool debug=false;

	//compute S based on ms=2S+1
	double S=(ms -1)/2.0;
	if(debug)
	{
		std::cout << "  S: " << S << ", degen: " << ms 
			  << ", N: " << N << "\n";
		std::cout << ".5*N-S: " << (int)((.5*N)-S) << "\n";
		std::cout << ".5*N-S-1: " << (int)((.5*N)-S-1) << "\n";
	}


	return nchoosek(N,(int)((.5*N)-S))-nchoosek(N, (int)((.5*N)-S-1));
}

std::vector<int> 
Frows_to_Fcols(std::vector<int> Frows)
{
	std::vector<int> C;
	C.clear();
	C.push_back(Frows.size());
	
	int N=0;
	for( int r : Frows)
		N=N+r;

	//only two columns and the total boxes = N
	C.push_back(N-Frows.size());
	return C;

	
}
bool
isvalid_young(std::vector<int> frame_rows, std::vector<int> tableau)
{
	bool debug=false;

	//check rows, must increase along row 
	for(int r=0; r<frame_rows.size(); r++)
	{
		
		if(frame_rows[r]==2)
		{

			if(debug)
			{
				std::cout << "row: " << r << " entries:\n";
				std::cout << "col0: " 
				          << tableau[ tableau_pos(r,0,frame_rows) ]
					  << "\n";
				std::cout << "col1: " 
				          << tableau[ tableau_pos(r,1,frame_rows) ]
					  << " (young::isvalid_young)\n";
			}


			if( ! (//if its not increasing along the row
			tableau[ tableau_pos(r,0,frame_rows) ]
			< tableau[ tableau_pos(r,1,frame_rows) ]
			)
			  )
				return false;
		}

	}
	
	//check columns
	auto frame_cols=Frows_to_Fcols(frame_rows);
	for(int c=0; c<frame_cols.size(); c++)
	{
		for( int r=0; r<frame_cols[c]-1; r++)
			if(  !(//if its not increasing along the column
					tableau[ tableau_pos(r,c,frame_rows) ]
					< tableau[ tableau_pos(r+1,c,frame_rows) ]
				))
					return false;

	}


	
	return true;
}
int
main()
{
	

	std::vector<int> Yrows;
	Yrows.push_back(2);
	Yrows.push_back(2);
	Yrows.push_back(1);

	std::vector<int> ytab={0,1,2,3,4};
	int ctr=0;
	for(int k=0; k<120; k++)
	{

		if(isvalid_young(Yrows,ytab))
		{

			std::cout << ctr << " : ";
			ctr++;
			for(auto item: ytab)
				std::cout << item << " ";
			std::cout << "\n";

		}
		else
		{

			
			if(false)
			{
				std::cout << "rejected : ";
				for(auto item: ytab)
					std::cout << item << " ";
				std::cout << "\n";
			}
		}

		next_permutation(ytab.begin(),ytab.end());

	}


	return 0;

	/*
	std::cout << "[21] -> " <<  num_young(3,2);
	std::cout << "[221] -> "  <<  num_young(5,2);
	std::cout << "[3] -> "  <<  num_young(3,2*1.5+1);
	*/
	/*
	int N=3;

	std::vector<int> ytab;
	std::vector<int> frows;
	frows.push_back(2);
	frows.push_back(1);
	get_init_ytableau(N,frows,ytab);

	print_tableau(frows,ytab);

	get_next_ytableau(frows,ytab);

	print_tableau(frows,ytab);

	get_next_ytableau(frows,ytab);

	print_tableau(frows,ytab);
*/
	/*
	std::vector<int> perm={0,1,2};
	do 
	{
		std::cout << sgn(perm) <<"  ";
		for(int pk : perm)
			std::cout << pk << " ";
		std::cout << "\n";

		std::cout << perm_matrix(perm);
		std::cout << "\n";


	}while(std::next_permutation(perm.begin(),perm.end()));
	*/


	/*
	perm_a P12;
	P12.perm={1,0,2};
	P12.coeff=.5;

	//std::cout << ".5 * P12\n";
	//print_perm(P12);

	perm_a id;
	id.perm={0,1,2};
	id.coeff=-1;

	//std::cout << "Id\n";
	//print_perm(id);

	std::vector<perm_a> pvec;
	pvec.clear();
	pvec.push_back(id);


	//print_perm(perm_multiply(P12,P12));
	*/
	
}

int
test_symmetrizer()
{

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

	print_tableau(frame_rows,young_tableau);

	std::vector<perm_a> S=Srow(frame_rows,young_tableau);
	std::cout << "Symmetrizer of rows:\n";
	for(int k=0; k<S.size(); k++)
	{
		print_perm(S[k]);
	}

	
	return 0;
}

std::vector<perm_a>
Srow(std::vector<int> frame,std::vector<int> tableau)
{
	bool debug=false;
	if(debug) std::cout << "here (young::Srow)\n";
	int N=0; //number of electrons
	for(int r=0; r<frame.size(); r++)
		N=N+frame[r];
	
	//initialize S as I
	std::vector<perm_a> S;
	S.clear(); S.push_back(perm_I(N));

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

		Sr.push_back(perm_I(N));
		Sr.push_back(P);
		
		if(debug)
		{ 
			std::cout << "Sr (young::Srow)\n";
			for( perm_a p : Sr)
				print_perm(p);
			std::cout << "\n";
			std::cout << "S  (young::Srow)\n";
			for( perm_a p : S)
				print_perm(p);
			std::cout << "\n";

		}

		S=multiply(S,Sr);
		if(debug) 
		{
			std::cout << "after multiply(S,Sr) (young::Srow)\n";
			for( perm_a p : S)
				print_perm(p);
			std::cout << "\n";
		}
	}
	return S;
}
