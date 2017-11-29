// Functions for the symmetric group approach
#include<iostream>
#include<assert.h>
#include<Eigen>
#include<vector>
#include<algorithm>
#include"young.h"

using std::cout; 

//test functions
int                   test_symmetrizer();
int                   test_antisymmetrizer();

void
print_perm(perm_a P)
{
	std::cout << P.coeff << " * ";
	for(auto p : P.perm)
		std::cout << p <<" ";
	std::cout << "\n";
	return;
}

void
print_bf(basis_func bf)
{
	if(bf.coeff>0)
		std::cout << "+";
	std::cout << bf.coeff << " | ";
	for(auto p : bf.orbs)
		std::cout << p <<" ";
	std::cout << " > \n";
	return;
}

void
print_wf(std::vector<basis_func> wf, bool one_line)
{
	for(auto bf : wf)
	{
		if(bf.coeff>0)
			std::cout << "+";
		std::cout << bf.coeff << " | ";
		for(auto p : bf.orbs)
			std::cout << p <<" ";
		std::cout << "> ";
		if(!one_line) cout << "\n";
	}
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
		p_result[j]=P2.perm[P1.perm[j]];
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

perm_a
multiply(const perm_a& P1, const perm_a& P2)
{
	return perm_multiply(P1,P2);
}

std::vector<int>
perm_multiply(const  std::vector<int>& p1, const std::vector<int>& p2)
{
	bool debug =false;
	std::vector <int> ans;
	ans.clear();

	for(int j=0; j<p1.size(); j++)
	{
		ans[j]=p2[p1[j]];
		if(debug)
		{
			std::cout << "ans element " << j << " "
				  << ans[j] ;
			std::cout << " (young::perm_multiply)\n";
		}
	}



	return ans;
}

basis_func
multiply(const perm_a P,const basis_func v)
{
	bool debug=false;
	using std::cout;
	if(P.perm.size()!=v.orbs.size())
	{
		std::cout << "Error in young::multiply(perm_a, basis_func), "
			  << "incompatible sizes.\n";
		return v;
	}
	basis_func pv;
	pv.orbs.resize(v.orbs.size());

	for(int j=0; j<v.orbs.size(); j++)
	{
		pv.orbs[j]=v.orbs[P.perm[j]];
	}

	pv.coeff=P.coeff*v.coeff;
	if(debug) std::cout << "pv.coeff = P.coeff * v.coeff " << pv.coeff << " = " << P.coeff << " * " << v.coeff << "\n";
	return pv;

}

bool
compare_perms(perm_a f1, perm_a f2)
{
	int n1=f1.perm.size();
	int n2=f2.perm.size();
	if(n1<n2)
		return n1<n2;
	if(n2>n1)
		return n1<n2;

	//from now on we'll assume n 
	int nelems=n1;

	for(int i=0; i<nelems; i++)
	{
		if(f1.perm[i] == f2.perm[i])
			continue;

		//now assume they're unequal and return
		return f1.perm[i]<f2.perm[i];
	}
	// in the case that they're all equal we'll end up here
	return true;

}

bool
compare_bfs(basis_func f1, basis_func f2)
{
	int n1=f1.orbs.size();
	int n2=f2.orbs.size();
	if(n1<n2)
		return n1<n2;
	if(n2>n1)
		return n1<n2;

	//from now on we'll assume n 
	int norbs=n1;

	for(int i=0; i<norbs; i++)
	{
		if(f1.orbs[i] == f2.orbs[i])
			continue;

		//now assume they're unequal and return
		return f1.orbs[i]<f2.orbs[i];
	}
	// in the case that they're all equal we'll end up here
	return true;

}

std::vector<perm_a>
multiply(std::vector<perm_a> A, std::vector<perm_a> B)
{
	bool debug=false;

	//multiply and sort a list of permutations and a list of basis functions
	std::vector<perm_a> AB;
	for(perm_a p: A)
		for( perm_a q : B)
		{
			//insert it into the list
			AB.push_back(multiply(p,q));
		}
        //first sort 
	std::sort(AB.begin(),AB.end(),compare_perms);

	if(debug)
	{
		std::cout << "(young::multiply(vector<perm_a>, vector<perm_a>))\n";
		for(auto p: AB)
		{
			print_perm(p);
			std::cout << "\n";
		}
	}
	

	//then check for repetitions
	int like_terms;
	int terms;
	do
	{
		like_terms=0;
		terms=AB.size();
		for(int j=0; j<terms-1; j++)
		{
			//the size of AB changes dynamically 
			//so we need to re-evaluate if we can execute the next line
			if( j+1 >= AB.size() )
				break;

			if(AB[j].perm == AB[j+1].perm || fabs(AB[j+1].coeff)<1e-10)
			{
				AB[j].coeff=AB[j].coeff+AB[j+1].coeff;
				AB.erase(AB.begin()+j+1);
				like_terms++;
			}


		}

	}while(like_terms!=0);
	
	if(debug)
	{
		std::cout << "(young::multiply(vector<perm_a>, vector<perm_a>))\n";
		for(auto p: AB)
		{
			print_perm(p);
			std::cout << "\n";
		}
	}
		
	return AB;

}

std::vector<basis_func>
multiply(std::vector<perm_a> O, std::vector<basis_func> wf)
{
	bool debug=false;

	//multiply and sort a list of permutations and a list of basis functions
	std::vector<basis_func> new_wf;
	for(perm_a p: O)
		for( basis_func bf : wf)
		{
			//insert it into the list
			new_wf.push_back(multiply(p,bf));
		}
	if(debug)
	{
		std::cout << "(young::multiply(vector<perm_a>, vector<basis_func>))\n";
		for(auto bf: new_wf)
		{
			print_bf(bf);
			std::cout << "\n";
		}
	}
	

        //first sort 
	std::sort(new_wf.begin(),new_wf.end(),compare_bfs);

	//then check for repetitions
	int like_terms;
	int terms;
	do
	{
		like_terms=0;
		terms=new_wf.size();
		for(int j=0; j<terms-1; j++)
		{
			//the size of new_wf changes dynamically 
			//so we need to re-evaluate if we can execute the next line
			if( j+1 >= new_wf.size() )
				break;

			if(new_wf[j].orbs == new_wf[j+1].orbs)
			{
				new_wf[j].coeff=new_wf[j].coeff+new_wf[j+1].coeff;
				new_wf.erase(new_wf.begin()+j+1);
				like_terms++;
			}

		}

	}while(like_terms!=0);
	
	if(debug)
	{
		std::cout << "(young::multiply(vector<perm_a>, vector<basis_func>))\n";
		for(auto bf: new_wf)
		{
			print_bf(bf);
			std::cout << "\n";
		}
	}
		
	return new_wf;

}

std::vector<basis_func>
multiply(std::vector<perm_a> O, basis_func f0)
{
	bool debug=false;

	//boot strap
	
	std::vector<basis_func> f={f0};

	std::vector<basis_func> Of = multiply(O,f);
	if(debug)
	{
		std::cout << "(young::multiply(vector<perm_a>, basis_func))\n";
		for(auto bf: Of)
		{
			print_bf(bf);
			std::cout << "\n";
		}
	}
	return Of;
}

std::vector<double>
multiply(perm_a P, std::vector<double> v)
{
	if(P.perm.size() != v.size())
	{
		std::cout << "Error in young::multiply(perm_a, vector<int>), "
			  << "incompatible sizes.\n";
		return v;
	}

	std::vector<double> pv;
	pv.resize(v.size());

	for(int j=0; j<P.perm.size(); j++)
	{
		pv[j] = v[P.perm[j]]; // don't use coefficient here. Just keep it as
				      // a rearranged list of orbitals; 
	}
	return pv;
}

std::vector<basis_func>
add(const std::vector<basis_func> A, const std::vector<basis_func> B)
{
	bool debug=true;
	std::vector<basis_func> C;

	for(auto a : A)
		C.push_back(a);
	for(auto b : B)
		C.push_back(b);

	return multiply({perm_I(A[0].orbs.size())},C);
	
}

std::vector<basis_func>
multiply(const double a, const std::vector<basis_func> F)
{
	std::vector<basis_func> G;
	basis_func bf;
	for(int k=0; k<F.size(); k++)
	{
		bf=F[k];
		bf.coeff=bf.coeff*a;
		
		G.push_back(bf);

	}

	return G;
}

std::vector<basis_func>
normalize(const std::vector<basis_func> wf)
{
	double norm2;
	for(int k=0; k<wf.size(); k++)
		norm2+=(wf[k].coeff)*std::conj(wf[k].coeff);

	return multiply(1.0/sqrt(norm2),wf);
}


//
/*  Here we don't check for like terms, better version is implemented
//
//

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
*/

/*
int
pos_ymax(const int pos, const int M, const std::vector<int> frame, const std::vector<int> tableau)
{
	int max_val=M-1;
	return max_val;
}
*/

int
pos_ymin(const int pos,const std::vector<int> frame, std::vector<int> tableau)
{

	int above_pos;
	int min_idx=0;
	int row,col;
	bool debug=false;

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

std::vector<perm_a> 
Ey(std::vector<int> frame_rows,std::vector<int> yT)
{
	std::vector<perm_a> S= Srow(frame_rows,yT);
	std::vector<perm_a> A= Acol(frame_rows,yT);
	return multiply(S,A);
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
get_next_young_tableau(const std::vector<int> frame, std::vector<int>& tableau)
{
	int ctr=0;
	while(1)
	{
	
		std::next_permutation(tableau.begin(),tableau.end());
		if(isvalid_young(frame,tableau))
			return 0;

		ctr++;
		if(ctr>1000)	
			break;
	}
	std::cout<<"warning: over 1000 permutations rejected (get_next_young_tableau)\n";
	return 1;

}

std::vector<int>
invperm(const std::vector<int> T)
{
    //make invT
    std::vector<int> invT;
    invT.clear();
    for (int i = 0 ; i < T.size() ; i++) 
	    invT.push_back(i);

    sort(invT.begin(),invT.end(),[&](const int& a, const int& b)
             {
	     return(T[a] < T[b]);
 	   });

    return invT;

}

int
test_young_tableau_generation()
{
	//test the computation of irreducible representation dimensions
	std::cout << "[21] -> " <<  num_young(3,2);
	std::cout << "[221] -> "  <<  num_young(5,2);
	std::cout << "[3] -> "  <<  num_young(3,2*1.5+1);
	
	//test the generation of young tableau in two ways for [221] young frame.
	std::vector<int> Yrows;
	Yrows.push_back(2);
	Yrows.push_back(2);
	Yrows.push_back(1);

	// 1 2 
	// 3 4
	// 5
	std::vector<int> ytab={0,1,2,3,4};

	std::cout << "[221] -> N=" <<  num_young(5,2) << "\n";

	for(int i=0; i<num_young(5,2); i++)
	{
		std::cout << i << ": ";
		for(auto item : ytab)
			std::cout << item << " ";
		std::cout << "\n";
		get_next_young_tableau(Yrows,ytab);

	}

	std::cout << "\n";
	int ctr=0;

	for(int k=0; k<200; k++)
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

		if(next_permutation(ytab.begin(),ytab.end()))
			continue;
		else
			break;

	}


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

std::vector<perm_a>
Acol(std::vector<int> frame,std::vector<int> tableau)
{
	bool debug=false;
	if(debug) std::cout << "here (young::Acol)\n";
	int N=0; //number of electrons
	for(int r=0; r<frame.size(); r++)
		N=N+frame[r];
	
	//initialize S as I
	std::vector<perm_a> A;
	A.clear(); A.push_back(perm_I(N));

	//get col information
	int n_col[2];
	n_col[0]=frame.size();
	n_col[1]=N-n_col[0];
	//total number of boxes = N = n_column_1 + n_column_2

	std::vector<perm_a> Ac;
	for(int c=0; c<2; c++)
	{

		Ac.clear();

		if(n_col[c]<2) 
			// nothing to do if the col length is 1 or 0
			continue;

		//get column
		std::vector<int> col_perm;
		col_perm.clear();
		for(int row=0; row < n_col[c]; row++)
		{
			col_perm.push_back(tableau[tableau_pos(row,c,frame)]);
		}
		std::vector<int> col_sorted=col_perm;
		std::sort (col_sorted.begin(),col_sorted.end());

		//make permutations over the column, keep items not in row in place
		perm_a P;
		P.perm.clear();
		int ctr=0;
		do
		{
			P.perm.clear();
			ctr=0;
			for(int k=0; k< N; k++)
			{
				if(col_sorted[ctr]==k)
				{
					//pull from permutation
					
					P.perm.push_back(col_perm[ctr]);
					//advance ctr
					ctr++;
				}
				else    //else just insert in place
					P.perm.push_back(k);
			}

			P.coeff=sgn(P.perm);

			Ac.push_back(P);

		}while( std::next_permutation(col_perm.begin(), col_perm.end()) );
		
		if(debug)
		{ 
			std::cout << "Ac (young::Acol)\n";
			for( perm_a p : Ac)
				print_perm(p);
			std::cout << "\n";
			std::cout << "A  (young::Acol)\n";
			for( perm_a p : A)
				print_perm(p);
			std::cout << "\n";

		}

		A=multiply(A,Ac);
		if(debug) 
		{
			std::cout << "after multiply(A,Ac) (young::Acol)\n";
			for( perm_a p : A)
				print_perm(p);
			std::cout << "\n";
		}
		
	}
	return A;
}

int
test_antisymmetrizer()
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

	std::cout << "frame=[221], tab=[0 1 /2 3 /4]\n";
	std::vector<perm_a> A=Acol(frame_rows,young_tableau);
	std::cout << "Anti-symmetrizer of cols:\n";
	for(int k=0; k<A.size(); k++)
	{
		print_perm(A[k]);
	}
	// --------------------------------------------------------------------
	frame_rows.clear();
	frame_rows.push_back(2);
	frame_rows.push_back(1);

	young_tableau.clear();
	young_tableau.push_back(0);
	young_tableau.push_back(1);
	young_tableau.push_back(2);

	print_tableau(frame_rows,young_tableau);

	std::cout << "frame=[21], tab=[01/2]\n";
	A=Acol(frame_rows,young_tableau);
	std::cout << "Anti-symmetrizer of cols:\n";
	for(int k=0; k<A.size(); k++)
	{
		print_perm(A[k]);
	}
	
	return 0;
}

int
test_invperm()
{
	perm_a P;
	P.perm={1,3,0,4,2};

	perm_a invP;

	invP.perm=invperm(P.perm);

	
	std::cout << "Permutation:\n";
	for(int p : P.perm)
		std::cout << p << " ";
	std::cout << "\nInverse Permutation:\n";
	for(int p : invP.perm)
		std::cout << p << " ";
	std::cout << "\n";

        auto ans=perm_multiply(P,invP);
	std::cout << "Product:\n";
	for(int p : ans.perm)
		std::cout << p << " ";
	std::cout << "\n";


	return 0;
}

double
dot(std::vector<basis_func> L,std::vector<basis_func> R)
{
	//this function can probably be improved to O(N) 
	double dot=0;
	for( basis_func a : L)
		for( basis_func b : R)
		{
			if(a.orbs==b.orbs)
				dot=dot+a.coeff*b.coeff;
		}

	return dot;
}

void
test_dot()
{
	using std::vector;

	vector<basis_func> wf;
	basis_func f;
	wf.clear();

	f.coeff= +1/sqrt(3);f.orbs.clear();
	f.orbs.push_back(1);f.orbs.push_back(0);f.orbs.push_back(0);
	wf.push_back(f);

	f.coeff= -2/sqrt(3);f.orbs.clear();
	f.orbs.push_back(0);f.orbs.push_back(1);f.orbs.push_back(0);
	wf.push_back(f);

	f.coeff= +1/sqrt(3);f.orbs.clear();
	f.orbs.push_back(0);f.orbs.push_back(0);f.orbs.push_back(1);
	wf.push_back(f);

	std::cout << " state 1 : \n"; // somewhat random initial state
	for(auto bf: wf)
		print_bf(bf);
	std::cout << "Norm^2 : " << dot(wf,wf) << "\n"; 

	vector<basis_func> wf2;
	wf2.clear();

	//how to pick the initial state to make sure it has projection on all irreps?
	f.coeff= +2/sqrt(3);f.orbs.clear();
	f.orbs.push_back(1);f.orbs.push_back(0);f.orbs.push_back(0);
	wf2.push_back(f);

	f.coeff= -2/sqrt(3);f.orbs.clear();
	f.orbs.push_back(0);f.orbs.push_back(1);f.orbs.push_back(0);
	wf2.push_back(f);

	f.coeff= +1/sqrt(3);f.orbs.clear();
	f.orbs.push_back(0);f.orbs.push_back(0);f.orbs.push_back(1);
	wf2.push_back(f);

	std::cout << " state 2 : \n"; // somewhat random initial state
	for(auto bf: wf2)
		print_bf(bf);
	std::cout << "Norm^2 : " << dot(wf2,wf2) << "\n"; 

	std::cout << " state 1 . state 2 : " << dot(wf,wf2) << "\n";
	std::cout << " (2/3) + (4/3) + 1/3 = " << 7/3. << "\n";

	
	return;
}

std::vector<basis_func> 
initial_state(const int N, const std::vector<int> occ_orbs)
{
	using std::vector;
	bool debug=false;

	if(debug)
	{
		cout << "in initial state...\nOccupied orbitals: ";
		for(auto o : occ_orbs)
			cout << o << " ";
		cout << "\n";
	}

	vector<basis_func> wf;
	wf.clear();

	std::vector<int> primitive;
	basis_func f;

	//sort input occupied orbital list
	primitive=occ_orbs;
	std::sort(primitive.begin(),primitive.end());



	/* RANDOM STATE */
	const std::vector<double> rand_nums={0.0402621714989,0.720446771216,0.113868284377,0.314375237251,0.368092160337,0.412386548763,0.0495483892997,0.814967948581,0.521897624288,0.487846464935,0.299025914147,0.455005750774};
	
	double norm2;

	for(int m=0; m<N; m++)
		norm2+=rand_nums[m]*rand_nums[m];



	for(int m=0; m<rand_nums.size(); m++)
	{


		f.coeff=rand_nums[m]/sqrt(norm2);
		f.orbs=primitive;

		wf.push_back(f);

		//advance primitive list, should we be looping over the number of permutations?
		if(!std::next_permutation(primitive.begin(),primitive.end()))
		{
			if(debug)
				cout << "quitting loop with m="<< m << "\n";
			break;
		}

	}

	if(debug)
	{
		cout << "wf= ";
	       	for( auto v : wf) print_bf(v);	
		cout << "\n\n";
	}
	return wf;

	/*
	a=a/sqrt(norm2);
	b=b/sqrt(norm2);
	c=c/sqrt(norm2);

	f.coeff= a;
	f.orbs.clear();
	f.orbs.push_back(1);f.orbs.push_back(0);f.orbs.push_back(0);
	wf.push_back(f);

	f.coeff= b;
	f.orbs.clear();
	f.orbs.push_back(0);f.orbs.push_back(1);f.orbs.push_back(0);
	wf.push_back(f);

	f.coeff= c;
	f.orbs.clear();
	f.orbs.push_back(0);f.orbs.push_back(0);f.orbs.push_back(1);
	wf.push_back(f);
	*/


	/* 100 INVARIANT UNDER P_{T_0 <- T_j} = P23
	 
	f.coeff= 1;
	f.orbs.clear();
	f.orbs.push_back(1);f.orbs.push_back(0);f.orbs.push_back(0);
	wf.push_back(f);

	 */

	/* 010 
	f.coeff= 1;
	f.orbs.clear();
	f.orbs.push_back(0);f.orbs.push_back(1);f.orbs.push_back(0);
	wf.push_back(f);
	 */

	/* 001  */
	f.coeff= 1;
	f.orbs.clear();
	f.orbs.push_back(0);f.orbs.push_back(0);f.orbs.push_back(1);
	wf.push_back(f);

	/* 100 + 010 + 001 
	f.coeff=1/sqrt(3);
	f.orbs.clear();f.orbs.push_back(0);f.orbs.push_back(0);f.orbs.push_back(1);
	wf.push_back(f);

	f.coeff=1/sqrt(3);
	f.orbs.clear();f.orbs.push_back(0);f.orbs.push_back(1);f.orbs.push_back(0);
	wf.push_back(f);

	f.coeff=1/sqrt(3);
	f.orbs.clear();f.orbs.push_back(1);f.orbs.push_back(0);f.orbs.push_back(0);
	wf.push_back(f);
	*/

	return wf;
}

/*
 Get the basis for irrep labelled by frame F projected from a given state
 */
std::vector<std::vector<basis_func>> 
get_irrep_basis(std::vector<int> F, std::vector<basis_func> wf ) 
{
	using std::vector;

	//get basis
	bool debug=false;
	if(debug)
	{
		cout << "\nIn young:get_irrep_basis(F= ";
		for( auto f : F) cout << f << " ";
		cout << ", wf= ";
	       	for( auto v : wf) print_bf(v);	
		cout << ")\n\n ";
	}

	
	Matrix S;
	int N=0;		//number of electrons 
	int num_ytabs;		//number of young tableaux
	int unpaired=0;		//unpaired boxes
	double s;		//spin quantum number
	int ms;			//degeneracy, 2s+1
	std::vector<int> perm;  //integer permutation variable
	perm_a invpT;		//permutation algebraic options
	int rank;		//rank of the overlap matrix

	for(int k=0; k<F.size(); k++)
	{
		//count number of boxes
		N=N+F[k];
		if(F[k]==1)
			unpaired++;
	}

	s  = unpaired*.5;
	ms = 2*s+1;
	num_ytabs=num_young(N,ms);

	//identity perm
	perm.clear();
	for(int i=0; i<N; i++)
		perm.push_back(i);

	//pull young operator
	auto E0= Ey(F,perm);

	if(debug)
	{
		std::cout << "E0:\n"; 
		for(auto p : E0)
			print_perm(p);
		std::cout << "\nwf:\n"; 
		for( basis_func bf : wf)
		 	print_bf(bf);
		std::cout << "\nE0*wf\n";
		for( basis_func bf : multiply(E0,wf))
		 	print_bf(bf);
		cout << "\n";


	}

	vector<vector<basis_func>> C;
	C.clear();

	//std::cout << "num+ytabs " << num_ytabs << "\n";

	int  ctr=0;
	do{
		invpT.perm=invperm(perm);
		C.push_back(multiply(multiply({invpT},E0),wf));

		if(debug)
		{
			std::cout << "iteration: "<< ctr++ << ", perm: ";
			print_perm(invpT);


			std::cout << "C matrix: \n";

			for(int j=0; j<C.size(); j++)
			{
				std::cout << j << ": \n";
				for(int k=0; k<C[j].size(); k++)
					print_bf(C[j][k]);
			}
			std::cout << "--\n\n";
		}


		if(C.size() ==1 )
			continue;

		S.resize(C.size(),C.size());
		for(int i=0; i<S.cols(); i++)
			for(int j=0; j<S.rows(); j++)
				S(i,j)=dot(C[i],C[j]);

		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(S);

		if (eigensolver.info() != Eigen::Success) 
		{
			std::cout << "Failed to solve eigensystem\n";
			abort();
		}

		Matrix s_evecs=eigensolver.eigenvectors();
		Matrix s_evals=eigensolver.eigenvalues();

		if(debug)
		{
			std::cout << "S: \n" << S << "\n";
			cout << "Eigen system:\n";
			std::cout << "S eigenvalues: \n" << s_evals << "\n";
			//std::cout << "S eigenvectors: \n" << s_evecs << "\n";
		}


		rank=0;
		for(int j=0; j<s_evals.rows(); j++)
			if( s_evals(j) < 1e-5 )
			{
				if(debug)
				{
					std::cout << "rejected  eigenvalue: " << s_evals(j) <<  " , j=" << j << "\n"; 

					for(int k=0; k<C[j].size(); k++)
						print_bf(C[j][k]);
					std::cout << "\n";
				}
				C.pop_back();

				break;
			}
			else
			{
				if(debug)
				{
					std::cout << "accepted j: " << j; 
					for(int k=0; k<C[j].size(); k++)
						print_bf(C[j][k]);
					std::cout << "\n";
				}
				rank=rank+1;
			}

		if(rank==num_ytabs)
			break;

	}while( std::next_permutation(perm.begin(),perm.end()) );


	if(rank<num_ytabs)
		cout << "Warning: input state does not span input irrep. (young::get_irrep_basis)\n";

	if(debug)
		std::cout << S << "\n";

	if(debug) std::cout << "Size of C: " << C.size() << "\n";

	return C;
}

//here we're using symmetric orthogonalization
Matrix 
get_ortho_matrix_X(std::vector<std::vector<basis_func>> C)
{
	bool debug=false;

	//overlap matrix
	Matrix S;
	S.resize(C.size(),C.size());

	for(int i=0; i<S.cols(); i++)
		for(int j=0; j<S.rows(); j++)
			S(i,j)=dot(C[i],C[j]);

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(S);

	if (eigensolver.info() != Eigen::Success) 
	{
		std::cout << "Failed to solve eigensystem\n";
		abort();
	}

	Matrix s=eigensolver.eigenvalues();
	Matrix U=eigensolver.eigenvectors();

	Matrix x=s;
	Matrix invs=s;

	for(int j=0; j<x.rows(); j++)
	{
		x(j)=1/sqrt(s(j));
		invs(j)=1/(s(j));
	}

	Matrix invS=U*invs.asDiagonal()*U.adjoint();

	if(debug)
	{
		std::cout << "S*inverse(S): \n"
		     << S*invS << "\n";

		std::cout << "s: ";
		for( int j=0; j<x.rows(); j++)
			std::cout << s(j) << " ";
		std::cout << "\n";

		std::cout << "x: ";
		for( int j=0; j<x.rows(); j++)
			std::cout << x(j) << " ";
		std::cout << "\n";
	}

	Matrix X=x.asDiagonal(); //*U.adjoint(); 

	if(debug)
	{
		std::cout << "X=[" << X.rows() << "," <<  X.cols() << "]\n";
		std::cout << "U=[" << U.rows() << "," <<  U.cols() << "]\n";
	}

	return U*x.asDiagonal()*U.adjoint();
}

Matrix
invert_matrix(Matrix X)
{
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(X);

	Matrix evals=eigensolver.eigenvalues();
	Matrix U=eigensolver.eigenvectors();
	Matrix invx=evals;

	for(int j=0; j<X.rows(); j++)
	{
		invx(j)=1/(evals(j));
	}

	Matrix invX= U * invx.asDiagonal() * U.adjoint();

	return invX;
}

std::vector<perm_a>
wigner_op(int i, int j,int Nelec, int gS, std::vector<std::vector<basis_func>> C)
{
	bool debug=false;
	using std::vector;
	using std::cout;

	int dl=num_young(Nelec,gS);
	int h = fac[Nelec];

	perm_a P;
	vector<perm_a> Wij;

	Matrix S;
	S.resize(C.size(),C.size());
	for(int i=0; i<S.cols(); i++)
		for(int j=0; j<S.rows(); j++)
			S(i,j)=dot(C[i],C[j]);

	Matrix invS=invert_matrix(S);

	P.coeff=1;
	P.perm={0,1,2};


	if(debug)
	{
		cout << "dl: " << dl << "\n";
		cout << "h: " << h << "\n";
		cout << "dl/h: " << (dl*1.0)/(1.0*h) << "\n";
		cout << "len(C): " << C.size() << "\n"; 
	}

	perm_a summand;

	do //loop over permutations
	{
		if(debug)
		{
			cout <<"Permutation: ";
			print_perm(P);
		}

		vector<vector<basis_func>> out=C;
		for(int i=0; i<C.size(); i++)
			out[i]=multiply({P},C[i]);

		//take dot product with each basis function to get D
		//i.e. C^+ (PC)

		Matrix out2=S; //just because S has the correct dimensions
		for(int i=0; i<out2.cols(); i++)
			for(int j=0; j<out2.rows(); j++)
				out2(i,j)=dot(C[i],out[j]);

		if(debug)
			cout << "C^+ (PC) \n" << out2 << "\n";

		Matrix D_P = invS*out2;

		if(debug)
			cout << "S^-1 C^+ (PC) \n" << D_P << "\n";

		//cout << "invP: \n";
		
	       	summand.perm=invperm(P.perm);

		// Wigner operator $W_{i,j}$ is proportional to $DP_{j,i}$, not
		// $DP_{i,j}$. This is a result of using the orthogonality 
		// theorem. That appears in W_{i,j} should be
		//
		// See elementary group theory notes equations (32-34)  
		summand.coeff= (dl*1.0)*D_P(j,i)/(1.0*h);

		if(fabs(summand.coeff)>1e-12)
			Wij.push_back(summand);

	}while( std::next_permutation(P.perm.begin() , P.perm.end()) );

	return Wij;
}

int 
test_wigner(std::vector<int> F, int ms)
{	
	using std::vector;
	using std::cout;

	/*
	std::vector<int> F;
	F.clear();
	F.push_back(2);
	F.push_back(1);
	*/

	
	std::cout << "Frame: ";
	for(int j=0; j<F.size(); j++)
		cout << F[j] << " " ;
	cout << "\n";

	std::cout << "spin multiplicity: " << ms << "\n";



	//std::vector<int> T;
	//T.clear();
	//T.push_back(0);
	//T.push_back(1);
	//T.push_back(2);

	int N=0;
	for(int k=0; k<F.size(); k++)
		//count number of boxes
		N=N+F[k];


	int num_ytabs=num_young(N,ms);

	std::cout << "number of young tableau: " << num_ytabs << "\n";

	//random state with orbital 0 occupied twice and orbital 1 occupied once
	vector<int> occ={0,0,1};
	vector<basis_func> wf = initial_state(N,occ);

	std::cout << "random state (fixed): \n"; 
	for(auto bf : wf)
		print_bf(bf);
	std::cout << "Norm^2 : " << dot(wf,wf) << "\n"; 

	auto C= get_irrep_basis(F,wf);

	assert(C.size()==num_ytabs);

	cout << "C.size():" << C.size() << "\n";

	auto X= get_ortho_matrix_X(C);
	
	Matrix S;
	S.resize(C.size(),C.size());
	for(int i=0; i<S.cols(); i++)
		for(int j=0; j<S.rows(); j++)
			S(i,j)=dot(C[i],C[j]);

	perm_a P;
	P.coeff=1;
	P.perm={0,1,2};




	int i=0,j=0;
	auto Wij=wigner_op(i,j,N,ms,C);

	std::cout << "\ni = " << i << ", j = " << j << "\n--\nW_{i,j}:\n";
	for( auto p : Wij )
		print_perm(p);
	cout << "\n";

	vector<perm_a> Wij2=multiply(Wij,Wij);

	/*
	cout << "W_{i,j}^2: \n";
	for(auto p : Wij2)
		print_perm(p);
	cout << "\n";
	*/

	i=1,j=0;
	Wij=wigner_op(i,j,N,ms,C);

	std::cout << "\ni = " << i << ", j = " << j << "\n--\nW_{i,j}:\n";
	for( auto p : Wij )
		print_perm(p);
	cout << "\n";

	Wij2=multiply(Wij,Wij);

	/*
	cout << "W_{i,j}^2: \n";
	for(auto p : Wij2)
		print_perm(p);
	cout << "\n";
	*/

	auto Wji=wigner_op(j,i,N,ms,C);

	std::cout << "\nj = " << j << ", i = " << i << "\n--\nW_{j,i}:\n";
	for( auto p : Wji )
		print_perm(p);
	cout << "\n";

	auto WijWji=multiply(Wij,Wji);
	auto WjiWij=multiply(Wji,Wij);

	/*
	cout << "W_{j,i}* W_{i,j}: \n";
	for(auto p : WjiWij)
		print_perm(p);
	cout << "\n";

	cout << "W_{i,j}* W_{j,i}: \n";
	for(auto p : WijWji)
		print_perm(p);
	cout << "\n";
	*/

	cout << "i = " << i << "\n--\nW_{i,i}:\n";
	for(auto p: wigner_op(i,i,N,ms,C))
		print_perm(p);
	cout << "\n";
	


	return 0;
}

int
alt_young_main()
{
	std::vector<int> F;
	F.clear();
	F.push_back(2);
	F.push_back(1);
	
	using std::vector; 

	std::vector<basis_func> wf= initial_state(3,{0,0,1});
	std::vector<std::vector<basis_func>> C = get_irrep_basis(F,wf);
	Matrix X= get_ortho_matrix_X(C);

	cout << "Basis functions for [1,1; 2]:\n"; 

	// G_j = sum_m  g_m X_{mj}  
	vector<basis_func> G0=add( multiply(X(0,0),C[0]) , multiply(X(1,0),C[1]) );
	vector<basis_func> G1=add( multiply(X(0,1),C[0]) , multiply(X(1,1),C[1]) );
	
	cout << "G0 : \n";
        for( auto f: G0 )
		print_bf(f);

	cout << "G1 : \n";
        for( auto f: G1 )
		print_bf(f);
	
	vector<perm_a> W00=wigner_op(0,0,3,2,C);
	vector<perm_a> W10=wigner_op(1,0,3,2,C);

	vector<basis_func> F0=multiply(W00,C[0]);
	vector<basis_func> F1=multiply(W10,C[0]);

	cout << "F0 : \n";
        for( auto f: F0 )
		print_bf(f);

	cout << "normalized F0 : \n";
	for(auto f: normalize(F0))
		print_bf(f);

	cout << "F1 : \n";
        for( auto f: F1 )
		print_bf(f);

	cout << "Basis functions for [1,1; 3]:\n"; 
	wf= initial_state(3,{0,0,2});
	C=  get_irrep_basis(F,wf);
	X= get_ortho_matrix_X(C);

	G0=add( multiply(X(0,0),C[0]) , multiply(X(1,0),C[1]) );
	G1=add( multiply(X(0,1),C[0]) , multiply(X(1,1),C[1]) );

	cout << "G0 : \n";
        for( auto f: G0 )
		print_bf(f);

	cout << "G1 : \n";
        for( auto f: G1 )
		print_bf(f);
	
	W00=wigner_op(0,0,3,2,C);
	W10=wigner_op(1,0,3,2,C);

	F0=multiply(W00,C[0]);
	F1=multiply(W10,C[0]);

	cout << "F0 : \n";
        for( auto f: F0 )
		print_bf(f);

	cout << "F1 : \n";
        for( auto f: F1 )
		print_bf(f);

	return 0;
}

