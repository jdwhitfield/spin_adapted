#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<assert.h>
#include<iomanip>
#include<random>
#include<chrono>
#include<algorithm>
//
/// fac[k] = k!
static constexpr std::array<int64_t,21> fac = {{1L, 1L, 2L, 6L, 24L, 120L, 720L, 5040L, 40320L, 362880L, 3628800L, 39916800L,
                                                    479001600L, 6227020800L, 87178291200L, 1307674368000L, 20922789888000L,
                                                    355687428096000L, 6402373705728000L, 121645100408832000L,
                                                    2432902008176640000L}};
 

int
PrimDiffs(int M, const std::vector<int>& primRef, const std::vector<int>& primK,
	  std::vector<int>&         Pr  , std::vector<int>&         Pk, 
	  double scaling_factor);

int 
perm_multiply(const std::vector<int>& P1, const std::vector<int>& P2,          
              std::vector<int>&      result);

int
nchoosek(int M,int N)
{
	return fac[M]/(fac[M-N]*fac[N]);
}

int
num_weyl(int M, int N, int ms)
{
	// Hall-Robinson formula
	// See slide 21 of generalized_operators.pdf
	
	if(M>21) 
	{
		std::cout << "M is too large (ci_matrix::num_weyl)\n"; return -999;
	}

	double S=(ms -1)/2.0;
	/*
	std::cout << "M: " << M << "\n";
	std::cout << "N: " << N << "\n";
	std::cout << "S: " << S << "\n";
	*/ 
	
        auto a=(ms)*1.0/(M+1.0);
        /*
	std::cout << "2S+1/M+1="<<ms <<"/"<<M+1.0 << ": "<< a << "\n ";
        */

	auto k1= .5*N+S+1;

        auto b=nchoosek(M+1,k1);
	/*
	std::cout << "nchoosek(M+1,N/2+S+1)=nchoosek(" 
		  << M+1 << "," << k1 << "): "<<b << "\n ";
        */
        
	auto k2= .5*N-S;
        auto c=nchoosek(M+1,k2);
	/*
	std::cout << "nchoosek(M+1,N/2-S)=nchoosek(" 
		  << M+1 << "," << k1 << "): "<<c << "\n ";
         */

        //int64_t a=(ms/(M+1.0))*fac[M]/fac[M-1-(N-2.0*S)/2.0];
        //int64_t b=fac[M+1]/fac[M+S-N/2.0];
        //int64_t c=fac[1+(N+2.0*S)/2.0]*fac[(N/2.0)-S];

return  (a*b*c);	
}


int 
antisymmeterizer(const std::vector<int>& idx_list, std::vector<int>& input_string)
{

	/*
	 123
	 132
	 213
	 231
	 312
	 321
	 */ 

	int n=idx_list.size();
	for(int j=0; j<fac[n]; j++)
	{
	   std::vector<int> permj(n);
	}


	return 0;
}

class tableaux
{
	public: 
		std::vector<int> shape;
};

void
print(const std::vector<int> tableau, const std::vector<int> frame)
{
	//print tableau
	for(auto c : tableau)
	   std::cout << c <<" ";
	std::cout << std::endl;
}

void
pos_to_row_col(const int pos,const std::vector<int> frame, int& row, int& col)
{
	int ctr=0;
	int val=frame[0];


	/*
	for(int s=0; s<frame.size(); s++)
	{

		int partial_s=0;
		for(int idx=0; idx<1+s; idx++)
			partial_s=partial_s+frame[idx];

		//std::cout << s << ": " <<  partial_s << "\n";

		if( partial_s-1 <pos) // -1 because we count pos from zero
			continue;

		//row=k;
		if(partial_s-1==pos) col=0;
		

		
	}
	*/

	while(1)
	{
		
		//val has partial sum of the frame rows
		//e.g. 221
		//val = 2
		//val +=2; val=4;
		//val +=1; val=5;
		
		//ctr contains the row number index
		
		//if pos is less than val then it's in this row
		if(pos<val)
		{
		   row=ctr;

		   if(frame[row]==1)
			   col=0;

		   if(frame[row]==2)
		   {
			   //std::cout << val << " and " << 2*row << "\n";

			   /*
			   int temp=0;
			   for(int m=0; m<row; m++)
			   {
				   temp=temp+frame[m];
			   }
			   */

			   
			   col=pos % 2; //if its odd its the second elem in the row
		   }

		   return;

		   break;
		}
		
		ctr++;

		if(ctr==frame.size())//out of bounds
		{
			std::cout << "out of bounds ( ::pos_to_row)\n";
			std::cout << "pos: "<< pos << " , sum_frame : ";
			int tot=0;
			for(int m=0; m<frame.size(); m++)
				tot=tot+frame[m];
			std::cout << tot << "\n";
			return;
			break;
		}

		//add in next row to total val
		val=val+frame[ctr];
	}

}

int 
main()
{
	using namespace std;

	int M=4;
	int N=5;
	int ctr=0;

	std::vector<int> frame{2,2,2,1,1};
	std::vector<int> tableau;

	//initialize tableau
	ctr=0;
	for(int row=0; row<frame.size(); row++)
	{
		if(frame[row]==2)
		{
			tableau.push_back(ctr);
			tableau.push_back(ctr);
			ctr++;
		}
		else
		{
			tableau.push_back(ctr);
			ctr++;
		}
	}

	print(tableau,frame);

	//int row,col;
	int row,col;
	pos_to_row_col(0,frame,row,col);
	cout<< "position 0 is in row " << row << ", col: "<< col << endl;
	pos_to_row_col(1,frame,row,col);
	cout<< "position 1 is in row " << row << ", col: "<< col << endl;
	pos_to_row_col(2,frame,row,col);
	cout<< "position 2 is in row " << row << ", col: "<< col << endl;
	pos_to_row_col(3,frame,row,col);
	cout<< "position 3 is in row " << row << ", col: "<< col << endl;
	pos_to_row_col(4,frame,row,col);
	cout<< "position 4 is in row " << row << ", col: "<< col << endl;
	pos_to_row_col(5,frame,row,col);
	cout<< "position 5 is in row " << row << ", col: "<< col << endl;
	pos_to_row_col(6,frame,row,col);
	cout<< "position 6 is in row " << row << ", col: "<< col << endl;
	pos_to_row_col(7,frame,row,col);
	cout<< "position 7 is in row " << row << ", col: "<< col << endl;


	return 0;

	/*
	int row=0;
	int col=0;
	*/
	int pos=N-1;
	ctr=0;
	bool reset=false;
	while(1)
	{
		//advance clock digit
		tableau[pos]++;

		//check if it went too far
		if(tableau[pos]==M)
		{

			//reset lower digits
			for(int p=pos+1; p<N; p++)
				tableau[p]=0;

			//lowest_val
			int pos_min=0;

			//figure what the minimal value should be

			//check above
			//get row from pos
			pos_to_row_col(pos,frame,row,col);

			//r=row(pos);
			


			

			pos_min=tableau[pos-1]+1;

			

			//reset
			tableau[pos]=pos_min;

			pos--;
			if(pos<0)
				break;
			else
				continue;
		}
		else
		{
			pos=N-1;
		}
		

		//print clock
		for(auto c : tableau)
		   cout << c <<" ";
	        cout << endl;




		ctr++;
		if(ctr>50)
			break;
	}




	







































	vector<int> Clock;
	//initialize clock
	for(int j=0; j<N; j++)
		Clock.push_back(0);

	//print clock
	for(auto c : Clock)
	   cout << c <<" ";
	cout << endl;


	pos=N-1;
	ctr=0;
	reset=false;
	while(1)
	{
		//advance clock digit
		Clock[pos]++;

		//check if it went too far
		if(Clock[pos]==M)
		{

			//reset lower digits
			for(int p=pos+1; p<N; p++)
				Clock[p]=0;

			//lowest_val
			int pos_min=0;

			//figure what the minimal value should be
			pos_min=Clock[pos-1]+1;

			

			//reset
			Clock[pos]=pos_min;

			pos--;
			if(pos<0)
				break;
			else
				continue;
		}
		else
		{
			pos=N-1;
		}
		

		//print clock
		for(auto c : Clock)
		   cout << c <<" ";
	        cout << endl;




		ctr++;
		if(ctr>50)
			break;
	}





	return 0;


	double S=.5;


	vector<int> occ;
	vector<int> virt;
	vector<int> Pr;
	vector<int> Pk;
	double c;

	/*
	 1 2
	 3 4
	 5


	 1 1
	 2

	 1 1
	 3

	 1 1
	 4

	*/


	int m=4;
	for(int i=0; i<m; i++)
		for(int j=0; j<2; j++)
		{
			
			cout << " counter: " << ctr 
			     << " w/ (i,j)=(" << i << "," << j <<") ->"
			     << 2*i+j << endl;
			ctr++;

		}

	return 0;

	cout << num_weyl(3,3,2) << "\n";
	cout << num_weyl(8,5,2) << "\n";
	cout << num_weyl(4,3,2) << "\n";

	//enumerate_weyl_tabs;
	std::vector<int> L{2,2,1};
	std::vector<int> tableaux;

	
	//initial tableaux
	tableaux.clear();
	for(int n=0; n<N; n++)
		tableaux.push_back(n);
	
	//loop over the tableaux
	/*
	int r=size.L;

	L
	2  , 2	, 1
	T
	1 1, 2 2, 3

	for(int t=0; t<num_weyl(M,N,2*S+1); t++)
	{
		
		if(L[row]==2)
		{
			if(row>0)
				r1min=tableaux[2*(r-1)+0]; //check the previous row
				
			//Anita Gupta, Production Buyer, Wonder Woman 2017
			for(int r1=r1min; r1<M; r1++)
		}
		
		for(int i=0; i<L.size(); i++)
			//loop over indices of each row
			for(int r1=0; r1<M; r1++)
				for(int r2=0; r2<M; r2++)
					;

	}
	*/

	return 0;

	int col1[] = {1,3,5};
	int nc1    = 3;
	int phase  = +1;

	//print antisymmetrizer
	do
	{

		if(phase>0)
		{
			std::cout << "+ ";
			phase=-1;
		}
		else
		{
			std::cout << "- ";
			phase=+1;
		}

		std::cout << col1[0] << ' ' 
		          << col1[1] << ' ' 
		          << col1[2] << ' ' 
		          << col1[3] << std::endl;

	}while(std::next_permutation(col1,col1+nc1)) ;

	//unit tests
	{
	vector<int> left  {1,5,3,7,2,2,6};
	vector<int> right {1,9,5,4,7,3,6};
	auto ndiff=PrimDiffs(10,left,right,Pr,Pk,c);
	std::cout<< "diffs: "<<ndiff;
        }
	/*
EXPECTED DEBUGGING OUTPUT

empty:0 8 
filled:1 3 5 6 7 
virt:4 9 
occ:2 2 
sorting permutation for ref: 5 4 0 2 1 6 3 
input ref: 1 5 3 7 2 2 6 
sorted ref: 2 2 1 3 5 6 7 

sorting permutation for K: 3 1 0 5 2 6 4 
input K: 1 9 5 4 7 3 6 
sorted K: 4 9 1 3 5 6 7 
	 */
	
        {
	vector<int> left  {1,9,3,7,2,2,6};
	vector<int> right {1,2,3,6,9,5,4};
	auto ndiff=PrimDiffs(10,left,right,Pr,Pk,c);
	std::cout<< "diffs: "<<ndiff;
	}
       /*
EXPECTED DEBUGGING OUTPUT

empty:0 8 
filled:1 3 6 9 
virt:4 5 
occ:2 7 
sorting permutation for ref: 5 4 3 0 2 6 1 
input ref: 1 9 3 7 2 2 6 
sorted ref: 2 2 7 1 3 6 9 

sorting permutation for K: 6 5 0 1 2 3 4 
input K: 1 2 3 6 9 5 4 
sorted K: 4 5 1 2 3 6 9 
*/

	{
	vector<int> left  {1,9,3,7,2,2,6};
	vector<int> right {0,2,3,6,9,5,4};
	PrimDiffs(10,left,right,Pr,Pk,c);
	}
	
	vector<int> Peye {0,1,2};
	vector<int> P12  {1,0,2};
	vector<int> P23  {0,2,1};

	vector<int> ans;

	/*
	perm_multiply(P12,P23,ans);

	for(auto i: ans)
		std::cout << " " << i; 
	std::cout << endl;
	*/



	return 0;
}

int 
perm_multiply(const std::vector<int>& P1, const std::vector<int>& P2,          
              std::vector<int>&      result)
{
	if(P1.size()!=P2.size())
	{
		std::cout << "Error in perm_multiply, perms must be same size" << std::endl;
		return 1;
	}

	result.clear();
	result.resize(P1.size());

	for(int j=0; j<P1.size(); j++)
		result[j]=P1[P2[j]];

	return 0;
}

/*
PRIMDIFFS takes 
   1. the total number of orbitals (or the maximally occupied index)
   2. the primative reference 
   3. the primative K
   4. at return Pr has the permuation to put things in order
          
 */ 
int // return excitation level
PrimDiffs(int M, const std::vector<int>& primRef, const std::vector<int>& primK,
	  std::vector<int>&         Pr  , std::vector<int>&         Pk, 
	  double scaling_factor) 
{

    //not using yet
    scaling_factor=1;

    //empty containers
    Pk.clear();
    Pr.clear();

    int     n_ex_lvl=0;
    int     left_occupancy=0;
    int N = primRef.size();
    bool debug=true;
    bool found;

    std::vector<int> unsortedRef;
    unsortedRef.resize(M,0); // resize and fill with val = 0

    const char* basis_fname="basis_funcs";

    libint2::BasisSet shells = parse_basisfile(basis_fname);

    for(auto s: shells)
	    std::cout << s << "\n";

    return 0;
}
