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

	//compute S based on ms=2S+1
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

//convert from row and col of tableau to a linear index
bool
row_col_to_pos(int& pos, const std::vector<int> frame, const int row, const int col)
{
	//input checks, other functions rely on the bool return value for incorrect accesses
	if(row<0 || col <0 || col>1 || row>frame.size()-1)
		return false;

	pos=0;

	for(int r=0; r<row; r++)
		pos=pos+frame[r];

	pos=pos+col;

	return true;
}

int
pos_min(const int pos,const std::vector<int> frame, std::vector<int> tableau)
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
pos_max(const int pos, const int M, const std::vector<int> frame, const std::vector<int> tableau)
{
	int max_val=M-1;
	return max_val;
}
	
int
reset_prev(const int pos0, const int M, const std::vector<int> frame, std::vector<int>& tableau)
{

	bool debug=false;


	if(pos0==tableau.size()) //already at the end, do nothing, no error
		return 0;


	for(int p=pos0; p<tableau.size(); p++)
	{
		//get the min value allowed here based on tableau
		int min_val=pos_min(p,frame,tableau);
		
		//check the value 
		if(min_val>M-1)			//if out of bounds
			return 1;
		else				//if reasonable
			tableau[p]=min_val;
	}
	return 0;
}

int
advance_weyl_tab(std::vector<int>& tableau, const std::vector<int> frame, const int M)
{



	bool debug=true;
	bool reset_flag=false;

	int pos=tableau.size()-1;
	for(int k=0; k < 1000; k++)
	{
		if( tableau[pos] < M-1)
		{
			tableau[pos]++;
			if(debug)
			{
				std::cout << "returning advanced pos " << pos 
				          << " to  T[ans]=" << tableau[pos] << "\n";

				break;
			}

			if(reset_flag)
			{
				;
				//for(int m=frame
			}
		}
		else
		{
			pos--;
			if(pos<0) 
			{
				std::cout <<
				"finished with pos " << pos << " and T[ans]= " << tableau[pos] << "\n";
				break;
			}
			reset_flag=true;
		}
	}

	if(pos<0)
		return 1;
	
	if(pos > -1)
		return 0;
}
int
advance_row(const int row,const int M, std::vector<int>& tableau,const std::vector<int> frame)
{
	bool debug=false;


	if(debug)
		std::cout << "(weyl::advance_row,";
	int pos,col=0;
 	row_col_to_pos(pos,frame,row,col);
	
	//maximum value
	int max_val=M;

	if(debug)
	{
		std::cout << "pos=" << pos << ",";
		std::cout << "row=" << row << ",";
		std::cout << "col=" << col <<",";
		std::cout << "M=" << M <<",";
		std::cout << "F[row]=" << frame[row] << ",";
		std::cout << "T[pos]=" << tableau[pos]
			  << ")";

	}

	if(tableau[pos]==max_val)
	{
		if(debug)
			std::cout << "(at max)";

		return 1;
	}

	if(frame[row]==2)
	{
		tableau[pos]++;
		tableau[pos+1]++;
	}
	else
	{
		tableau[pos]++;
	}


	return 0;
}

int
reset_row(const int row,  const std::vector<int> frame, std::vector<int>& tableau)
{

	int pos, col=0;
 	row_col_to_pos(pos,frame,row,col);
	if(frame[row]==2)
	{
		tableau[pos]=0;
		tableau[pos+1]=0;
	}
	else
	{
		tableau[pos]=0;
	}

	return 0;
}

int 
main()
{
	using namespace std;

	bool debug=true;

	//input variables
	int M=4;
	int N=5;
	std::vector<int> frame{2,2,1};

	std::cout << "Num expected :" <<  num_weyl(M,N,2) << "\n";


	//dummy variables
	int prev_row,row,col,pos,row_idx;
	int ctr,ntabs;

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

	ntabs=0;

	print(tableau,frame);
	ntabs++;

	/*
	std::cout << "tableau.size(): "<<tableau.size();
	std::cout << ",pos0: " << tableau.size()-1 << "\n";
	*/
	//start from right most digit

	/*
	for(int k=0; k<9; k++)
		std::cout << "k :"<<k <<"\n";
	for(int k=9-1; k>0-1; k--)
		std::cout << "k :"<<k <<"\n";
	std::cout << "\n";
	*/

	ctr=0;
	for(pos=tableau.size()-1; pos>0-1; pos--)
	{
		if(tableau[pos]==pos_max(pos,M,frame,tableau))
		{
			continue;
		}
		else
		{
			//advance tableau position
			tableau[pos]++;


			int first_pos_to_be_reset=pos+1;
			auto out=reset_prev(first_pos_to_be_reset,M,frame,tableau);
			

			if(0==out) //everything worked
			{
			       	print(tableau,frame);
			 
				//go to front position
				pos=tableau.size();

				ntabs++;
			}
			
		}
		if(ctr++>100)
			break;
	}
	std::cout << "total num of tableaux: " << ntabs;

	return 0;

	for(int j=0; j<15; j++)
	{
		cout << "out :" << advance_weyl_tab(tableau,frame,M) << "\n";
		print(tableau,frame);
	}


	return 0; 

	row=frame.size()-1;
		// if row < max
		//  advance row
		// else
		//  while(upper row)
		//    try to advance upper row
		//  reset to minimal values
		// end
		//
		// print tableaux



	auto out=advance_row(row,M,tableau,frame);
	
	print(tableau,frame);

	print(tableau,frame);

	for(int k=0; k<10; k++)
	{
	;	
		/* <2017-09-25>, M=3, F={1,1,1}
0 1 2 
output var: 0
0 1 3 
output var: 1
0 1 3 
output var: 1
0 1 3 
		*/
	}


	int r0,r1;
	while(0)
	{
		// if row < max
		//  advance row
		// else
		//  while(upper row)
		//    try to advance upper row
		//  reset to minimal values
		// end
		//
		// print tableaux

		
		// if row < max
		col=0;
		row_col_to_pos(pos,frame,row,col);
		//try to advance
		if(advance_row(row,M,tableau,frame)!=0)
		{
			//if it doesn't work, try to advance the previous row
			prev_row=row;

			ctr=0;
			while(prev_row>0)
			{
				prev_row--;
				if(prev_row==-1)
					break;

			 	if(advance_row(prev_row,M,frame,tableau)==0)
				{ 	//success
					break;
				}

			}

			if(prev_row==-1)
				break;

			//reset lower digits
			for(int r=prev_row+1;r<frame.size(); r++)
			{
				reset_row(r,frame,tableau);	
	
			}

		}	
		else
			if(debug)
			{
				std::cout << "(advanced on the first try)";
			}
		print(tableau,frame);
	}

	return 0;


/*
	{
	int row,col;
	std::vector<int> frame{2,2,1};
        

	pos_to_row_col(0,frame,row,col);
	std::assert
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


	row_col_to_pos(pos,frame,0,0);
	cout << "row " << 0 << " col " << 0 << " : " 
		<< pos << endl;

	row_col_to_pos(pos,frame,0,1);
	cout << "row " << 0 << " col " << 1 << " : " 
		<< pos << endl;

	row_col_to_pos(pos,frame,1,0);
	cout << "row " << 1 << " col " << 0 << " : " 
		<< pos << endl;

	row_col_to_pos(pos,frame,1,1);
	cout << "row " << 1 << " col " << 1 << " : " 
		<< pos << endl;

	row_col_to_pos(pos,frame,0,1);
	row_col_to_pos(pos,frame,1,0);
	row_col_to_pos(pos,frame,1,1);
	row_col_to_pos(pos,frame,2,0);
	row_col_to_pos(pos,frame,2,1);
	row_col_to_pos(pos,frame,3,0);
	}
*/


	pos=N-1;
	ctr=0;
	bool reset=false;
	while(1)
	{
		//advance clock digit
		tableau[pos]++;

		//check if it went too far
		if(tableau[pos] >= M)
		{

			
			//reset lower digits
			for(int p=pos+1; p<N; p++)
			{
				int min_val=pos_min(p,frame,tableau);
				if(min_val==M)
				{
					//not standard
					//restore previous tableau	
					
					min_val=0;
				}


					
			
				tableau[p]=0;
			}



  		        /*******************/	
			//figure what the minimal value should be
  		        /*******************/	
		
  		        /*******************/	


			//reset
			tableau[pos]=pos_min(pos,frame,tableau);

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




	



	return 0;




































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

    //put occupied orbitals into unsorted vector
    for(int j=0; j<N; j++)
	    unsortedRef[primRef[j]]++;

    std::vector<int> unsortedK;
    unsortedK.resize(M,0);

    //put virtual orbitals into unsorted vector
    for(int j=0; j<N; j++)
	    unsortedK[primK[j]]++;


    //these are the orbitals that are occupied in ref but not in K
    std::vector<int> occ;
    //these are the orbitals that are occupied in K but not in ref
    std::vector<int> virt;

    //these are the orbitals that are empty in both
    std::vector<int> empty;
    //these are the orbitals that are occupied in both
    std::vector<int> filled;

    //arrays occ and virt are the orbitals whose occupancies have changed
    
    for(int j=0; j<M; j++)
    {
	    if(unsortedK[j]==unsortedRef[j]) //either filled or empty in both
	    {
		    if(unsortedK[j])
			    filled.push_back(j);
		    else
			    empty.push_back(j);
		    continue;
	    }

	    //at this point the occupancies must differ
	    auto diff=unsortedRef[j]-unsortedK[j];


	    switch(diff)
	    {
		    case -2:
			virt.push_back(j);
		    case -1:
		    	virt.push_back(j);
			break;
		    case 2:
			occ.push_back(j);
		    case 1:
		    	occ.push_back(j);
			break;
		    case 0:
			    std::cerr << "broken logic in PrimDiffs\n"; 
	    }

    }


    if(debug)
    {
	    std::cout << "\n";
	    std::cout << "empty:";
            for(auto i : empty) std::cout << i << " ";
	    std::cout << "\n";

	    std::cout << "filled:";
            for(auto i : filled) std::cout << i << " ";
	    std::cout << "\n";

	    std::cout << "virt:";
            for(auto i : virt) std::cout << i << " ";
	    std::cout << "\n";

	    std::cout << "occ:";
            for(auto i : occ) std::cout << i << " ";
	    std::cout << "\n";
    }

    if(occ.size()>2) // no need to continue as the matrix element is zero
	return 3;


    //sort reference vector, save permutation P_r
    Pr.resize(N,0);

    for (int i = 0 ; i != Pr.size() ; i++) 
	    Pr[i] = i;

    sort(Pr.begin(),Pr.end(),[&](const int& a, const int& b)
  			                   {
						   for(int j=0;j<occ.size();j++)
						   {
						       //sort marked elements to one end
						       if(primRef[a]==occ[j])
							  return true;
						       if(primRef[b]==occ[j])
							  return false;
						   }
					   return (primRef[a] < primRef[b]);
					   });

    //sort K vector, save permutation P_k
    Pk.resize(N,0);
    for (int i = 0 ; i != Pk.size() ; i++) 
	    Pk[i] = i;

    sort(Pk.begin(),Pk.end(),[&](const int& a, const int& b)
  			                   {
						   for(int j=0;j<virt.size();j++)
						   {
						       if(primK[a]==virt[j])
							  return true;
						       if(primK[b]==virt[j])
							  return false;
						   }
						   return (primK[a] < primK[b]);
					   });

    //DEBUGGING/TESTING
    if(debug)
    {
	std::cout << "sorting permutation for ref: ";
	for (auto s : Pr)
		std::cout << s << " ";
	std::cout << std::endl;

	std::cout << "input ref: " ;
	for(int i=0; i!=primRef.size(); i++)
	   std::cout << primRef[i] << " ";
	std::cout << std::endl;

	std::cout << "sorted ref: ";
	for(int i=0; i!=Pr.size(); i++)
	std::cout << primRef[Pr[i]] << " ";
	std::cout << std::endl;

	std::cout << std::endl;

	std::cout << "sorting permutation for K: ";
	for (auto s : Pk)
		std::cout << s << " ";
	std::cout << std::endl;

	std::cout << "input K: " ;
	for(int i=0; i!=primK.size(); i++)
	   std::cout << primK[i] << " ";
	std::cout << std::endl;

	std::cout << "sorted K: ";
	for(int i=0; i!=Pk.size(); i++)
		std::cout << primK[Pk[i]] << " ";
	std::cout << std::endl;

    }

    return occ.size(); // this is the excitation level


}


