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
get_next_tableau(const int M, const std::vector<int> frame, std::vector<int>& tableau)
{
	int ctr=0;
	for(int pos=tableau.size()-1; pos>0-1; pos--)
	{
		if(tableau[pos]==pos_max(pos,M,frame,tableau))
		{//at the max allowed, go to next pos
			continue;
		}
		else
		{//not at the max allowed
			//advance tableau position
			tableau[pos]++;

			//try to reset lower positions
			int first_pos_to_be_reset=pos+1;
			auto out=reset_prev(first_pos_to_be_reset,M,frame,tableau);
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
get_init_tableau(const int M, const std::vector<int> frame, std::vector<int>&tableau)
{
	//empty tableau
	tableau.clear();
	

	int ctr=0;
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

	if(ctr>=M)
	{
		cout << "Error. M too small.\n";
		return 1;
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
	//initial tableau


	ntabs=0;

	std::cout << ntabs <<" : ";
	print(tableau,frame);


	while(get_next_tableau(M,frame,tableau)==0)
	{
		std::cout << ++ntabs <<" : ";
	       	print(tableau,frame);
	}
	std::cout << "total num of tableaux: " << ++ntabs;

	return 0;
}






