#include<stdio.h>
#include<iostream>
#include<algorithm>
#include"weyl.h"

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
	print_tableau(T,f);

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
		//print_tableau(T,f);


	

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

int 
main(int argc, char *argv[])
{
    bool debug=true;
    if(debug)
    {
	    std::cout.precision(10);
	    std::cout << std::scientific;
    }

    int M,N;
    M=4; N=2; 
    std::cout << "M:" << M <<", N:" << N <<"\n";

    int multiplicity;

    
    if( N % 2 ) multiplicity=2;//doublet
    else        multiplicity=1;//triplet

    int nweyl=num_weyl(M/2,N,multiplicity);
    //int ndets=nchoosek(M,N);
    //int D=nchoosek(M,N);

    if(debug)
    {
	    std::cout << "Number of Weyl tableau: " << nweyl << std::endl;
	    //std::cout << "Number of determinants: " << D << std::endl;
    }

<<<<<<< HEAD

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
    print_tableau(T,frame_rows);
    /*
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
    */

    std::cout << "\n";
    gelfand(M, T, frame_rows);
    

    return 0;
}


    
