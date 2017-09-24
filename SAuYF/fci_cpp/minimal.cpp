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

int 
main()
{
	using namespace std;


	vector<int> occ;
	vector<int> virt;
	vector<int> Pr;
	vector<int> Pk;
	double c;

	/*
	 1 2
	 3 4
	 5
	*/

	cout << num_weyl(3,3,2) << "\n";
	cout << num_weyl(8,5,2) << "\n";
	cout << num_weyl(4,3,2) << "\n";

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



