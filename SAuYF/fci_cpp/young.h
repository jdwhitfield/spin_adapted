#ifndef YOUNG_ALG
#define YOUNG_ALG

#include<Eigen>
#include<vector>
#include"weyl.h"

class basis_func
{
	public:
		double coeff=1;
		std::vector<int> orbs;
}; // this is the same as the perm_a class right now.

class perm_a
{
	public:
		std::vector<int> perm;
		double coeff=1;
		int    n=0;
	// TODO: overload the multiply operator
	
};


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        Matrix;  // import dense, dynamically sized Matrix type from Eigen;

//BASIC ALGEBRA
perm_a                   perm_multiply(const perm_a& P1, const perm_a& P2);
perm_a                   multiply(const perm_a& P1, const perm_a& P2);
basis_func               multiply(const perm_a P, const basis_func v);
std::vector<perm_a>      multiply(std::vector<perm_a> A, std::vector<perm_a> B);
std::vector<basis_func>  multiply(std::vector<perm_a> O, std::vector<basis_func> wf);
std::vector<basis_func>  multiply(std::vector<perm_a> O, basis_func f0);
std::vector<basis_func>  multiply(const double a, const std::vector<basis_func> F);
std::vector<basis_func>  add(const std::vector<basis_func>,
		             const std::vector<basis_func> );
std::vector<double>      multiply(perm_a P, std::vector<double> v);
std::vector<int>         perm_multiply(const  std::vector<int>& p1, const std::vector<int>& p2);
double                   dot(std::vector<basis_func> L,std::vector<basis_func> R);

//Utilities
std::vector<basis_func>  normalize(const std::vector<basis_func> wf);
perm_a                   perm_I(int n);
bool                     compare_bfs(basis_func f1, basis_func f2);
bool                     compare_perms(perm_a f1, perm_a f2);
int 			 sgn(std::vector<int> perm);
Matrix                   perm_matrix(std::vector<int> perm);
Matrix                   invert_matrix(Matrix X);
std::vector<int>         Frows_to_Fcols(std::vector<int> Frows);
std::vector<int>         invperm(const std::vector<int> T);
std::vector<basis_func>  initial_state(const int N, const std::vector<int> occ_orbs);
	
//IRREP CONSTRUCTION
std::vector<std::vector<basis_func>> 
                         get_irrep_basis(std::vector<int>, std::vector<basis_func>);
Matrix                   get_ortho_matrix_X(std::vector<std::vector<basis_func>> C);
std::vector<perm_a>      wigner_op(int i, int j, int Nelec, int gS, 
		                   std::vector<std::vector<basis_func>> C);

//SYMMETRIC GROUP
std::vector<perm_a>   Ey(std::vector<int> frame_rows,std::vector<int> yT);
std::vector<perm_a>   Acol(std::vector<int> frame,std::vector<int> tableau);
std::vector<perm_a>   Srow(std::vector<int> frame,std::vector<int> tableau);

int                   num_young(int,int);
int                   reset_prev_young(const int pos0, const int M, 
		                       const std::vector<int> frame, 
				       std::vector<int>& tableau);
int                   get_next_ytableau(const std::vector<int> frame, 
		                        std::vector<int>& tableau);
int                   get_next_young_tableau(const std::vector<int> frame,
	           	                std::vector<int>& tableau);
int                   get_init_ytableau(const int M, const std::vector<int> frame,
	        		        std::vector<int>&tableau);
bool                  isvalid_young(std::vector<int> frame_rows,
	                	    std::vector<int> tableau);

//PRINT FUNCTIONS
void  print_perm(perm_a);
void  print_bf(basis_func);
void  print_wf(std::vector<basis_func> wf, bool one_line=true);



#endif
