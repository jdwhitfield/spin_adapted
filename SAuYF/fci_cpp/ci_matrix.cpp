/******************************************************************************
 * 
 * ci_matrix.cpp
 *
 * This file has functions to form a Hamiltonian and compute its eigenvalues.  
 * 
 * JDWhitfield 
 * Dartmouth, 2016-2017
 * 
 *
 * Parts of the code are borrowed from the open-source projects 
 * Libints and PyQuante under the GNU General Public License.
 *
 ***
 *
 *      This program is free software: you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation, either version 2 of the License, or
 *      (at your option) any later version.
 * 
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 * 
 *      You should have received a copy of the GNU General Public License
 *      along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 *****************************************************************************/


#include<cmath>
#include<stdio.h>
#include<iostream>
#include<string>
#include<vector>
#include<algorithm>
#include<libint2.hpp>
#include"ci_matrix.h"
#include"parser.h"

#include<lapacke.h>

using std::vector;

/// fac[k] = k!, for k < 21
static constexpr std::array<int64_t,21> fac = {{1L, 1L, 2L, 6L, 24L, 120L, 720L, 5040L, 40320L, 362880L, 3628800L, 39916800L, // 12
                                                    479001600L, 6227020800L, 87178291200L, 1307674368000L, 20922789888000L, //5
                                                    355687428096000L, 6402373705728000L, 121645100408832000L, // 3
                                                    2432902008176640000L}}; // 1 
 
size_t max_nprim(const std::vector<Shell>& shells);
int max_l(const std::vector<Shell>& shells);
size_t nbasis(const std::vector<Shell>& shells);

// ******* DEBUGGING ROUTINES
double
print_datum(Matrix M,int i, int j);

double
print_upper_triangular_matrix(int D, double *array);

double
print_datum(Matrix M,int i, int j)
{
    return M(i,j);
}

double
print_upper_triangular_matrix(int D, double *array)
{
    
   using std::cout;

   cout.precision(5);
   cout << std::scientific;

   /*
     if UPLO = 'U', for 1<=i<=j
     AP(i-1 + (j)*(j-1)/2) = A(i,j) 
    */

    int j=1,i=1;

    for(int ctr=0; ctr<D*(D+1)/2; ctr ++)
    {
       cout <<  array[ctr] << "\t";
       //<< "[" << i << "," << j << "] -> ctr:"<< ctr 
       //<< "& idx:" << (i-1 + (j)*(j-1)/2) 
       
       i++;
       if(i>j)
       {
           i=1;
           j++;
           cout << "\n";
       }
       

    }
    return 0;
}
// ******* DEBUGGING ROUTINES

//returns excitation level
int
SlaterDiffs(const vector<int>& detRef, const vector<int>& detK, 
            vector<int>&       virt  , vector<int>&        occ, 
            int& phase)
{
    //empty containers
    virt.clear();
    occ.clear();
    phase=1;

    bool    found;
    int     n_ex_lvl=0;
    int     left_occupancy=0;
    
    //number of electrons
    int N=detRef.size();


    // %%%%%%%% find first orb in reference not in K %%%%%%%%
    for(int occ_orb_idx=0; occ_orb_idx<N; occ_orb_idx++) //loop through occ ref orbs
    {
        //pull occ. orbital label
        auto occ_orb=detRef[occ_orb_idx];

        //flag to say if it was found
        found=false;

        for(int K_orb_idx=0; K_orb_idx<N; K_orb_idx++) //see if its in the K det
        {
            if(occ_orb==detK[K_orb_idx]) // found orb from ref in det K
            {
                found=true;
                break;
            }
        }

        if(!found) // found orbital in ref not in K
        {
            //increase excitation level
            n_ex_lvl++;

            //L_occ1=left_occupancy=orb_index
            left_occupancy=occ_orb_idx;
            
            //% keep looking for second orb in ref not in K
            //%   L_occ2=left_occupancy=orb_index - 1 // extra sub for occ1
            if(n_ex_lvl==2)
                left_occupancy=occ_orb_idx-1;



            //% keep looking for third orb in ref not in K
            //%   if you find it, stop and return as H_JK=0
            if(n_ex_lvl==3)
                return 3;
            
            //save orb 
            occ.push_back(occ_orb);

            if(left_occupancy%2) //if left_occupancy is odd...
                //flip phase
                phase*=-1;
        }
    }	
    
    //no need to go through K if all are orbs the same
    if(n_ex_lvl==0)
        return 0;


    //int nvirts=0;

    // %%%%%%%% find orb in K not in reference %%%%%%%%
    for(int K_orb_idx=0; K_orb_idx<N; K_orb_idx++)
    {
        //pull orbital number
        auto K_orb=detK[K_orb_idx];

        found=false;
        for(int occ_orb_idx=0; occ_orb_idx<N; occ_orb_idx++)
        {
            if(K_orb==detRef[occ_orb_idx])
            {
                found=true;
                break;
            }
        }

        if(!found)
        {

            //save data
            virt.push_back(K_orb);
            
            //taking into account previously found
            left_occupancy=K_orb_idx - (virt.size()-1);

            //accumulate phase
            if(left_occupancy%2)
                phase*=-1;

            //stop looking if we only need one
            if(n_ex_lvl==1)
                return n_ex_lvl;
            
            //stop looking if we found two
            if(virt.size()==2)
                return n_ex_lvl;

        }

    }	
    
 
    printf("\n\nSomething went wrong in SlaterDiffs, logic broken\n");
    throw;
    return -9;

}	

double*
ci_matrix_packed_upper(int M, int N, std::vector<double> MOInts, Matrix h2)
{
    int    D=nchoosek(M,N);
    int    diffs;
    //int    xocc,xvirt,xorbl,xorbr; // spatial orbitals
    //int    xorba,xorbA,xorbb,xorbB;// spatial orbitals
    //int    sA,sa,sB,sb;            // spin variables
    int    reordering_phase;       // reordering phase
    double elem=0;
    vector<int> virt; //orbs in Ldet not in Rdet
    vector<int> occ;  //orbs in Rdet not in Ldet

    //request and clear memory
    double* H=new double[D*(D+1)/2];
    memset(H, 0, D*(D+1)/2 * sizeof(double));
    
    
    for(int Ridx=1; Ridx< D+1; Ridx++)
    {
        //left determinant
        auto detR=CboI(Ridx,M,N);

        /*
        printf("\ndet %i: ",Lidx);
        for( auto l : detL)
            printf("%i ",l);
        printf("\n ");
        */

        for(int Lidx=1; Lidx<Ridx+1; Lidx++)
        {

            //left determinant
            auto detL=CboI(Lidx,M,N);

            //number of differences, needed to evaluate Slater-Condon rules
            diffs=SlaterDiffs(detL,detR,virt,occ,reordering_phase);

            elem=0;

            switch(diffs)
            {
                case 0:
                    elem=  zero_excitation(detL, detR, virt, occ, MOInts, h2);
                    break;
                case 1:
                    elem=single_excitation(detL, detR, virt, occ, MOInts, h2);
                    break;
                case 2:
                    elem=double_excitation(detL, detR, virt, occ, MOInts, h2);
                    break;

                default:
                    //any other number of differences, we can just go to the next
                    continue;
            }

            if(std::abs(elem)<1e-13)
                elem=0;
 
            //include reordering phase
            elem*=reordering_phase;

            //std::cout << "(" << Lidx << "," << Ridx << ")= ["<< Lidx-1+(Ridx)*(Ridx-1)/2 << "]= " << elem << "\n";
            H[Lidx-1+(Ridx)*(Ridx-1)/2]=elem;
            /*if UPLO = 'U', for 1<=i<=j
             * AP(i-1 + (j)*(j-1)/2) = A(i,j) 
             * WARNING: Data is changed during the computation.
             */
        }
    }
    return H;
}

double*
ci_eigenvals(int M, int N, std::vector<double> MOInts, Matrix h2, int DEBUG)
{

    //LAPACK routine vars
    int dim=nchoosek(M,N);
    double* eigenvalues= new double[dim];  
    double* eigenvectors=new double[dim*dim];

    //timing variables
    std::chrono::duration<double> time_elapsed;

    auto start_time = std::chrono::high_resolution_clock::now();
    double* matrix_array = ci_matrix_packed_upper(M,N,MOInts,h2);
    time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
    if(DEBUG)
    {
        std::cout << "\n\n----------------\n";
        std::cout << "Making packed CI matrix took " << time_elapsed.count() << "ms.\n";
        if(DEBUG>1)
        {
            std::cout << "Printing CI matrix:\n";
            print_upper_triangular_matrix(dim,matrix_array);
            std::cout << "\n\n";
        }
    }

 
    start_time = std::chrono::high_resolution_clock::now();
    int out=LAPACKE_dspev(
            LAPACK_COL_MAJOR, 
            'N',            // for eigvals only or 'V' for eigvals and vects
            'U',            // input is upper triangular or 'L' for lower tri
            dim,            // matrix dimensions
            matrix_array,   // pointer to matrix array
            eigenvalues,    // on success, its filled with ordered eigvals 
            eigenvectors,   // on success, ordered eigvects (if requested)
            1               // The leading dimension of the eigvect array >= 1
                            // if JOBZ = 'V', LDZ >= max(1,N).
             );
    time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
    if(DEBUG)
    {

        std::cout << "LAPACK took " << time_elapsed.count() 
                  << "ms to get eigenvalues (return value "  << out << " ).\n";

        for(int kk=0; kk<dim; kk++)
        {
            std::cout << eigenvalues[kk] << "\n";
        }

        if(DEBUG>2)
        {
            std::cout << "\nReference against Eigen code\n";

            start_time = std::chrono::high_resolution_clock::now();
            Matrix EigenH= ci_matrix(M,N,MOInts,h2);
            time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
            std::cout << "Making Eigen CI matrix took " << time_elapsed.count() << "ms.\n";

            start_time = std::chrono::high_resolution_clock::now();
            Eigen::SelfAdjointEigenSolver<Matrix> eigsys(EigenH);
            time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
            std::cout << "getting Eigen eigenvalues took " << time_elapsed.count() << "ms.\n";
            std::cout << eigsys.eigenvalues();
            std::cout << "\n";
        }

    }
    if(out)
        std::cout << "Something went wrong. INFO code: "<< out << "\n";

    delete[] eigenvectors;
    delete[] matrix_array;

    return eigenvalues;
}

/*
int
num_weyl(int M, int N, int ms)
{
	// Hall-Robinson formula
	// See slide 21 of generalized_operators.pdf
	
	double S=(ms -1)/2.0;

	
        auto a=(ms)*1.0/(M+1.0);

	auto k1= .5*N+S+1;
        auto b=nchoosek(M+1,k1);
        
	auto k2= .5*N-S;
        auto c=nchoosek(M+1,k2);

        return  (a*b*c);	
}
*/



Matrix
uga_matrix(int M, int N, int multiplicity, std::vector<double> MOInts, Matrix h2)
{
    // don't support M > 21  
    if(M>21) std::cout << "warning M> 16 is not advisable\n";
    
	
    // 2S+1= g, g-1=2S, S=(g-1)/2
    double S=(multiplicity-1)/2;

    // S= (g-1) / 2

    //number of spin functions for M, N, S= 
    int    D=nchoosek(M,N);
    int    diffs;
    int    reordering_phase;       // reordering phase
    double elem=0;
    vector<int> virt; //orbs in Ldet not in Rdet
    vector<int> occ;  //orbs in Rdet not in Ldet

    Matrix H(D,D);
    H.fill(0);

    
    for(int Lidx=1; Lidx<D+1; Lidx++)
    {
        //left determinant
        auto detL=CboI(Lidx,M,N);

        /*
        printf("\ndet %i: ",Lidx);
        for( auto l : detL)
            printf("%i ",l);
        printf("\n ");
        */

        for(int Ridx=Lidx; Ridx< D+1; Ridx++)
        {
            //left determinant
            auto detR=CboI(Ridx,M,N);

            //number of differences, needed to evaluate Slater-Condon rules
            diffs=SlaterDiffs(detL,detR,virt,occ,reordering_phase);

            elem=0;

            switch(diffs)
            {
                case 0:
                    elem=  zero_excitation(detL, detR, virt, occ, MOInts, h2);
                    break;
                case 1:
                    elem=single_excitation(detL, detR, virt, occ, MOInts, h2);
                    break;
                case 2:
                    elem=double_excitation(detL, detR, virt, occ, MOInts, h2);
                    break;

                default:
                    //any other number of differences, we can just go to the next
                    continue;
            }

            if(std::abs(elem)<1e-13)
                elem=0;
            
            //include reordering phase
            elem*=reordering_phase;

            H(Lidx-1,Ridx-1)=elem;
            if(Lidx!=Ridx)
                H(Ridx-1,Lidx-1)=elem;
        }
    }
    return H;
}

Matrix
ci_matrix(int M, int N, std::vector<double> MOInts, Matrix h2)
{
    int    D=nchoosek(M,N);
    int    diffs;
    int    reordering_phase;       // reordering phase
    double elem=0;
    vector<int> virt; //orbs in Ldet not in Rdet
    vector<int> occ;  //orbs in Rdet not in Ldet

    Matrix H(D,D);
    H.fill(0);
    
    for(int Lidx=1; Lidx<D+1; Lidx++)
    {
        //left determinant
        auto detL=CboI(Lidx,M,N);

        /*
        printf("\ndet %i: ",Lidx);
        for( auto l : detL)
            printf("%i ",l);
        printf("\n ");
        */

        for(int Ridx=Lidx; Ridx< D+1; Ridx++)
        {
            //left determinant
            auto detR=CboI(Ridx,M,N);

            //number of differences, needed to evaluate Slater-Condon rules
            diffs=SlaterDiffs(detL,detR,virt,occ,reordering_phase);

            elem=0;

            switch(diffs)
            {
                case 0:
                    elem=  zero_excitation(detL, detR, virt, occ, MOInts, h2);
                    break;
                case 1:
                    elem=single_excitation(detL, detR, virt, occ, MOInts, h2);
                    break;
                case 2:
                    elem=double_excitation(detL, detR, virt, occ, MOInts, h2);
                    break;

                default:
                    //any other number of differences, we can just go to the next
                    continue;
            }

            if(std::abs(elem)<1e-13)
                elem=0;
            
            //include reordering phase
            elem*=reordering_phase;

            H(Lidx-1,Ridx-1)=elem;
            if(Lidx!=Ridx)
                H(Ridx-1,Lidx-1)=elem;
        }
    }
    return H;
}


// Chemist's notation: [ij|kl]
int 
term4(int i,int j,int k,int l)
{
    //Based on PyQuante implementation
    //similar to algorithm 17.4.2.2 of Cook "Handbook of Computational Quantum Chemistry"
    if(i<j)
    { //swap i and j
        auto temp=i;
        i=j;
        j=temp;
    }

    if(k<l)
    { //swap k and l
        auto temp=k;
        k=l;
        l=temp;
    }
 
    int ij = i*(i+1)/2+j;
    int kl = k*(k+1)/2+l;
    if(ij < kl)
    {
       auto temp= ij;
       ij=kl;
       kl=temp;
    }

    return ij*(ij+1)/2+kl;
}

int
make_matrix(Matrix h2, std::vector<double> MOInts)
{
    int    M=6;
    int    N=3;
    int    D=nchoosek(M,N);
    int    diffs;
    int    xocc,xvirt,xorbl,xorbr; // spatial orbitals
    int    xorba,xorbA,xorbb,xorbB; // spatial orbitals
    int    sA,sa,sB,sb;            // spin variables
    int    phase;
    double elem=0;
    vector<int> virt; //orbs in Ldet not in Rdet
    vector<int> occ;//orbs in Rdet not in Ldet
    
    for(int Lidx=1; Lidx<D+1; Lidx++)
    {
        //left determinant
        auto detL=CboI(Lidx,M,N);


        for(int Ridx=Lidx; Ridx< D+1; Ridx++)
        {
            //left determinant
            auto detR=CboI(Ridx,M,N);


            //get the differences between determinants
            SlaterDiffs(detL,detR,virt,occ,phase);


            /* TEST FOR SlaterDiffs()
            std::cout << "< ";
            for(auto l : detL)
                std::cout<< l << " ";
            std::cout << "| ";
            for(auto r : detR)
                std::cout<< r << " ";
            std::cout << "> \t ( ";
            for(auto orbr : occ)
                std::cout << orbr << " ";
            std::cout << ") ";

            std::cout << "->\t ( ";
            for(auto orbl : virt)
                std::cout << orbl << " ";
            std::cout << ")\n";
            */
            
            //number of differences, needed to evaluate Slater-Condon rules
            diffs=virt.size();

            elem=0;

            switch(diffs)
            {
                case 0:
                    //  H = \sum_m^{occ} h_{mm}+.5*\sum_{nm}^{occ} ([mm|nn]-[mn|nm])

                    for( auto orbl : detL )
                    {
                        //convert to spatial orbital
                        xorbl=ceil(orbl/2);

                        elem+=h2(xorbl,xorbl);

                       
                        for(auto orbr :detR)
                        {
                            xorbr=ceil(orbr/2);

                            //coulomb term
                            //elem+=.5*[xorbl,xorbl|xorbr,xorbr];
                            elem+=.5*MOInts[term4(xorbl,xorbl,xorbr,xorbr)];
                            
                            //exchange, if the spins agree
                            //elem-=.5*[xorbl,xorbr|xorbl,xorbr];
                            if(orbl%2 ==orbr%2) 
                                elem-=.5*MOInts[term4(xorbl,xorbr,xorbl,xorbr)];

                        }
                    }

                    break;
                case 1:
                    //  H = h_{Aa}+\sum_n^{occR} ( [Aa|nn]-[An|na] )
                    
                    if(virt[0]%2 != occ[0]%2)
                    {
                        elem=0;
                        continue;
                    }

                    xocc =ceil(occ[0] /2);
                    xvirt=ceil(virt[0]/2);

                    //one-body part
                    elem+=h2(xvirt,xocc);

                    //two-body part
                    for(auto orbr : detR)
                    {
                        //occupied orbitals
                        xorbr=ceil(orbr/2);

                        //Coulomb analog
                        elem+=MOInts[term4(xocc,xvirt,xorbr,xorbr)];
                        

                        //exchange analog, if spins agree
                            //elem-= [xocc,xorbr|xvirt,xorbr];
                        if(orbr%2 == occ[0]%2)
                            elem-=MOInts[term4(xocc,xorbr,xvirt,xorbr)];
                    }
                    break;
                case 2:
                    //  H = [Aa|Bb]-[Ab|Ba]
                    
                    //check spins
                    xorbA=ceil(virt[0]/2);
                    sA   =virt[0]%2;

                    xorbB=ceil(virt[1]/2);
                    sB   =virt[1]%2;

                    xorba=ceil(occ[0]/2);
                    sa   =occ[0]%2;

                    xorbb=ceil(occ[1]/2);
                    sb   =occ[1]%2;

                    /*
                    elem=[xorbA,xorba|xorbB,xorbb] * (sA==sa && sB==sb)
                        -[xorbA,xorbb|xorbB,xorba] * (sA==sb && sB==sa)
                    */
                    elem=MOInts[term4(xorbA,xorba,xorbB,xorbb)] * (sA==sa && sB==sb)
                        -MOInts[term4(xorbA,xorbb,xorbB,xorba)] * (sA==sb && sB==sa);

                    break;
            }

            if(std::abs(elem)>1e-10)
                printf("%d , %d , %lf\n",Lidx,Ridx,elem);
        }
    }
    return 0;
}

//simple nchoosek function for integers
int 
nchoosek(double n, double k)
{

    if(k>n)
    {
        printf("\n\nInvalid k in {n \\choose k} function, n=%lf, k=%lf \n\n",n,k);
        throw; 
    }

	if(k==0)
		return 1;
	if(k<0)
	{

        printf("\n\nWarning: k in {n \\choose k} function,  k=%lf<0. Returning 0.\n\n",k);
		return 0;
	}


    if(n<21) //use lookup table
        return(round(fac[n]/(fac[n-k]*fac[k])));

    if(k>n/2) // use smaller k if available
        k=n-k;

    if(k<=1)
        return(pow(n,k));
    else
    {
        double c=1;
        double denom=1;
        for(int numer=n-k+1; numer<n+1; numer++)
        {
            c*=numer/denom;
            denom++;

        }
        return round(c);
    }
}

double
zero_excitation(const vector<int> detL,  const vector<int> detR,
                const vector<int> virt,  const vector<int> occ,
                std::vector<double> MOInts, Matrix h2,
                bool spatial_orbs)
{
    //  H = \sum_m^{occ} h_{mm}+.5*\sum_n^{occ} ([mm|nn]-[mn|nm])
    double elem=0;

    if(spatial_orbs)
    {
        for( auto orbl : detL )
        {
            //convert to spatial orbital
            auto xorbl=spin2spatial(orbl);

            elem+=h2(xorbl,xorbl);

           
            for( auto orbr :detR)
            {
                auto xorbr=spin2spatial(orbr);

                //coulomb term
                //elem+=.5*[xorbl,xorbl|xorbr,xorbr];
                elem+=.5*MOInts[term4(xorbl,xorbl,xorbr,xorbr)];
                
                //exchange, if the spins agree
                //elem-=.5*[xorbl,xorbr|xorbl,xorbr];
                if(orbl%2 ==orbr%2) 
                    elem-=.5*MOInts[term4(xorbl,xorbr,xorbl,xorbr)];

            }
        }
    }
    else
    {
        std::cout<< "Spin orbitals only has not be implemented. Use spatial orbitals!";
    }
    return elem;
}

double
single_excitation(const vector<int> detL,  const vector<int> detR,
                  const vector<int> virt,  const vector<int> occ,
                  std::vector<double> MOInts, Matrix h2,
                  bool spatial_orbs)
{
    //  H = h_{Aa}+\sum_n^{occR} ( [Aa|nn]-[An|na] )
    double elem=0;

    if(!spatial_orbs)
    {
        std::cout << "Spin orbitals only not implemented. Use spatial orbitals!\n";
        throw;
    }

    if(virt[0]%2 != occ[0]%2)
        return 0;


    auto xocc =spin2spatial(occ[0]);
    auto xvirt=spin2spatial(virt[0]);

    //one-body part
    elem+=h2(xvirt,xocc);

    //two-body part
    for(auto orbr : detR)
    {
        //occupied orbitals
        auto xorbr=spin2spatial(orbr);

        //Coulomb analog
        elem+=MOInts[term4(xocc,xvirt,xorbr,xorbr)];
        
        //exchange analog, if spins agree
            //elem-= [xocc,xorbr|xvirt,xorbr];
        if(orbr%2 == occ[0]%2)
	{
            elem-=MOInts[term4(xocc,xorbr,xvirt,xorbr)];
	}
    }

    return elem;
}

double
double_excitation(const vector<int> detL,  const vector<int> detR,
                  const vector<int> virt,  const vector<int> occ,
                  std::vector<double> MOInts, Matrix h2,
                  bool spatial_orbs)
{

    if(!spatial_orbs)
    {
        std::cout << "Spin orbitals only not implemented. Use spatial orbitals!\n";
        throw;
    }

    //check spins
    auto xorbA=spin2spatial(virt[0]);
    auto xorba=spin2spatial( occ[0]);
    auto xorbB=spin2spatial(virt[1]);
    auto xorbb=spin2spatial( occ[1]);
    auto sA   =virt[0]%2;
    auto sB   =virt[1]%2;
    auto sa   =occ[0]%2;
    auto sb   =occ[1]%2;
    
    /*
    elem=[xorbA,xorba|xorbB,xorbb] * (sA==sa && sB==sb)
        -[xorbA,xorbb|xorbB,xorba] * (sA==sb && sB==sa)
    */
    return MOInts[term4(xorbA,xorba,xorbB,xorbb)] * (sA==sa && sB==sb)
          -MOInts[term4(xorbA,xorbb,xorbB,xorba)] * (sA==sb && sB==sa);


   

}



void
SlaterDiffs_OLD(const vector<int> detL,const vector<int> detR, vector<int>& virt, vector<int>& occ)
{
    //empty containers
    virt.clear();
    occ.clear();
    bool found;
    
    //find the number and location of differences
    for(unsigned int l=0; l<detL.size();l++)
    {
        //pull orbital label
        auto orbl=detL[l];

        //flag to say if it was found
        found=false;

        //loop through the orbs of detR
        for(unsigned int r=0; r<detL.size();r++)
            if(orbl==detR[r])
            {
                //when orbl is in detR
                //mark it and quit searching
                found=true;
                break;
            }
        //if it wasn't found save it
        if(!found)
            virt.push_back(orbl);

        //if were more than 2 differences, we can stop because Helem=0
        if(virt.size()>2)
            return;

    }


    
    //we need to find the occupied orbitals
    for(unsigned int r=0; r<detR.size(); r++)
    {
        found=false;

        //pull orbital label
        auto orbr=detR[r];

        //loop over orbs of detL
        for(unsigned int l=0; l<detL.size(); l++)
            if(orbr==detL[l])   
            {
               found=true;
               break;
            }
        if(!found)
            occ.push_back(orbr);
    }

    //check for canonical ordering
    if(virt.size()==2)
    {
        //of virtual
        if(virt[0]<virt[1])
        {
            double temp=virt[0];
            virt[0]=virt[1];
            virt[1]=temp;
        }
        //and of occupied orbitals
        if(occ[0]>occ[1])
        {
            double temp=occ[0];
            occ[0]=occ[1];
            occ[1]=temp;
        }
    }


    return;
}

/*
%CBOI returns the combo from the index
%the inputs are the index, v is an integer and M is the number of items
% when v is an integer 1 to v is used as the ensemble to select from
%
% Examples CboI(1,10,3) = [1,1,1]
%          CboI(105,10,3) = [5,7,8]
%          CboI(1,[1,2,5,8],3) = [1,2,5]
%
%based on stackoverflow.com/questions/13730772 by Origin

%VCQ
%JDWhitfield 2012-2013, 2015
*/
std::vector<int>
CboI(int idx, int N, int M)
{

    //variables
    int k=idx;
    int l=0;
    int C;
    
    //combination
    std::vector<int> cbo;
    cbo.reserve(M);

    for(int i=1; i<M+1; i++)
    {
        // %determine combination(i)
        while(1)
        {
            l++;
            C=nchoosek(N-l,M-i);

            if(k>C)
            {
                k=k-C;
            }
            else
            {
                cbo.push_back(l);
                break;
            }
        }
    }

    return cbo;
}


int
idx(int n,int m, int M)
{
    return (n)*M/2+m;
}

int 
spin2spatial(int spinorb)
{

    /********************************
     *  orb       CIbo      spatial 
     *   1a        1           0
     *   1b        2           0
     *   2a        3           1
     *   2b        4           1
     *   .         .           .
     *   .         .           .
     *   .         .           .
     *********************************/
    return(std::floor((spinorb-1)/2));
}

// Helper functions copied from libint's implementation in hartree-fock.cc
size_t nbasis(const std::vector<Shell>& shells) 
{
  size_t n = 0;
  for (const auto& shell: shells)
    n += shell.size();
  return n;
}

size_t max_nprim(const std::vector<Shell>& shells) 
{
  size_t n = 0;
  for (auto shell: shells)
    n = std::max(shell.nprim(), n);
  return n;
}

int max_l(const std::vector<Shell>& shells) 
{
  int l = 0;
  for (auto shell: shells)
    for (auto c: shell.contr)
      l = std::max(c.l, l);
  return l;
}

/// return the map from basis function index to shell index
std::vector<size_t> map_basis_functions_to_shell(const std::vector<Shell>& shells)
{
    std::vector<size_t> result;
    result.reserve(nbasis(shells));

    size_t n=0;
    for(auto s : shells)
    {
        for(unsigned int i=0; i<s.size(); i++)
            result.push_back(n);
        n++;
    }
    return result;
}

/// @return the map from shell index to index of the first basis function from this shell
/// \note basis functions are ordered as shells, i.e. shell2bf[i] >= shell2bf[j] iff i >= j
std::vector<size_t> map_shell_to_basis_function(const std::vector<Shell>& shells) 
{

  std::vector<size_t> result;

  //reserves space in the vector for a certain number of items
  result.reserve(shells.size());

  // size_t  type is the unsigned integer type that is the result of the sizeof 
  // operator (and the offsetof operator), so it is guaranteed to be big enough 
  // to contain the size of the biggest object your system can handle (e.g.,
  // a static array of 8Gb). 
  //
  //  from http://stackoverflow.com/questions/131803/unsigned-int-vs-size-t 
  size_t n = 0;
  for (auto shell: shells) {
    result.push_back(n);
    n += shell.size();
  }

  return result;
}


std::vector<double>
transform4(const int M, const Matrix& C, const vector<double>& AOInts, int debug)
{
    
    if(C.isIdentity())
        return AOInts; //nothing to do

    Matrix Cdag=Matrix(C.adjoint());

    //expanded Int arrays
    //
    //indexing: 
    //  [pq|rs] = Ints(idx(p,q,M),idx(r,s,M))
    Matrix Ints1(M*M/4,M*M/4);
    Matrix Ints2(M*M/4,M*M/4);
    Ints1.fill(0);
    Ints2.fill(0);

    //number of 2e- integrals
    int m=(M/2) - 1;
    int nints=term4(m,m,m,m)+1;

    std::vector<double> MOInts(nints, 0.0);
    double intval;

   
    //AOInts -> Ints1 = [p jj | kk ll]
    //printf("\nTransform 1st index\n");
    for(int p=0; p<M/2; p++)
        for(int jj=0; jj<M/2; jj++)
            for(int kk=0; kk<M/2; kk++)
                for(int ll=0; ll <M/2; ll++)
                {
                    intval=0;
                    for(int ii=0; ii<M/2; ii++)
                    {
                        intval+=Cdag(p,ii)*AOInts[term4(ii,jj,kk,ll)];
                    }
                    Ints1(idx(p,jj,M),idx(kk,ll,M))=intval;
                }

    // Ints1 -> Ints2 = [p q | kk ll]
    //printf("\nTransform 2nd index\n");
    for(int p=0; p<M/2; p++)
        for(int q=0; q<M/2; q++)
            for(int kk=0; kk<M/2; kk++)
                for(int ll=0; ll <M/2; ll++)
                {
                    intval=0;
                    for(int jj=0; jj<M/2; jj++)
                    {
                        intval+=
                            C(jj,q)*Ints1(idx(p,jj,M),idx(kk,ll,M));
                        }
                    Ints2(idx(p,q,M),idx(kk,ll,M))=intval;
                }
 
    // Ints2 -> Ints1 = [p q | r ll]
    //printf("\nTransform 3rd index\n");
    for(int p=0; p<M/2; p++)
        for(int q=0; q<M/2; q++)
            for(int r=0; r<M/2; r++)
                for(int ll=0; ll <M/2; ll++)
                {
                    intval=0;
                    for(int kk=0; kk<M/2; kk++)
                    {
                        intval+=
                            Cdag(r,kk)*Ints2(idx(p,q,M),idx(kk,ll,M));
                        
                    }
                    Ints1(idx(p,q,M),idx(r,ll,M))=intval;
                }

    // Ints1 -> MO = [p q | r s]
    //printf("\nTransform 4th index\n");
    for(int p=0; p<M/2; p++)
        for(int q=0; q<=p; q++)
            for(int r=0; r<=p; r++)
                for(int s=0; s <= (p==r ? q : r); s++)
                {
                    intval=0;
                    for(int ll=0; ll<M/2; ll++)
                    {
                        intval+=
                            C(ll,s)*Ints1(idx(p,q,M),idx(r,ll,M));
                    }
                    if(debug)
                        printf("%i: [ %i %i | %i %i ] = %lf\n",term4(p,q,r,s),p,q,r,s,intval);
                    MOInts[term4(p,q,r,s)]=intval;
                }


    return MOInts;
}



// OLD MAIN ROUTINE
/*
int
main(int argc, char *argv[])
{
    //timing variables
    std::chrono::system_clock::time_point start_time;
    std::chrono::duration<double> time_elapsed;

    // ****************************
    // * Parse commandline inputs *
    // ****************************
    if (argc > 1)
    {	
	    if(!strcmp(argv[1],"-h"))
	    {
            std::cout << "ci_matrix [basis_func] [nuc_field]";
            std::cout << "\n\tDefault files are used when called with no parameters.";
		    return 0;
	    }
    }

    const auto basis_fname=(argc>1) ? argv[1] : "basis_funcs";
    const auto nuc_fname  =(argc>2) ? argv[2] : "nuc_field";

    // ****************************
    // * Parse data files         *
    // ****************************
    auto shells = parse_basisfile(basis_fname);

    printf("Basis set\n");

    for( auto s : shells)
        std::cout << s << "\n";

    assert(max_l(shells)==0);
    
    //number of spin orbitals in the basis set
    int M=2*nbasis(shells);
    
    //number of electrons
    int N=-99;

    //number of 2e- integrals
    int m=(M/2) - 1;
    int nints=term4(m,m,m,m)+1;

    auto atoms = parse_nucfile(nuc_fname,N);

    if(N==-99)
    {
        std::cout << "Nelec not set\n";
        throw;
    }

    
    Matrix S, h;
    std::vector<double> AOInts;

    //libints
    get_libints(shells,atoms,S,h,AOInts);


    std::cout << "S\n" << S << "\n";
    std::cout << "h\n" << h << "\n";
    std::cout << "AO 2body integrals\n";
    for(auto p=0; p!=M/2; p++) //unique integral labels, looping scheme from libint
        for(auto q=0; q<=p; q++)
            for(auto r=0; r<=p; r++)
                for(auto s=0; s<= (p==r ? q : r) ; s++)
                    if(std::abs(AOInts[term4(p,q,r,s)]) > 1e-5)
                        printf("%i: [%i%i|%i%i] = %lf\n",term4(p,q,r,s),p,q,r,s,AOInts[term4(p,q,r,s)]);



    Matrix h2;
    std::vector<double> MOInts;
    if(Matrix(S).isIdentity(1e-9))
    {
        //nothing changes if S is the identity matrix
        h2=h;
        MOInts=AOInts;
        std::cout<<"We need to orthogonalize\n";
    }
    else
    {// *     ORTHOGONALIZATION    *


        / *********************************************************** /
        // Transformation matricies
        //
        //Basis transform to an orthogonal space
        //C = C(W) = XW for any WW^\dag =\id
        //

        Eigen::SelfAdjointEigenSolver<Matrix> Eigensystem(S);
        Matrix s=Matrix(Eigensystem.eigenvalues());
        //check which eigenvalues are non-trivial to avoid linear dependence
        //then compute s^{-1/2}
        Matrix shalf;
        shalf.resize(M/2,M/2);    // set size
        shalf.fill(0);            // make sure all elements are zero
        //make sure its a column vector
        if(s.rows()<s.cols())
            s.transposeInPlace();

        for(int j=0; j<s.rows(); j++)
        {
            if(std::abs(s(j,0))>1e-6)
                shalf(j,j)=1/sqrt(s(j,0));
            else
                shalf(j,j)=0;
        }

        //orthogonalization transform
        //Here we're just using W=\id
        Matrix C;
        bool CANONICAL=true;
        if(CANONICAL)
            C=Matrix(Eigensystem.eigenvectors()*shalf);
        else
            C=Matrix(Eigensystem.eigenvectors()*shalf*Eigensystem.eigenvectors().adjoint());

        //time the transformations
        start_time = std::chrono::high_resolution_clock::now();

        h2       =Matrix(C.adjoint())*h*C;   //One-body transform
        MOInts   = transform4(M,C,AOInts);//two-body transform

        time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
        std::cout << "Transformation to orthogonal basis took " << time_elapsed.count() << "ms \n";
        / *
        std::cout << "Ux\n" << Eigensystem.eigenvectors();
        std::cout << "\nLx\n" << Eigensystem.eigenvalues();
        printf("\n");

        std::cout << "X ";
        std::cout <<  CANONICAL 
                      ? "via canonical orthogonalization\n" 
                      : "X via symmetric orthogonalization\n";
        std::cout <<  X << "\n";
        * /

    }


 


    std::cout << h2 << "\n";
    for(auto m : MOInts)
        std::cout << m << "\n";

    Matrix H= ci_matrix(M,N,MOInts,h2);

    printf("Hamiltonian:\n");
    std::cout << H << "\n\n";

    printf("Diagonal elements\n");
    std::cout << H.diagonal() << "\n";
    

    Eigen::SelfAdjointEigenSolver<Matrix> eigsys(H);
    Matrix eigvals=Matrix(eigsys.eigenvalues());
    printf("\n\nEigenvalues: \n");
    for(int i=0; i<eigvals.size(); i++)
        printf("%15.10lf\n",eigvals(i));



    return 0;


}
*/
