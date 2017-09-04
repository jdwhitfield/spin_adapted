//Complile command
//g++ -std=c++11 test_eigen.cpp -I/usr/local/include/eigen -llapack -llapacke -o test_eigen
#include<stdio.h>
#include<lapacke.h>
#include<fstream>
#include<Eigen/Eigen>
#include<stdlib.h>
#include<iostream>
#include<assert.h>
#include<iomanip>
#include<random>
#include<chrono>
#include<algorithm>


int 
main()
{
    using namespace std;

    int dim;
    double* matrix_col_major_upper;
    Eigen::MatrixXd matrix;	
    bool debug=true;

    ifstream file("matrix.dat");

    //read the dimension of the matrix in the file
    file >> dim;
    
    //declare space to store values
    matrix_col_major_upper=new double[dim*(dim+1)/2];
    matrix.resize(dim,dim);
    

    int ctr=0;
    double trace;
    for(int i=0; i < dim; i++)         //loop conventions for upper tri part of 
	    for(int j=i; j < dim; j++) //a matrix stored in col-major form
	    {
    		file >> matrix_col_major_upper[ctr];
		matrix(i,j)=matrix_col_major_upper[ctr];
		matrix(j,i)=matrix_col_major_upper[ctr];
		if(i==j)
			trace+=matrix_col_major_upper[ctr];
    		ctr++;
	    }

    
    file.close();
    //done reading file

    //confirm file was read correctly
    if(debug)
    {
	cout << "packed matrix : ";
	for(int j=0; j < dim*(dim+1)/2; j++)
		cout << matrix_col_major_upper[j]<< " ";

	cout << endl;
	
	cout << "unpacked matrix :\n";
	cout << matrix;
        cout << std::endl;

	cout << "trace: "<< trace;
	cout << endl;
    }


     //<< "& idx:" << (i-1 + (j)*(j-1)/2) 
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //LAPACK
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //diagonalization vars
    double* eigenvalues= new double[dim];  
    double* eigenvectors=new double[dim*dim];

    //timing variables
    std::chrono::duration<double> time_elapsed;

    auto start_time = std::chrono::high_resolution_clock::now();

    int out=LAPACKE_dspev(
            LAPACK_ROW_MAJOR, 
            'V',            // 'V' for eigvals and vects or 'N' for just eigvals
            'U',            // input is 'U' upper tri or 'L' lower tri
            dim,            // matrix dimensions
            matrix_col_major_upper, // pointer to matrix array
	    /*
	     * (input/output) DOUBLE PRECISION array, dimension
	     * On entry, the upper or lower triangle of the sym-
	     * metric matrix A, packed columnwise in a linear
	     * array.  The j-th column of A is stored in the array
	     * AP as follows: if UPLO = 'U', AP(i + (j-1)*j/2) =
	     * A(i,j) for 1<=i<=j; if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) 
	     * = A(i,j) for j<=i<=n.
	     */
            eigenvalues,    // on success, its filled with ordered eigvals 
            eigenvectors,   // on success, ordered eigvects (if requested)
            dim             // The leading dimension of the eigvect array >= 1
                            // if JOBZ = 'V', LDZ >= max(1,N).
             );
 
    time_elapsed = std::chrono::high_resolution_clock::now() - start_time;

    cout << "Row-major, U\n";
    for (int j=0; j<dim; j++)
	    cout << eigenvalues[j] << " ";
    cout << "\n";


    
     //convert to row-major upper triangular array from matrix
     for(int i=0; i<dim; i++)
	     for (int j=i;j<dim; j++)
	     {

	     /*Packed columnwise, then The j-th column of A is stored in the array
	     * AP as follows: if UPLO = 'U', AP(i + (j-1)*j/2) =
	     * A(i,j) for 1<=i<=j; if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) 
	     * = A(i,j) for j<=i<=n.
	     */
		    auto idx=(i+j*(j+1)/2) ;
		    matrix_col_major_upper[idx]=matrix(i,j);
	     }


    out=LAPACKE_dspev(
            LAPACK_ROW_MAJOR, 
            'V',            // 'V' for eigvals and vects or 'N' for just eigvals
            'L',            // input is 'U' upper tri or 'L' lower tri
            dim,            // matrix dimensions
            matrix_col_major_upper, // pointer to matrix array
	    /*
	     * (input/output) DOUBLE PRECISION array, dimension
	     * On entry, the upper or lower triangle of the sym-
	     * metric matrix A, packed columnwise in a linear
	     * array.  The j-th column of A is stored in the array
	     * AP as follows: if UPLO = 'U', AP(i + (j-1)*j/2) =
	     * A(i,j) for 1<=i<=j; if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) 
	     * = A(i,j) for j<=i<=n.
	     */
            eigenvalues,    // on success, its filled with ordered eigvals 
            eigenvectors,   // on success, ordered eigvects (if requested)
            dim             // The leading dimension of the eigvect array >= 1
                            // if JOBZ = 'V', LDZ >= max(1,N).
             );
 
    cout << "Row-major, L\n";
    for (int j=0; j<dim; j++)
	    cout << eigenvalues[j] << " ";
    cout << "\n";

    
     //convert to row-major upper triangular array from matrix
     for(int i=0; i<dim; i++)
	     for (int j=i;j<dim; j++)
	     {
		    auto idx=(i+j*(j+1)/2) ;
		    matrix_col_major_upper[idx]=matrix(i,j);
	     }


    out=LAPACKE_dspev(
            LAPACK_COL_MAJOR, 
            'V',            // 'V' for eigvals and vects or 'N' for just eigvals
            'L',            // input is 'U' upper tri or 'L' lower tri
            dim,            // matrix dimensions
            matrix_col_major_upper, // pointer to matrix array
	    /*
	     * (input/output) DOUBLE PRECISION array, dimension
	     * On entry, the upper or lower triangle of the sym-
	     * metric matrix A, packed columnwise in a linear
	     * array.  The j-th column of A is stored in the array
	     * AP as follows: if UPLO = 'U', AP(i + (j-1)*j/2) =
	     * A(i,j) for 1<=i<=j; if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) 
	     * = A(i,j) for j<=i<=n.
	     */
            eigenvalues,    // on success, its filled with ordered eigvals 
            eigenvectors,   // on success, ordered eigvects (if requested)
            dim             // The leading dimension of the eigvect array >= 1
                            // if JOBZ = 'V', LDZ >= max(1,N).
             );
 
    
     //convert to row-major upper triangular array from matrix
     for(int i=0; i<dim; i++)
	     for (int j=i;j<dim; j++)
	     {
		    auto idx=(i+j*(j+1)/2) ;
		    matrix_col_major_upper[idx]=matrix(i,j);
	     }


    cout << "col-major, L\n";
    for (int j=0; j<dim; j++)
	    cout << eigenvalues[j] << " ";
    cout << "\n";

    out=LAPACKE_dspev(
            LAPACK_COL_MAJOR, 
            'V',            // 'V' for eigvals and vects or 'N' for just eigvals
            'U',            // input is 'U' upper tri or 'L' lower tri
            dim,            // matrix dimensions
            matrix_col_major_upper, // pointer to matrix array
	    /*
	     * (input/output) DOUBLE PRECISION array, dimension
	     * On entry, the upper or lower triangle of the sym-
	     * metric matrix A, packed columnwise in a linear
	     * array.  The j-th column of A is stored in the array
	     * AP as follows: if UPLO = 'U', AP(i + (j-1)*j/2) =
	     * A(i,j) for 1<=i<=j; if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) 
	     * = A(i,j) for j<=i<=n.
	     */
            eigenvalues,    // on success, its filled with ordered eigvals 
            eigenvectors,   // on success, ordered eigvects (if requested)
            dim             // The leading dimension of the eigvect array >= 1
                            // if JOBZ = 'V', LDZ >= max(1,N).
             );
 
    
     //convert to row-major upper triangular array from matrix
     for(int i=0; i<dim; i++)
	     for (int j=i;j<dim; j++)
	     {
		    auto idx=(i+j*(j+1)/2) ;
		    matrix_col_major_upper[idx]=matrix(i,j);
	     }


    cout << "col-major, U\n";
    for (int j=0; j<dim*(dim+1)/2; j++)
	    cout << matrix_col_major_upper[j];
    cout << " -> ";
    for (int j=0; j<dim; j++)
	    cout << eigenvalues[j] << " ";
    cout << "\n";





    std::cout << "LAPACK took : " << time_elapsed.count() << " ms " << std::endl;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //Eigen
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    
    start_time = std::chrono::high_resolution_clock::now();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(matrix);
    time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
    if (eigensolver.info() != Eigen::Success) abort();


    std::cout << "Eigen took : " << time_elapsed.count() << std::endl;

    auto evs=eigensolver.eigenvalues();
    
    std::vector<double> evs_eigen(dim);
    std::vector<double> evs_lapack(dim);

    for(int i=0; i<dim; i++)
    {
	    evs_eigen[i]=evs[i];
	    evs_lapack[i]=eigenvalues[i];
    }
    sort(evs_eigen.begin(),evs_eigen.end());
    sort(evs_lapack.begin(),evs_lapack.end());
    

    double diff=0;
    double trace_eigen=0;
    double trace_lapack=0;
    cout << "eigen,lapack\n";
    for(int i=0; i<dim; i++)
    {
	//    diff+=abs(evs_eigen[i]-evs_lapack[i]);
	    cout<< i << " : " << evs_eigen[i]<< "," << evs_lapack[i] << endl;
	    trace_eigen+=evs_eigen[i];
	    trace_lapack+=evs_lapack[i];
    }
		   
    //std::cout << "L_1 Discrepancy: " << diff << endl;
    std::cout << "Trace Eigen :" << trace_eigen << endl;
    std::cout << "Trace LAPACK: " << trace_lapack << endl;

    return 0;
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



