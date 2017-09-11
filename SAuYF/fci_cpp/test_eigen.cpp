//Complile command
//g++ -std=c++11 test_eigen.cpp -I/usr/local/include/eigen -llapack -llapacke -o test_eigen
#include<stdio.h>
#include<fstream>
#include<stdlib.h>
#include<iostream>
#include<assert.h>
#include<iomanip>
#include<random>
#include<chrono>
#include<algorithm>
#include<lapacke.h>
#define EIGEN_USE_LAPACKE
#include<Eigen/Eigen>

void
reset_array(int dim, Eigen::MatrixXd matrix, double*& array);

int 
main()
{
    using namespace std;

    int dim;
    double* matrix_col_major_upper;
    Eigen::MatrixXd matrix;	
    bool debug=false;

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
		if(debug)
			cout << "read "<< ctr << " as " << matrix_col_major_upper[ctr] << " into (" << i << ","<<j<< ") with ";
		matrix(i,j)=matrix_col_major_upper[ctr];
		matrix(j,i)=matrix_col_major_upper[ctr];

		if(i==j)
			trace+=matrix_col_major_upper[ctr];

		if(debug)
	           cout << "\t" << ctr << " : " << (i*dim+j-i*(i+1)/2) << "\n";
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


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //LAPACK
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //diagonalization vars
    double* eigenvalues= new double[dim];  
    double* eigenvectors=new double[dim*dim];
    int out;

    //timing variables
    std::chrono::duration<double> time_elapsed;

    if(debug)
    {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    reset_array(dim, matrix, matrix_col_major_upper);
    out=LAPACKE_dspev(
            LAPACK_ROW_MAJOR, 
            'V',            // 'V' for eigvals and vects or 'N' for just eigvals
            'U',            // input is 'U' upper tri or 'L' lower tri
            dim,            // matrix dimensions
            matrix_col_major_upper, // pointer to matrix array
            eigenvalues,    // on success, its filled with ordered eigvals 
            eigenvectors,   // on success, ordered eigvects (if requested)
            dim             // The leading dimension of the eigvect array >= 1
                            // if JOBZ = 'V', LDZ >= max(1,N).
             );
 

    cout << "Row-major, U\n";
    for (int j=0; j<dim; j++)
	    cout << eigenvalues[j] << " ";
    cout << "\n";
    
    reset_array(dim, matrix, matrix_col_major_upper);
   
  
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    out=LAPACKE_dspev(LAPACK_ROW_MAJOR, 'V',
		    'L', dim,  matrix_col_major_upper, eigenvalues, eigenvectors, dim);
 
    cout << "Row-major, L\n";
    for (int j=0; j<dim; j++)
	    cout << eigenvalues[j] << " ";
    cout << "\n";

    reset_array(dim, matrix, matrix_col_major_upper);
  
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    out=LAPACKE_dspev(LAPACK_COL_MAJOR, 'V',
		    'L', dim,  matrix_col_major_upper, eigenvalues, eigenvectors, dim);
 
    cout << "Col-major, L\n";
    for (int j=0; j<dim; j++)
	    cout << eigenvalues[j] << " ";
    cout << "\n";

    reset_array(dim, matrix, matrix_col_major_upper);
  
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    out=LAPACKE_dspev(LAPACK_COL_MAJOR, 'V',
		    'U', dim,  matrix_col_major_upper, eigenvalues, eigenvectors, dim);
 
    cout << "Col-major, U\n";
    for (int j=0; j<dim; j++)
	    cout << eigenvalues[j] << " ";
    cout << "\n";

    reset_array(dim, matrix, matrix_col_major_upper);
  
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    cout << "\n";
    
    }

    auto start_time = std::chrono::high_resolution_clock::now();

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
 
    time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
    std::cout << "LAPACK took : " << time_elapsed.count() << " ms " << std::endl;
    
    if(debug)
    {
    //convert to row-major upper triangular array from matrix
    for(int i=0; i<dim; i++)
	     for (int j=i;j<dim; j++)
	     {
		    auto idx=(i+j*(j+1)/2) ;
		    matrix_col_major_upper[idx]=matrix(i,j);
	     }


    cout << "row-major, U\n";
    //for (int j=0; j<dim*(dim+1)/2; j++)
    //	    cout << matrix_col_major_upper[j];
    //cout << " -> ";
    for (int j=0; j<dim; j++)
	    cout << eigenvalues[j] << " ";
    cout << "\n";
    }





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
	    diff+=abs(evs_eigen[i]-evs_lapack[i]);
	    if(debug)
		    cout<< i << " : " << evs_eigen[i]<< "," << evs_lapack[i] << endl;
	    trace_eigen+=evs_eigen[i];
	    trace_lapack+=evs_lapack[i];
    }
		   
    std::cout << "L_1 Discrepancy: " << diff << endl;

    return 0;
}

void
reset_array(int dim, Eigen::MatrixXd matrix, double* &array)
{
//convert to row-major upper triangular array from matrix
    for(int i=0; i<dim; i++)
	     for (int j=i;j<dim; j++)
	     {
		    auto idx=(i*dim+j-i*(i+1)/2);
		    array[idx]=matrix(i,j);
	     }
return;
}
    

