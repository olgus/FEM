#pragma once
#include "Solver_definitions.h"
#include "Matrix.h"
#include <mkl_spblas.h>
#include <string.h>

#include <mkl.h>

const double SM_PRECIS = 1e-15;
const double MORE_PRECIS = 1e-30;
const double LIL = 1e-50;
const int ITERMAX = 5000;
const double LESS_PRECIS = 1e-7;

class compressed_matrix // compressed_matrix, non-symmetrical
{
private:
	int size; // square matrix, btw
	int n_entries; // amount of potentially non-zero entries

	// second set goes for factorization, so i don't have to change original matrix's data
	MathVector * ggl, *ggu, *di, *di_inv; // l, u, di entries
	MathVector * L, *U, *Ld; // L, U, Ld entries
	MathVector *x; // x vectors
	int * ig, *jg; // ig and jg arrays 

	// CSR format
	bool use_CSR;
	int * ia, * ja;
	double * a;
	void refresh_CSR ();

	void mult_U_v (const MathVector & v, MathVector * res); // multiply U-section of the matrix by v and put into res
	void inverse_di ();

	void CSRC_add_to_entry (int I, int J, double value); // adds value to entry I, J (as in dense matrix)
	void CSRC_set_entry (int I, int J, double value); // sets value to entry I, J (as in dense matrix)
	void CSRC_multiply_entry (int I, int J, double value); // multiplies entry I, J (as in dense matrix) to value
	double CSRC_get_entry (int I, int J); // returns value of the entry I, J (as in dense matrix)

	void CSR_add_to_entry (int I, int J, double value); // adds value to entry I, J (as in dense matrix)
	void CSR_set_entry (int I, int J, double value); // sets value to entry I, J (as in dense matrix)
	void CSR_multiply_entry (int I, int J, double value); // multiplies entry I, J (as in dense matrix) to value
	double CSR_get_entry (int I, int J); // returns value of the entry I, J (as in dense matrix)
public:
	MathVector *x0, *f; // f, x0 vectors
	void Inner_copy (); // copy di, ggl and ggu into Ld, L, U respectively

	compressed_matrix (); // constructor 
	compressed_matrix (int N, int N_L); // constructor that sets amount of elements
	compressed_matrix (const compressed_matrix & crm); // constructor-copy
	~compressed_matrix (); // destructor

	compressed_matrix & operator= (compressed_matrix & crm); // assignment operator

	void set_matrix_size (int Size, int N_entries); // set size and n_entries
	int Size (); // returns dense size of the matrix
	void get_from_files (); // set matrix from files
	void get_from_files (char * folder); // set matrix from files
	void set_x0 (double value); // sets x0 to vector of values
	bool set_ig_jg (int * Ig, int * Jg); // set ig and kg from arrays
	void Clear (); // fills ggl, ggu, di and f with zeroes

	// sequential realisations
	void mult_A_v (const MathVector & v, MathVector * res); // multiply matrix by v and put into res		
	void calc_LvR (const MathVector & r, MathVector * v); // calc L * v = r. v is the result of the function
	void calc_URV (const MathVector & v, MathVector * r); // calc U * r = v. r is the result of the function
	void calc_LTRV (const MathVector & v, MathVector * r); // calc LT * r = v. r is the result of the function
	void calc_DvR (const MathVector & r, MathVector * v); // calc v = D-1 R 

	// parallel realisations, personal
	void mult_A_v (const MathVector & v, MathVector * res, double * res_omp, int par_num_threads); // multiply matrix by v and put into res		
	void calc_DvR (const MathVector & r, MathVector * v, int par_num_threads); // calc v = D-1 R 

	// MKL realisations
	void mult_A_v_MKL (const sparse_matrix_t SMA, const MathVector & v, MathVector * res); // multiply matrix by v and put into res		
	void calc_DvR_MKL (const MathVector & r, MathVector * v); // calc v = D-1 R 
	void calc_LvR_MKL (const sparse_matrix_t SMA, const MathVector & r, MathVector * v); // calc L * v = r. v is the result of the function
	void calc_URV_MKL (const sparse_matrix_t SMA, const MathVector & v, MathVector * r); // calc U * r = v. r is the result of the function

	// sequential methods
	int solve_GMRES (int d_type, int max_depth); // approved for both
	int solve_LOS (int d_type); // approved for both
	int solve_CGM_symm (int d_type); // approved for both
	void solve_LU (); // approved 
	void solve_LL (); // approved 

	// parallel methods, personal functions
	int solve_LOS (int d_type, int par_num_threads); // approved for both
	int solve_GMRES (int d_type, int max_depth, int par_num_threads); // approved for both
	int solve_CGM_symm (int d_type, int par_num_threads); // approved for both

	// use CSR
	void convert_to_CSR (); // sets ia, ja, a and csr flag
	void make_lut (int *Lia, int *Lja, double * Lt, int *Uia, int *Uja, double * Ut); // saves factorized data as upper triangle with 1-diagonal and lower tirangle with di-diagonal
	int solve_LOS_MKL (int d_type, int par_num_threads); // LOS for CSR with MKL functions
	void solve_pardiso ();
	void solve_pardiso (int par_num_threads);

	void print (); // for test purposes
	void fprint (); // for test purposes
	void fprint (char * file_name); // for test purposes
	void fprint_dif_files (char * folder); // for test purposes

	void add_to_entry (int I, int J, double value); // adds value to entry I, J (as in dense matrix)
	void set_entry (int I, int J, double value); // sets value to entry I, J (as in dense matrix)
	void add_to_f_entry (int I, double value); // adds value to f[I]
	void set_f_entry (int I, double value); // sets value to f[I]
	void multiply_entry (int I, int J, double value); // multiplies entry I, J (as in dense matrix) to value
	
	double get_entry (int I, int J); // returns value of the entry I, J (as in dense matrix)

	void get_solution (MathVector * solution); // returns x
	void set_starting_point (MathVector * start_x); // sets starting vector for iterator
	void clear_row (int iF); // clears iF-row

	void factorize_Cholesky (); // factorization, symmetrical
	void factorize (); // factorization of the matrix. Non-symmetrical

	void test_full (int Size, bool symm);
	void test_sparse (int Size, bool symm);

	int solve (int * param); // param: method, decomp_type, depth (for gmres), parallel, number of threads, mkl functions
};
