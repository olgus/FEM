#pragma once
#include <omp.h>

#define BLOCK 100

namespace P_BLAS
{
	void mult_A_v (int n, int * ig, int * jg, double * di, double * ggl, double * ggu, double * x, double * res, double * res_omp, int par_num_threads);
	double Scalar_Product (int n, double * a, double * b, int par_num_threads); // return scalar product of this and c vectors, parallel
	void Linear_Combination (int n, double a, double * c, double b, double * d, int par_num_threads); // add the vector b * c to current with a modifiera
	void Linear_Combination (int n, double a, double * c, double b, double * d, double * f, int par_num_threads); // add the vector b * c to current with a modifiera
	void calc_DvR (int size, double * di_inv, double * x, double * res, int par_num_threads); // calc v = D-1 R 
}