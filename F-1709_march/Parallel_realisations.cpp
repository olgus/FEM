#include "Parallel_realisations.h"

namespace P_BLAS
{
	void mult_A_v (int n, int * ig, int * jg, double * di, double * ggl, double * ggu, double * x, double * res, double * res_omp, int par_num_threads)
	{
		int i, j;
		int start, end, row;
		int rank, addr;

#pragma omp parallel shared (ig, jg, di, ggl, ggu, x, res, res_omp) private (i, j, start, end, row, rank, addr) num_threads (par_num_threads)
		{
#pragma omp for
			for (i = 0; i < n; i++)
				res[i] = di[i] * x[i];
#pragma omp for
			for (i = 0; i < n * (par_num_threads - 1); i++)
				res_omp[i] = 0.0;

#pragma omp for schedule (dynamic, BLOCK)
			for (i = 0; i < n; i++)
			{
				rank = omp_get_thread_num ();

				if (rank == 0)
				{
					start = ig[i];
					end = ig[i + 1];
					for (j = start; j < end; j++)
					{
						row = jg[j];
						//умножение на матрицу А
						res[i] += ggl[j] * x[row];
						res[row] += ggu[j] * x[i];
					}
				}
				else
				{
					addr = (rank - 1) * n;
					start = ig[i];
					end = ig[i + 1];
					for (j = start; j < end; j++)
					{
						row = jg[j];
						//умножение на матрицу А
						res_omp[addr + i] += ggl[j] * x[row];
						res_omp[addr + row] += ggu[j] * x[i];
					}
				}
			}
#pragma omp for
			for (i = 0; i < n; i++)
			{
				for (j = 0; j < par_num_threads - 1; j++)
				{
					res[i] += res_omp[j * n + i];
				}
			}
		}

	}

	double Scalar_Product (int n, double * a, double * b, int par_num_threads)
	{
		double r = 0.0;
#pragma omp parallel shared (a, b) num_threads(par_num_threads)
		{
#pragma omp for reduction (+:r)
			for (int i = 0; i < n; i++)
			{
				r += a[i] * b[i];
			}
		}
		return r;
	}

	void Linear_Combination (int n, double a, double * c, double b, double * d, int par_num_threads)
	{
#pragma omp parallel shared (c, d) num_threads(par_num_threads)
		{
#pragma omp for
			for (int i = 0; i < n; i++)
			{
				c[i] = a * c[i] + b * d[i];
			}
		}
	}

	void Linear_Combination (int n, double a, double * c, double b, double * d, double * f, int par_num_threads)
	{
#pragma omp parallel shared (c, d) num_threads(par_num_threads)
		{
#pragma omp for
			for (int i = 0; i < n; i++)
			{
				f[i] = a * c[i] + b * d[i];
			}
		}
	}

	void calc_DvR (int size, double * di_inv, double * x, double * res, int par_num_threads)
	{
#pragma omp parallel shared (size, di_inv, x, res) num_threads (par_num_threads)
#pragma omp for
		for (int i = 0; i < size; i++)
		{
			res[i] = x[i] * di_inv[i];
		}
	}
}