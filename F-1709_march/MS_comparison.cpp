#include "MS_comparison.h"

double compare_solutions (int k_system, Task_pointer * comp_task, Task_pointer * rel_task)
{
	// get solution from comp_task on a mesh for rel_task
	// build M_matrix and f from local matrixes

	compressed_matrix * C = new compressed_matrix ();
	(*C) = rel_task->eq_matrixes[k_system];
	 
	int dim = rel_task->mesh_pointer->get_dimentionality ();
	int n_functions = rel_task->mesh_pointer->get_amount_non_zero_functions (0);
	int n_functions_cur;
	Matrix * M = new Matrix (n_functions, n_functions);
	MathVector * B = new MathVector (n_functions);

	int iF, jF;
	for (int k_element = 0, k_element_end = rel_task->mesh_pointer->get_n_elements(); k_element < k_element_end; k_element++)
	{
		n_functions_cur = rel_task->mesh_pointer->get_amount_non_zero_functions (k_element);
		if (n_functions != n_functions_cur)
		{
			delete B;
			delete M;

			n_functions = n_functions_cur;

			M = new Matrix (n_functions, n_functions);
			B = new MathVector (n_functions);
		}

		// get local M
		rel_task->mesh_pointer->get_M_local_matrix (k_element, M);
		// build local B
		{
			double jac = 0.0; // jacobian
			int n_integr_points = rel_task->mesh_pointer->amount_of_integration_points (k_element);
			double * weigths = new double[n_integr_points];
			double ** points = new double *[n_integr_points];
			for (int i = 0; i < n_integr_points; i++)
			{
				points[i] = new double[dim];
			}

			rel_task->mesh_pointer->integration_points (k_element, points, weigths, &jac);

			double val;
			double value;
			for (int i = 0; i < n_functions; i++)
			{
				// get vector of integrals f * function_i
				val = 0.0;
				// get theta values * basis_functions in those points
				for (int k = 0; k < n_integr_points; k++)
				{
					comp_task->get_solution_in_point (k_system, points[k], &value);
					// sum them multiplying by weigths
					val += rel_task->mesh_pointer->get_basis_function_value (k_element, i, points[k]) * value * weigths[k];
				}
				// multiply by jac
				val *= jac;
				B->setElem (i, val);
			}

			delete[] weigths;
			for (int i = 0; i < n_integr_points; i++)
			{
				delete[] points[i];
			}
			delete[] points;
		}

		// go by element's functions
		for (int i = 0; i < n_functions; i++)
		{
			iF = rel_task->get_function_global_number (k_element, i);
			if (iF != -1)
			{
				for (int j = 0; j < n_functions; j++)
				{
					jF = rel_task->get_function_global_number (k_element, j);

					if (iF != -1)
					{
						// add M into respective places (by def_nodes)
						C->add_to_entry (iF, jF, M->Elem (i, j));
					}
				}

				// put B into respective places (by def_nodes)
				C->add_to_f_entry (iF, B->getElem (i));
			}
		}

		//// get local G
		//rel_task->mesh_pointer->get_G_local_matrix (k_element, M);
		//for (int i = 0; i < n_functions; i++)
		//{
		//	iF = rel_task->get_function_global_number (k_element, i);
		//	if (iF != -1)
		//	{
		//		for (int j = 0; j < n_functions; j++)
		//		{
		//			jF = rel_task->get_function_global_number (k_element, j);

		//			if (iF != -1)
		//			{
		//				// add M into respective places (by def_nodes)
		//				C->add_to_entry (iF, jF, 1e-11 * M->Elem (i, j));
		//			}
		//		}
		//	}
		//}
	}
	delete M;
	delete B;
	//C->fprint ();

	MathVector * QC = new MathVector (C->Size ());
	MathVector * QR = new MathVector (C->Size ());
	MathVector * Q = new MathVector (C->Size ());
	// solve it
	//C->solve_LOS (0);
	C->solve_pardiso (4);
	// save solution
	C->get_solution (QC);

	rel_task->get_full_solution (k_system, QR);
	// get a vector of difference between solutions
	Q->Copy (*QR);
	Q->Substract (*QC);
	// multiply M on it
	C->mult_A_v (*Q, QR);
	// multiply vectors
	double sp = QR->Scalar_Product (*Q);

	delete C;
	delete Q;
	delete QR;
	delete QC;

	return sp;
}
