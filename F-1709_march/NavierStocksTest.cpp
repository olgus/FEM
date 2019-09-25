#include "NavierStocksTest.h"

void NavierStocksTest::set_n_systems ()
{
	n_systems = 3;
	eq_matrixes = new compressed_matrix[n_systems];
	printf ("systems settled:\t%i\n", n_systems);
	fprintf (log_file, "systems settled:\t%i\n", n_systems);
}

void NavierStocksTest::build_system (int k_system)
{
	switch (k_system)
	{
	case 0:
	{
		Matrix * G;
		Matrix * M;
		Matrix * DDFx;
		Matrix * DDFy;
		MathVector * B;
		MathVector * Q;
		MathVector * MQ;
		MathVector * Time_coef;
		MathVector * Psi;

		int n_functions;
		int n_functions_cur;
		int iF, jF;
		int dim = mesh_pointer->get_dimentionality ();

		Time_coef = new MathVector (time_sampling);
		get_time_coefficients (Time_coef);

		n_functions = mesh_pointer->get_amount_non_zero_functions (0);
		G = new Matrix (n_functions, n_functions);
		M = new Matrix (n_functions, n_functions);
		B = new MathVector (n_functions);
		Psi = new MathVector (n_functions);
		DDFx = new Matrix (n_functions, n_functions);
		DDFy = new Matrix (n_functions, n_functions);
		Q = new MathVector (n_functions);
		MQ = new MathVector (n_functions);

		for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
		{
			n_functions_cur = mesh_pointer->get_amount_non_zero_functions (k_element);
			if (n_functions != n_functions_cur)
			{
				delete B;
				delete G;
				delete M;
				delete Q;
				delete MQ;
				delete Psi;
				delete DDFx;
				delete DDFy;

				n_functions = n_functions_cur;

				G = new Matrix (n_functions, n_functions);
				M = new Matrix (n_functions, n_functions);
				B = new MathVector (n_functions);
				Psi = new MathVector (n_functions);
				DDFx = new Matrix (n_functions, n_functions);
				DDFy = new Matrix (n_functions, n_functions);
				Q = new MathVector (n_functions);
				MQ = new MathVector (n_functions);
			}

			// get G, M matrices 
			mesh_pointer->get_G_local_matrix (k_element, G);
			mesh_pointer->get_M_local_matrix (k_element, M);

			// get B-vector
			get_local_B (k_system, k_element, B);
			
			// get psi values
			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (k_element, i);
				if (iF != -1)
				{
					Psi->setElem (i, previous_time_layers_solutions[2][0].getElem (iF));
				}
			}
			// go by element's functions
			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (k_element, i);
				if (iF != -1)
				{
					for (int j = 0; j < n_functions; j++)
					{
						jF = get_function_global_number (k_element, j);
						if (jF != -1)
						{
							// add G into respective places (by def_nodes)
							eq_matrixes[k_system].add_to_entry (iF, jF, G->Elem (i, j));
							// add M (by t) into respective places (by def_nodes)
							eq_matrixes[k_system].add_to_entry (iF, jF, M->Elem (i, j) * Time_coef->getElem (0));
						}
					}

					// put B into respective places (by def_nodes)
					eq_matrixes[k_system].add_to_f_entry (iF, B->getElem (i));
				}
			}

			// add derivatives
			for (int k = 0; k < n_functions; k++)
			{
				mesh. elements[k_element]->get_DDF (k, 1, 0, DDFy);
				mesh. elements[k_element]->get_DDF (k, 0, 1, DDFx);

				for (int i = 0; i < n_functions; i++)
				{
					iF = get_function_global_number (k_element, i);

					if (iF != -1)
					{
						for (int j = 0; j < n_functions; j++)
						{
							jF = get_function_global_number (k_element, j);
							if (jF != -1)
							{
								eq_matrixes[k_system].add_to_entry (iF, jF, -DDFx->Elem (i, j) * Psi->getElem (k));
								eq_matrixes[k_system].add_to_entry (iF, jF, DDFy->Elem (i, j) * Psi->getElem (k));
							}
						}
					}
				}
			}

			// go by amount of time layers that count
			for (int t = 1; t < time_sampling; t++)
			{
				// get q's
				for (int i = 0; i < n_functions; i++)
				{
					iF = get_function_global_number (k_element, i);
					if (iF != -1)
					{
						Q->setElem (i, previous_time_layers_solutions[k_system][t - 1].getElem (iF));
					}
				}
				M->MultiplyMatrixByVector (*Q, MQ);
				// go by element's functions
				for (int i = 0; i < n_functions; i++)
				{
					iF = get_function_global_number (k_element, i);
					if (iF != -1)
					{
						eq_matrixes[k_system].add_to_f_entry (iF, MQ->getElem (i) * Time_coef->getElem (t));
					}
				}
			}

		}

		delete Time_coef;
		delete B;
		delete G;
		delete M;
		delete Q;
		delete MQ;
		delete Psi;
		delete DDFx;
		delete DDFy;
		break;
	}
	case 1:
	{
		Matrix * G;
		Matrix * M;
		Matrix * DDFx;
		Matrix * DDFy;
		Matrix * D;
		MathVector * B;
		MathVector * Q;
		MathVector * MQ;
		MathVector * Time_coef;
		MathVector * Psi;
		MathVector * T;
		MathVector * TDF;

		int n_functions;
		int n_functions_cur;
		int iF, jF;
		int dim = mesh_pointer->get_dimentionality ();

		Time_coef = new MathVector (time_sampling);
		get_time_coefficients (Time_coef);

		n_functions = mesh_pointer->get_amount_non_zero_functions (0);
		G = new Matrix (n_functions, n_functions);
		M = new Matrix (n_functions, n_functions);
		B = new MathVector (n_functions);
		Psi = new MathVector (n_functions);
		DDFx = new Matrix (n_functions, n_functions);
		DDFy = new Matrix (n_functions, n_functions);
		Q = new MathVector (n_functions);
		MQ = new MathVector (n_functions);
		T = new MathVector (n_functions);
		TDF = new MathVector (n_functions);
		D = new Matrix (n_functions, n_functions);

		for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
		{
			n_functions_cur = mesh_pointer->get_amount_non_zero_functions (k_element);
			if (n_functions != n_functions_cur)
			{
				delete B;
				delete G;
				delete M;
				delete Q;
				delete MQ;
				delete Psi;
				delete DDFx;
				delete DDFy;
				delete T;
				delete TDF;
				delete D;

				n_functions = n_functions_cur;

				G = new Matrix (n_functions, n_functions);
				M = new Matrix (n_functions, n_functions);
				B = new MathVector (n_functions);
				Psi = new MathVector (n_functions);
				DDFx = new Matrix (n_functions, n_functions);
				DDFy = new Matrix (n_functions, n_functions);
				Q = new MathVector (n_functions);
				MQ = new MathVector (n_functions);
				T = new MathVector (n_functions);
				TDF = new MathVector (n_functions);
				D = new Matrix (n_functions, n_functions);
			}
			// get G, M matrices 
			mesh_pointer->get_G_local_matrix (k_element, G);
			mesh_pointer->get_M_local_matrix (k_element, M);

			// get B-vector
			get_local_B (k_system, k_element, B);

			// get psi values
			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (k_element, i);
				if (iF != -1)
				{
					Psi->setElem (i, previous_time_layers_solutions[2][0].getElem (iF));
				}
			}

			// get T values
			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (k_element, i);
				if (iF != -1)
				{
					T->setElem (i, previous_time_layers_solutions[0][0].getElem (iF));
				}
			}

			// get dT/dx
			mesh. elements[k_element]->get_D (0, D);
			D->MultiplyMatrixByVector (*T, TDF);

			// go by element's functions
			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (k_element, i);
				if (iF != -1)
				{
					for (int j = 0; j < n_functions; j++)
					{
						jF = get_function_global_number (k_element, j);
						if (jF != -1)
						{
							// add G into respective places (by def_nodes)
							eq_matrixes[k_system].add_to_entry (iF, jF, G->Elem (i, j));
							// add M (by t) into respective places (by def_nodes)
							eq_matrixes[k_system].add_to_entry (iF, jF, M->Elem (i, j) * Time_coef->getElem (0));
						}
					}

					// put B into respective places (by def_nodes)
					eq_matrixes[k_system].add_to_f_entry (iF, B->getElem (i));

					// add dt/dx
					eq_matrixes[k_system].add_to_f_entry (iF, TDF->getElem (i));
				}
			}

			// add derivatives
			for (int k = 0; k < n_functions; k++)
			{
				mesh. elements[k_element]->get_DDF (k, 1, 0, DDFy);
				mesh. elements[k_element]->get_DDF (k, 0, 1, DDFx);

				for (int i = 0; i < n_functions; i++)
				{
					iF = get_function_global_number (k_element, i);
					if (iF != -1)
					{
						for (int j = 0; j < n_functions; j++)
						{
							jF = get_function_global_number (k_element, j);
							if (jF != -1)
							{
								eq_matrixes[k_system].add_to_entry (iF, jF, -DDFx->Elem (i, j) * Psi->getElem (k));
								eq_matrixes[k_system].add_to_entry (iF, jF, DDFy->Elem (i, j) * Psi->getElem (k));
							}
						}
					}
				}
			}

			// go by amount of time layers that count
			for (int t = 1; t < time_sampling; t++)
			{
				// get q's
				for (int i = 0; i < n_functions; i++)
				{
					iF = get_function_global_number (k_element, i);
					if (iF != -1)
					{
						Q->setElem (i, previous_time_layers_solutions[k_system][t - 1].getElem (iF));
					}
				}
				M->MultiplyMatrixByVector (*Q, MQ);
				// go by element's functions
				for (int i = 0; i < n_functions; i++)
				{
					iF = get_function_global_number (k_element, i);
					if (iF != -1)
					{
						eq_matrixes[k_system].add_to_f_entry (iF, MQ->getElem (i) * Time_coef->getElem (t));
					}
				}
			}

		}

		delete Time_coef;
		delete B;
		delete G;
		delete M;
		delete Q;
		delete MQ;
		delete Psi;
		delete DDFx;
		delete DDFy;
		delete T;
		delete TDF;
		delete D;
		break;
	}
	case 2:
	{
		Matrix * G;
		Matrix * M;
		MathVector * B;
		MathVector * W;
		MathVector * WM;
		int n_functions;
		int n_functions_cur;
		int iF, jF;
		int dim = mesh_pointer->get_dimentionality ();

		n_functions = mesh_pointer->get_amount_non_zero_functions (0);
		G = new Matrix (n_functions, n_functions);
		M = new Matrix (n_functions, n_functions);
		B = new MathVector (n_functions);
		W = new MathVector (n_functions);
		WM = new MathVector (n_functions);
		
		for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
		{
			n_functions_cur = mesh_pointer->get_amount_non_zero_functions (k_element);
			if (n_functions != n_functions_cur)
			{
				delete B;
				delete G;
				delete M;
				delete W;
				delete WM;

				n_functions = n_functions_cur;

				G = new Matrix (n_functions, n_functions);
				M = new Matrix (n_functions, n_functions);
				B = new MathVector (n_functions);
				W = new MathVector (n_functions);
				WM = new MathVector (n_functions);
			}
			// get G, M matrices 
			mesh_pointer->get_G_local_matrix (k_element, G);
			mesh_pointer->get_M_local_matrix (k_element, M);
			
			// get B-vector
			get_local_B (k_system, k_element, B);

			// get W
			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (k_element, i);
				if (iF != -1)
				{
					W->setElem (i, previous_time_layers_solutions[1][0].getElem (iF));
				}
			}
			M->MultiplyMatrixByVector (*W, WM);

			// go by element's functions
			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (k_element, i);
				if (iF != -1)
				{
					for (int j = 0; j < n_functions; j++)
					{
						jF = get_function_global_number (k_element, j);
						if (jF != -1)
						{
							// add G into respective places (by def_nodes)
							eq_matrixes[k_system].add_to_entry (iF, jF, G->Elem (i, j));
						}
					}

					// put B into respective places (by def_nodes)
					eq_matrixes[k_system].add_to_f_entry (iF, B->getElem (i));
					// put W into respective places (by def_nodes)
					eq_matrixes[k_system].add_to_f_entry (iF, WM->getElem (i));
				}
			}

		}
		delete B;
		delete G;
		delete M;
		delete W;
		delete WM;
		break;

	}
	}
}

double NavierStocksTest::function_f (int k_system, double * coordinates, int area)
{
	double r = 0.0;
	double x = coordinates[0];
	double y = coordinates[1];
	double t = time_layers[current_time_layer];

	switch (k_system)
	{
	case 0:
	{
		r += exp (x + y + t) * ((3.0 * pow (x, 3.0) * pow (y, 2.0)) - (3.0 * pow (x, 2.0) * pow (y, 3.0)) - 1.0);
		//r += 2.0;
		break;
	}
	case 1:
	{
		r += sin (x + y + t) * ((3.0 * pow (x, 2.0) * pow (y, 3.0)) - (3.0 * pow (x, 3.0) * pow (y, 2.0)) - 1.0);
		r += 2.0 * cos (x + y + t);
		r -= exp (x + y + t);
		//r += - x - y - t - 1.0;
		break;
	}
	case 2:
	{
		r -= cos (x + y + t);
		r -= ((6.0 * x * pow (y, 3.0)) + (6.0 * pow (x, 3.0) * y));
		//r += (x + y) * t;
		break;
	}
	}
	return r;
}

double NavierStocksTest::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	double r = 0.0;
	double x = coordinates[0];
	double y = coordinates[1];
	double t = time_layers[current_time_layer];

	switch (k_system)
	{
	case 0:
	{
		r += exp (x + y + t);
		//r += x + y + t;
		break;
	}
	case 1:
	{
		r += cos (x + y + t);
		//r += -(x + y) * t;
		break;
	}
	case 2:
	{
		r += pow (x, 3.0) * pow (y, 3.0);
		//r += x * x + y * y;
		//r = cos (x + y);
		//r += x + 2.0 * y + 3.0 * t;
		break;
	}
	}
	return r;
}

NavierStocksTest::NavierStocksTest ()
{
}

NavierStocksTest::NavierStocksTest (const NavierStocksTest & nst)
{
}

NavierStocksTest::~NavierStocksTest ()
{
}

//void NavierStocks::set_starting_conditions ()
//{
//	for (int i = 0; i < time_sampling; i++)
//	{
//		previous_time_layers_solutions[0][i].setSize (eq_matrixes[0].Size ()); // default starting solution point is 0
//		previous_time_layers_solutions[0][i].Initialize_const (phys_T0);
//	}
//
//	for (int k_system = 1; k_system < n_systems; k_system++)
//	{
//		for (int i = 0; i < time_sampling; i++)
//		{
//			previous_time_layers_solutions[k_system][i].setSize (eq_matrixes[k_system].Size ()); // default starting solution point is 0
//		}
//	}
//	printf ("\tstarting conditions are settled\n");
//}
