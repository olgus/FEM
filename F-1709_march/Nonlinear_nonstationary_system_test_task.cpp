#include "Nonlinear_nonstationary_system_test_task.h"

bool Nonlinear_nonstationary_system_task::solve_task (int solver_param[][5])
{
	{
		if (!prepared)
		{
			printf ("ERROR: missing preparation\n");
			return false;
		}

		// n_systems is set in preparation to 1 by default
		non_linear_layers_solutions = new MathVector[n_systems];
		MathVector * prev_non_linear_layers_solutions = new MathVector (eq_matrixes[0].Size ());
		for (int k_system = 0; k_system < n_systems; k_system++)
		{
			non_linear_layers_solutions[k_system].setSize (eq_matrixes[k_system].Size ()); // default starting solution point is 0
		}

		double sum_non_linear_dif;
		double * non_linear_dif = new double[n_systems];
		double w_rel = 1.0;

		// 0 non_linear solution is a start solution
		set_starting_conditions ();

		int solver_iterations;
		bool found = false;

		// TODO time cycle
		bool stable = false; // indicator of reaching stable condition
		for (current_time_layer = 1; current_time_layer < n_time_layers && !stable; current_time_layer++)
		{
			printf ("time layer: %.7lf\n", time_layers[current_time_layer]);
			found = false;
			for (int k_system = 0; k_system < n_systems; k_system++)
			{
				non_linear_layers_solutions[k_system].Copy ((previous_time_layers_solutions[k_system][0]));
			}

			for (int k_nl_iter = 0; k_nl_iter < MAX_NONLINEAR_ITER && !found; k_nl_iter++)
			{
				if ((k_nl_iter + 2) % 8 == 0)
				{
					system ("cls");
					printf ("time layer: %.7lf\n", time_layers[current_time_layer]);
				}

				printf ("k_nl_iter: %i\n", k_nl_iter);
				sum_non_linear_dif = 0.0;

				for (int k_system = 0; k_system < n_systems; k_system++)
				{
					// save prev non linear solution
					prev_non_linear_layers_solutions->Copy (non_linear_layers_solutions[k_system]);
					// set starting point for the solver
					eq_matrixes[k_system].set_starting_point (&(non_linear_layers_solutions[k_system])); // set last layer solution as x0 for solver
																										 // build system
					build_system (k_system);
					// apply boundary conditions
					apply_boundary_conditions (k_system); // apply boundary conditions
														  // solve it
					solver_iterations = eq_matrixes[k_system].solve (solver_param[k_system]);

					// TODO relax
					// for now just copy the solution from the matrix
					//eq_matrixes[k_system].get_solution (&(non_linear_layers_solutions[k_system]));
					//relax (k_system, 1.0, *prev_non_linear_layers_solutions);
					//w_rel = relax (k_system, *prev_non_linear_layers_solutions);
					w_rel = relax (k_system, 0.1, 1.0, *prev_non_linear_layers_solutions);

					// non_linear_layers_solutions contains entire _new_ solution 
					prev_non_linear_layers_solutions->Substract (non_linear_layers_solutions[k_system]);
					non_linear_dif[k_system] = prev_non_linear_layers_solutions->Norm () / non_linear_layers_solutions[k_system].Norm ();

					printf ("\t\tk_system: %i\tdif: %.3e\tw_rel: %.5lf\n", k_system, non_linear_dif[k_system], w_rel);
					sum_non_linear_dif += non_linear_dif[k_system];

					print_solutions (k_nl_iter);
				}
				if ((k_nl_iter > 0) && (sum_non_linear_dif < NLNSSTT_NL_PRECIS))
					found = true;
			}

			// move time layers
			for (int k_system = 0; k_system < n_systems; k_system++)
			{
				for (int j = time_sampling - 1; j > 0; j--)
				{
					previous_time_layers_solutions[k_system][j] = previous_time_layers_solutions[k_system][j - 1];
				}
				previous_time_layers_solutions[k_system][0] = non_linear_layers_solutions[k_system];
			}
			print_solutions ();
		}

		delete[] non_linear_dif;
		delete prev_non_linear_layers_solutions;
		delete[] non_linear_layers_solutions;
		non_linear_layers_solutions = NULL;
		return true;
	}
}

double Nonlinear_nonstationary_system_task::relax (int k_system, double w, const MathVector & prev_non_linear_layer_solution)
{
	eq_matrixes[k_system].get_solution (&(non_linear_layers_solutions[k_system])); // u k
	non_linear_layers_solutions[k_system].Linear_Combination (w, prev_non_linear_layer_solution, 1.0 - w);
	return 0.0;
}

double Nonlinear_nonstationary_system_task::relax (int k_system, double w_left, double w_right, const MathVector & prev_non_linear_layer_solution)
{
	double w = w_right;
	double f_prev, f;
	int i_pos = 0;
	MathVector * sol = new MathVector (eq_matrixes[k_system].Size ()); // matrix solution
	eq_matrixes[k_system].get_solution (sol);
	MathVector * mv = new MathVector (eq_matrixes[k_system].Size ());

	int par_num_threads = 4;
	double * res_omp = NULL;
	if (par_num_threads > 1)
		res_omp = new double[(par_num_threads - 1) * eq_matrixes[k_system].Size ()];

	f_prev = 1e+10;
	bool found = false;
	for (int i = 0; i < 10 && w > w_left && !found; i++)
	{
		non_linear_layers_solutions[k_system].Linear_Combination (*sol, w, prev_non_linear_layer_solution, 1.0 - w);
		build_system (k_system);
		apply_boundary_conditions (k_system);
		eq_matrixes[k_system].mult_A_v (non_linear_layers_solutions[k_system], mv, res_omp, par_num_threads);
		mv->Substract (*(eq_matrixes[k_system].f));
		f = mv->Norm ();

		if (f > f_prev)
		{
			i_pos = i - 1;
			found = true;
		}
		f_prev = f;
		w = w / 1.25;
		//printf ("%i", i);
	}
	if (i_pos == -1)
		i_pos = 0;

	w = w_right / pow (1.25, i_pos);
	non_linear_layers_solutions[k_system].Linear_Combination (*sol, w, prev_non_linear_layer_solution, 1.0 - w);

	delete sol;
	if (par_num_threads > 1)
		delete[] res_omp;
	delete mv;

	return w;
}

double Nonlinear_nonstationary_system_task::relax (int k_system, const MathVector & prev_non_linear_layer_solution)
{
	double ms;
	clock_t t1, t2;

	double c = (1.0 + sqrt (5.0)) / 2.0; // coef
	MathVector * sol = new MathVector (eq_matrixes[k_system].Size ()); // matrix solution
	MathVector * mv = new MathVector (eq_matrixes[k_system].Size ());
	eq_matrixes[k_system].get_solution (sol);

	// w boundaries
	double w;
	double a = 0.0;
	double b = 1.0;
	// current boundaries
	double x1 = b - (b - a) / c;
	double x2 = a + (b - a) / c;

	double y1, y2; // function values at xs

	int par_num_threads = 1;
	double * res_omp = NULL;
	if (par_num_threads > 1)
		res_omp = new double[(par_num_threads - 1) * eq_matrixes[k_system].Size ()];

	// left
	non_linear_layers_solutions[k_system].Linear_Combination (*sol, x1, prev_non_linear_layer_solution, 1.0 - x1);
	build_system (k_system);
	apply_boundary_conditions (k_system);
	eq_matrixes[k_system].mult_A_v (non_linear_layers_solutions[k_system], mv, res_omp, par_num_threads);
	mv->Substract (*(eq_matrixes[k_system].f));
	y1 = mv->Norm ();
	//printf ("%lf %.3e\n", x1, y1);
	// right
	non_linear_layers_solutions[k_system].Linear_Combination (*sol, x2, prev_non_linear_layer_solution, 1.0 - x2);
	build_system (k_system);
	apply_boundary_conditions (k_system);
	eq_matrixes[k_system].mult_A_v (non_linear_layers_solutions[k_system], mv, res_omp, par_num_threads);
	mv->Substract (*(eq_matrixes[k_system].f));
	y2 = mv->Norm ();
	//printf ("%lf %.3e\n", x2, y2);

	for (int i = 0; (i < 100) && (fabs (b - a) > 2e-2); i++)
	{
		//printf ("%i", i);
		if (y1 >= y2)
		{
			a = x1;
			x1 = x2;
			y1 = y2;

			x2 = a + (b - a) / c;
			non_linear_layers_solutions[k_system].Linear_Combination (*sol, x2, prev_non_linear_layer_solution, 1.0 - x2);

			//t1 = clock ();

			build_system (k_system);

			//t2 = clock ();
			//ms = (double)(t2 - t1);
			//printf ("time spend building system: %lf s\n", ms / 1000);

			//t1 = clock ();

			apply_boundary_conditions (k_system);

			//t2 = clock ();
			//ms = (double)(t2 - t1);
			//printf ("time spend applying boundary cond: %lf s\n", ms / 1000);

			//t1 = clock ();

			eq_matrixes[k_system].mult_A_v (non_linear_layers_solutions[k_system], mv, res_omp, par_num_threads);

			//t2 = clock ();
			//ms = (double)(t2 - t1);
			//printf ("time spend multiplying: %lf s\n", ms / 1000);

			mv->Substract (*(eq_matrixes[k_system].f));
			y2 = mv->Norm ();
			//printf ("%lf %.3e\n", x2, y2);
		}
		else
		{
			b = x2;
			x2 = x1;
			y2 = y1;

			x1 = b - (b - a) / c;
			non_linear_layers_solutions[k_system].Linear_Combination (*sol, x1, prev_non_linear_layer_solution, 1.0 - x1);
			//t1 = clock ();

			build_system (k_system);

			//t2 = clock ();
			//ms = (double)(t2 - t1);
			//printf ("time spend building system: %lf s\n", ms / 1000);

			//t1 = clock ();

			apply_boundary_conditions (k_system);

			//t2 = clock ();
			//ms = (double)(t2 - t1);
			//printf ("time spend applying boundary cond: %lf s\n", ms / 1000);

			//t1 = clock ();

			eq_matrixes[k_system].mult_A_v (non_linear_layers_solutions[k_system], mv, res_omp, par_num_threads);

			//t2 = clock ();
			//ms = (double)(t2 - t1);
			//printf ("time spend multiplying: %lf s\n", ms / 1000);

			mv->Substract (*(eq_matrixes[k_system].f));
			y1 = mv->Norm ();
			//printf ("%lf %.3e\n", x1, y1);
		}
	}


	// test purposes
	//{
	//	FILE * file = fopen ("Test//rel.txt", "w");
	//	double x;
	//	for (int i = 0; i < 101; i++)
	//	{
	//		printf ("%i ", i);
	//		x = i * 0.01;
	//		non_linear_layers_solutions[k_system].Linear_Combination (*sol, x, prev_non_linear_layer_solution, 1.0 - x);
	//		build_system (k_system);
	//		apply_boundary_conditions (k_system);
	//		eq_matrixes[k_system].mult_A_v (non_linear_layers_solutions[k_system], mv);
	//		mv->Substract (*(eq_matrixes[k_system].f));
	//		y1 = mv->Norm ();
	//		fprintf (file, "%lf %.13lf\n", x, y1);
	//	}
	//	fclose (file);
	//	printf ("\n");
	//}


	w = (a + b) / 2.0;
	if (w > 0.99)
		w = 1.0;

	non_linear_layers_solutions[k_system].Linear_Combination (*sol, w, prev_non_linear_layer_solution, 1.0 - w);

	delete sol;
	if (par_num_threads > 1)
		delete[] res_omp;
	delete mv;

	return w;
}

double Nonlinear_nonstationary_system_task::get_function_value (int k_system, int k_function, double * coordinates, int area)
{
	double r = 0.0;
	switch (k_function)
	{
	case 0:
		// f function
		r = function_f (k_system, coordinates, area);
		break;
	}
	return r;
}

double Nonlinear_nonstationary_system_task::function_f (int k_system, double * coordinates, int area)
{
	double t = time_layers[current_time_layer];
	double x = coordinates[0];
	double y = coordinates[1];

	double r = 0.0;
	switch (k_system)
	{
	//case 0:
	//	r = 1.0;
	//	break;
	//case 1:
	//	r = -1.0;
	//	break;
	//case 2:
	//	r = x + y - t;
	//	break;

	//case 0:
	//	r = -3.0 - t;
	//	break;
	//case 1:
	//	r = t + 1.0;
	//	break;
	//case 2:
	//	r = x + y - t; 
	//	break;

	//case 0:
	//	r = x * x + y * y + exp (x + y + t) * 2.0 * t * (x - y) - 4.0 * t;
	//	break;
	//case 1:
	//	r = x * x - 2.0 * y * y + exp (x + y + t) * 2.0 * t * (x + 2.0 * y) + 4.0 * t - 2.0 * x * t;
	//	break;
	//case 2:
	//	r = -2.0 * exp (x + y + t) - x * x * t + 2.0 * t * y * y;
	//	break;

	case 0:
		r = x * x + y * y - 4.0 * t + 4.0 * x * x * y * t * t - 2.0 * y * y * y * t * t;
		break;
	case 1:
		r = 5.0 - 2.0 * x * t - 4.0 * x * x * y * t + 2.0 * y * y * y * t;
		break;
	case 2:
		r = -2.0 * x * t + x * x + y * y - t;
		break;
	}
	return r;
}

double Nonlinear_nonstationary_system_task::function_starting_condition (int k_system, double * coordinates, int area)
{
	return function_FCondition (k_system, coordinates, area, 0);
}

double Nonlinear_nonstationary_system_task::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	double x = coordinates[0];
	double y = coordinates[1];
	double t = time_layers[current_time_layer];

	double r = 0.0;
	switch (k_system)
	{
	//case 0:
	//	r = 2.0 * x + y + t;
	//	break;
	//case 1:
	//	r = - x - y + t;
	//	break;
	//case 2:
	//	r = x * t - 2.0 * y;
	//	break;

	//case 0:
	//	r = (x * x + y * y) * t;
	//	break;
	//case 1:
	//	r = x * x * t - 2.0 * t * y * y;
	//	break;
	//case 2:
	//	r = exp (x + y + t);
	//	break;

	case 0:
		r = (x * x + y * y) * t;
		break;
	case 1:
		r = - x * x - y * y + t;
		break;
	case 2:
		r = x * y * y * t;
		break;
	}
	return r;
}

void Nonlinear_nonstationary_system_task::get_function_coefficients (int k_system, int function, MathVector * coef)
{
	double coord[2];
	int area;
	for (int i = 0; i < mesh_pointer->get_n_nodes (); i++)
	{
		mesh_pointer->get_node_coordinates (i, coord);
		area = mesh_pointer->get_area (coord);
		coef->setElem (i, get_function_value (k_system, function, coord, area));
	}
}

void Nonlinear_nonstationary_system_task::build_system (int k_system)
{
	eq_matrixes[k_system].Clear ();
	int n_functions	= mesh_pointer->get_amount_non_zero_functions (0);
	int iF, jF, kF;

	switch (k_system)
	{
	case 0:
	{
		Matrix * G;
		Matrix * M;
		MathVector * DDx;
		MathVector * DDy;
		MathVector * Time_coef;
		MathVector * T_prev;
		MathVector * Psi;
		MathVector * F; // TEST

		G = new Matrix (n_functions, n_functions);
		M = new Matrix (n_functions, n_functions);
		DDx = new MathVector (n_functions * n_functions * n_functions);
		DDy = new MathVector (n_functions * n_functions * n_functions);
		T_prev = new MathVector (n_functions);
		Psi = new MathVector (n_functions);
		F = new MathVector (eq_matrixes[k_system].Size ()); // TEST

		Time_coef = new MathVector (time_sampling);
		get_time_coefficients (Time_coef);

		// get F cooeficients
		get_function_coefficients (k_system, 0, F);

		for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
		{
			// get G matrix
			mesh_pointer->get_G_local_matrix (k_element, G);
			// get M matrix
			mesh_pointer->get_M_local_matrix (k_element, M);
			// build Psi
			for (int i = 0; i < n_functions; i++)
			{
				Psi->setElem (i, non_linear_layers_solutions[2].getElem (get_function_global_number (k_element, i)));
			}
			// get DDs
			mesh.elements[k_element]->get_DDF (1, 0, DDx);
			mesh.elements[k_element]->get_DDF (0, 1, DDy);

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
							// add G
							eq_matrixes[k_system].add_to_entry (iF, jF, G->Elem (i, j));

							// add M
							eq_matrixes[k_system].add_to_entry (iF, jF, Time_coef->getElem (0) * M->Elem (i, j));
						
							// add and substract partial derivatives
							for (int k = 0; k < n_functions; k++)
							{
								eq_matrixes[k_system].add_to_entry (iF, jF, DDx->getElem (i * n_functions * n_functions + j * n_functions + k) * Psi->getElem (k));
								eq_matrixes[k_system].add_to_entry (iF, jF, -DDy->getElem (i * n_functions * n_functions + j * n_functions + k) * Psi->getElem (k));
							}
						}
					}

					// add f TEST
					for (int j = 0; j < n_functions; j++)
					{
						eq_matrixes[k_system].add_to_f_entry (iF, F->getElem (get_function_global_number (k_element, j)) * M->Elem (i, j));
					}
				}
			}

			// add prev time layers data
			{
				// go by amount of time layers that count
				for (int t = 1; t < time_sampling; t++)
				{
					// get previous time layers solutions
					for (int i = 0; i < n_functions; i++)
					{
						iF = get_function_global_number (k_element, i);
						if (iF != -1)
						{
							T_prev->setElem (i, previous_time_layers_solutions[k_system][t - 1].getElem (iF));
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
								eq_matrixes[k_system].add_to_f_entry (iF, T_prev->getElem (j) * Time_coef->getElem (t) * M->Elem (i, j));
							}
						}
					}
				}
			}
		}

		delete G;
		delete M;
		delete DDx;
		delete DDy;
		delete Time_coef;
		delete T_prev;
		delete Psi;
		delete F; // TEST
		break;
	}
	case 1:
	{
		Matrix * G;
		Matrix * M;
		Matrix * D;
		MathVector * DDx;
		MathVector * DDy;
		MathVector * Time_coef;
		MathVector * W_prev;
		MathVector * T;
		MathVector * F; // TEST
		MathVector * Psi;

		G = new Matrix (n_functions, n_functions);
		M = new Matrix (n_functions, n_functions);
		D = new Matrix (n_functions, n_functions);
		DDx = new MathVector (n_functions * n_functions * n_functions);
		DDy = new MathVector (n_functions * n_functions * n_functions);
		W_prev = new MathVector (n_functions);
		Psi = new MathVector (n_functions);
		T = new MathVector (n_functions);
		F = new MathVector (eq_matrixes[k_system].Size ()); // TEST

		Time_coef = new MathVector (time_sampling);
		get_time_coefficients (Time_coef);

		// get F cooeficients
		get_function_coefficients (k_system, 0, F);

		for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
		{
			// get G matrix
			mesh_pointer->get_G_local_matrix (k_element, G);
			// get M matrix
			mesh_pointer->get_M_local_matrix (k_element, M);
			// get bf * dbf/dx
			mesh.elements[k_element]->get_D (0, D);

			// build Psi
			for (int i = 0; i < n_functions; i++)
			{
				Psi->setElem (i, non_linear_layers_solutions[2].getElem (get_function_global_number (k_element, i)));
			}
			// build T
			for (int i = 0; i < n_functions; i++)
			{
				T->setElem (i, non_linear_layers_solutions[0].getElem (get_function_global_number (k_element, i)));
			}

			// get DDs
			mesh.elements[k_element]->get_DDF (1, 0, DDx);
			mesh.elements[k_element]->get_DDF (0, 1, DDy);

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
							// add G
							eq_matrixes[k_system].add_to_entry (iF, jF, G->Elem (i, j));

							// add M
							eq_matrixes[k_system].add_to_entry (iF, jF, Time_coef->getElem (0) * M->Elem (i, j));
						
							// add D (dT/dx) into f
							eq_matrixes[k_system].add_to_f_entry (iF, D->Elem (i, j) * T->getElem (j));

							// add and substract partial derivatives
							for (int k = 0; k < n_functions; k++)
							{
								eq_matrixes[k_system].add_to_entry (iF, jF, DDx->getElem (i * n_functions * n_functions + j * n_functions + k) * Psi->getElem (k));
								eq_matrixes[k_system].add_to_entry (iF, jF, -DDy->getElem (i * n_functions * n_functions + j * n_functions + k) * Psi->getElem (k));
							}					
						}
					}

					// add f TEST
					for (int j = 0; j < n_functions; j++)
					{
						eq_matrixes[k_system].add_to_f_entry (iF, F->getElem (get_function_global_number (k_element, j)) * M->Elem (i, j));
					}
				}
			}

			// add prev time layers data
			{
				// go by amount of time layers that count
				for (int t = 1; t < time_sampling; t++)
				{
					// get previous time layers solutions
					for (int i = 0; i < n_functions; i++)
					{
						iF = get_function_global_number (k_element, i);
						if (iF != -1)
						{
							W_prev->setElem (i, previous_time_layers_solutions[k_system][t - 1].getElem (iF));
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
								eq_matrixes[k_system].add_to_f_entry (iF, W_prev->getElem (j) * Time_coef->getElem (t) * M->Elem (i, j));
							}
						}
					}
				}
			}
		}

		delete G;
		delete M;
		delete D;
		delete DDx;
		delete DDy;
		delete Time_coef;
		delete W_prev;
		delete Psi;
		delete T;
		delete F; // TEST
		break;
	}
	case 2:
	{
		Matrix * G;
		Matrix * M;
		MathVector * W;
		MathVector * F; // TEST

		G = new Matrix (n_functions, n_functions);
		M = new Matrix (n_functions, n_functions);
		W = new MathVector (n_functions);
		F = new MathVector (eq_matrixes[k_system].Size ()); // TEST

		// get F cooeficients
		get_function_coefficients (k_system, 0, F);
		
		for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
		{
			// get G matrix
			mesh_pointer->get_G_local_matrix (k_element, G);
			// get M matrix
			mesh_pointer->get_M_local_matrix (k_element, M);
			// get W vector
			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (k_element, i);
				if (iF != -1)
				{
					W->setElem (i, non_linear_layers_solutions[1].getElem (iF));
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

							// add W * M into f
							eq_matrixes[k_system].add_to_f_entry (iF, W->getElem (j) * M->Elem (i, j));
						}
					}

					// add f TEST
					for (int j = 0; j < n_functions; j++)
					{
						eq_matrixes[k_system].add_to_f_entry (iF, F->getElem (get_function_global_number(k_element, j))  * M->Elem (i, j));
					}
				}
			}

		} 

		delete G;
		delete M;
		delete W;
		delete F; // TEST
		break;
	}
	}
}

void Nonlinear_nonstationary_system_task::print_solutions (int knliter)
{
	// picture for every non-linear solution 
	std::vector<double> local_time_stamps;
	bool save_solution = false;

	if (use_time_mapping)
	{
		// save solution in current time
		save_solution = true;
		local_time_stamps.push_back (time_layers[current_time_layer]);
	}
	else
	{
		// if there is some scheme
		// save which times fall between last one and current
		for (size_t i = 0, i_end = time_stamps.size (); i < i_end; i++)
		{
			if (time_layers[current_time_layer - 1] + 1e-7 < time_stamps[i] && time_stamps[i] < time_layers[current_time_layer] + 1e-7)
			{
				local_time_stamps.push_back (time_stamps[i]);
			}
		}
		if (local_time_stamps.size () != 0)
		{
			save_solution = true;
		}
	}
	if (save_solution)
	{
		FILE * file;
		MathVector * c = new MathVector (time_sampling);

		char name[64];
		for (int k_system = 0; k_system < n_systems; k_system++)
		{
			// fprint solution
			for (size_t i = 0, i_end = local_time_stamps.size (); i < i_end; i++)
			{
				get_time_approx (local_time_stamps[i], c);
				sprintf (name, "Result//NLNSSTT//s%i_t_%.5lf_nliter_%i.txt", k_system, local_time_stamps[i], knliter);
				file = fopen (name, "w");
				non_linear_layers_solutions[k_system].FPrint (file);
				fclose (file);
			}
		}
		delete c;
	}

}

void Nonlinear_nonstationary_system_task::print_solutions ()
{
	// picture for every non-linear solution 
	std::vector<double> local_time_stamps;
	bool save_solution = false;

	if (use_time_mapping)
	{
		// save solution in current time
		save_solution = true;
		local_time_stamps.push_back (time_layers[current_time_layer]);
	}
	else
	{
		// if there is some scheme
		// save which times fall between last one and current
		for (size_t i = 0, i_end = time_stamps.size (); i < i_end; i++)
		{
			if (time_layers[current_time_layer - 1] + 1e-7 < time_stamps[i] && time_stamps[i] < time_layers[current_time_layer] + 1e-7)
			{
				local_time_stamps.push_back (time_stamps[i]);
			}
		}
		if (local_time_stamps.size () != 0)
		{
			save_solution = true;
		}
	}
	if (save_solution)
	{
		FILE * file;
		MathVector * c = new MathVector (time_sampling);

		char name[64];
		for (int k_system = 0; k_system < n_systems; k_system++)
		{
			// fprint solution
			for (size_t i = 0, i_end = local_time_stamps.size (); i < i_end; i++)
			{
				get_time_approx (local_time_stamps[i], c);
				printf ("\tsaving solutions into files, time layer:\t%.13lf\n", local_time_stamps[i]);
				sprintf (name, "Result//NLNSSTT//s%i_t_%.5lf.txt", k_system, local_time_stamps[i]);
				file = fopen (name, "w");
				double val;
				for (int j = 0, j_end = previous_time_layers_solutions[k_system][0].getSize (); j < j_end; j++)
				{
					val = 0;
					for (int t = 0; t < time_sampling; t++)
					{
						val += c->getElem (t) * previous_time_layers_solutions[k_system][t].getElem (j);
					}
					fprintf (file, "%.13lf\n", val);
				}
				fclose (file);
			}

			char pic_name[128];
			wchar_t wtemp[128];

			// draw solution
			for (size_t i = 0, i_end = local_time_stamps.size (); i < i_end; i++)
			{
				painter->set_max_resolution (1000);
				painter->set_axis_names ("X", "Y");

				painter->set_iso_line_width (1.0);
				painter->set_axes_scale (10, 10);

				painter->draw_field (k_system, COLOR_SCALE_RED_YELLOW);
				painter->draw_contour_lines (k_system, 15);

				sprintf (name, "Pictures//NLNSSTT//s%i_t_%.5lf.png", k_system, local_time_stamps[i]);

				mbstowcs (wtemp, name, strlen (name) + 1);
				std::wstring w_name = wtemp;

				painter->draw_to_file (w_name);
				painter->reset ();
			}

		}
		delete c;

		// TEST slices
		for (int k_system = 0; k_system < n_systems; k_system++)
		{
			for (size_t i = 0, i_end = local_time_stamps.size (); i < i_end; i++)
			{
				{
					char name[64];
					sprintf (name, "Result//NLNSSTT//svert%i_t_%.5lf.txt", k_system, local_time_stamps[i]);
					save_slice (k_system, "Source Files//NLNSSTT//slice_vert.txt", name, false);
				}
				{
					char name[64];
					sprintf (name, "Result//NLNSSTT//shor%i_t_%.5lf.txt", k_system, local_time_stamps[i]);
					save_slice (k_system, "Source Files//NLNSSTT//slice_hor.txt", name, false);
				}
			}
		}
	}
}

void Nonlinear_nonstationary_system_task::set_n_systems ()
{
	n_systems = 3;
	eq_matrixes = new compressed_matrix[n_systems];
	printf ("systems settled:\t%i\n", n_systems);
	fprintf (log_file, "systems settled:\t%i\n", n_systems);
}
