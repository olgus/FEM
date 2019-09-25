#include "Nonlinear_nonstationary_test_task.h"

void Nonlinear_nonstationary_task::build_BM (int k_system)
{
	// copy portrait
	BM[k_system] = eq_matrixes[k_system];
	BM[k_system].Clear ();

	Matrix * G;

	int n_functions;
	int area;
	double lambda;
	int iF, jF;
	int dim = mesh_pointer->get_dimentionality ();

	n_functions = mesh_pointer->get_amount_non_zero_functions (0);
	G = new Matrix (n_functions, n_functions);

	for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
	{
		// get G matrix
		mesh_pointer->get_G_local_matrix (k_element, G);

		// get average lambda on the element
		area = mesh.elements[k_element]->get_area ();
		lambda = 0.0;
		for (int k_function = 0; k_function < n_functions; k_function++)
		{
			lambda += Lambda (k_system, k_function, area);
		}
		lambda /= (double)n_functions;

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
						BM[k_system].add_to_entry (iF, jF, lambda * G->Elem (i, j));
					}
				}
			}
		}
	}

	delete G;
}

double Nonlinear_nonstationary_task::get_function_value (int k_system, int k_function, double * coordinates, int area)
{
	double r = 0.0;
	switch (k_function)
	{
	case 0:
		// sigma
		r = function_sigma (k_system, coordinates, area);
		break;
	case 1:
		r = function_f (k_system, coordinates, area);
		break;
	}
	return r;
}

double Nonlinear_nonstationary_task::function_f (int k_system, double * coordinates, int area)
{
	double x = coordinates[0];
	double y = coordinates[1];
	double t = time_layers[current_time_layer];
	double u = 0.0;
	get_non_linear_solution_in_point (k_system, coordinates, &u);
	//return u;
	//return (x*x + y*y) * u - 4.0 * t;
	//return -4.0 + exp (u) * x * y;
	return 50.0 * u * x * y - 4.0;
}

double Nonlinear_nonstationary_task::function_starting_condition (int k_system, double * coordinates, int area)
{
	return function_FCondition (k_system, coordinates, area, 0);
}

double Nonlinear_nonstationary_task::function_sigma (int k_system, double * coordinates, int area)
{
	double u = 0.0;
	get_non_linear_solution_in_point (k_system, coordinates, &u);
	//return u;
	//return exp (u);
	return 50 * u;
}

double Nonlinear_nonstationary_task::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	double x = coordinates[0];
	double y = coordinates[1];
	double t = time_layers[current_time_layer];
	//return x + y + t;
	//return (x*x + y*y) * t;
	return 1.0 + x * x + y * y + t * x * y;
}

void Nonlinear_nonstationary_task::get_function_coefficients (int k_system, int function, MathVector * coef)
{
	// let's do full one to see if it saves some time because holy shit it sucks all the time


	//double coord[2];
	//int area;
	//for (int i = 0; i < mesh_pointer->get_n_nodes (); i++)
	//{
	//	mesh_pointer->get_node_coordinates (i, coord);
	//	area = mesh_pointer->get_area (coord);
	//	coef->setElem (i, get_function_value (k_system, function, coord, area));
	//}

	{
		Matrix * M;
		MathVector * B;
		int n_functions;
		int n_functions_cur;
		int iF, jF;
		int dim = mesh_pointer->get_dimentionality ();

		n_functions = mesh_pointer->get_amount_non_zero_functions (0);
		M = new Matrix (n_functions, n_functions);
		B = new MathVector (n_functions);

		// integration details
		double jac = 0.0; // jacobian
		int n_integr_points = mesh.elements[0]->amount_of_integration_points ();
		double * weigths = new double[n_integr_points];
		double ** points = new double *[n_integr_points];
		for (int i = 0; i < n_integr_points; i++)
		{
			points[i] = new double[dim];
		}
		double func, val;
		int area;

		// uses eq_matrixes[0] memory !
		eq_matrixes[0].Clear ();
		for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
		{
			n_functions_cur = mesh_pointer->get_amount_non_zero_functions (k_element);
			if (n_functions != n_functions_cur)
			{
				delete M;
				delete B;

				n_functions = n_functions_cur;

				M = new Matrix (n_functions, n_functions);
				B = new MathVector (n_functions);
			}
			// get integration points
			mesh.elements[k_element]->integration_points (points, weigths, &jac);
			area = mesh_pointer->get_area (k_element);

			// get M matrix 
			mesh_pointer->get_M_local_matrix (k_element, M);
			for (int i = 0; i < n_functions; i++)
			{
				val = 0.0;
				// integrate them * basis_functions on the element
				for (int j = 0; j < n_integr_points; j++)
				{
					// sum them multiplying by weigths
					func = get_function_value (k_system, function, points[j], area);
					val += mesh.elements[k_element]->get_basis_function_value (i, points[j]) * func * weigths[j];
				}
				// multiply by jac
				val *= jac;
				B->setElem (i, val);
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

						if (iF != -1)
						{
							// add M into respective places (by def_nodes)
							eq_matrixes[0].add_to_entry (iF, jF, M->Elem (i, j));
						}
					}

					// put B into respective places (by def_nodes)
					eq_matrixes[0].add_to_f_entry (iF, B->getElem (i));
				}
			}
		}

		eq_matrixes[0].solve_LOS (SOLVER_DECOMP_TYPE_LU, 1);
		eq_matrixes[0].get_solution (coef);
	}
}

//void Nonlinear_nonstationary_task::build_system (int k_system)
//{
//	Matrix * G;
//	MathVector * S;
//	MathVector * U_prev;
//	MathVector * F;
//	MathVector * Time_coef;
//	MathVector * DDD;
//
//	int n_functions;
//	int area;
//	double lambda;
//	int iF, jF;
//	int dim = mesh_pointer->get_dimentionality ();
//
//	n_functions = mesh_pointer->get_amount_non_zero_functions (0);
//	G = new Matrix (n_functions, n_functions);
//	S = new MathVector (n_functions);
//	F = new MathVector (n_functions);
//	U_prev = new MathVector (n_functions);
//	DDD = new MathVector (n_functions * n_functions * n_functions);
//
//	int sigma_integrate_coef[] = { 3, 0 };
//	int f_integrate_coef[] = { 1 };
//
//	Time_coef = new MathVector (time_sampling);
//	get_time_coefficients (Time_coef);
//
//	eq_matrixes[k_system].Clear ();
//
//	for (int i = 0; i < n_functions; i++)
//	{
//		S->setElem (i, 1.0);
//		F->setElem (i, 1.0);
//	}
//
//	for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
//	{
//		// get G matrix
//		mesh_pointer->get_G_local_matrix (k_element, G);
//		// get DDD matrix 
//		mesh_pointer->get_DDD_local_matrix (k_element, DDD);
//
//		// get sigma cooeficients
//		get_function_coefficients (k_element, 0, S);
//
//		// get F cooeficients
//		get_function_coefficients (k_element, 1, F);
//
//		// get average lambda on the element
//		area = mesh.elements[k_element]->get_area ();
//		lambda = 0.0;
//		for (int k_function = 0; k_function < n_functions; k_function++)
//		{
//			lambda += Lambda (k_system, k_function, area);
//		}
//		lambda /= (double)n_functions;
//
//		// go by element's functions
//		for (int i = 0; i < n_functions; i++)
//		{
//			iF = get_function_global_number (k_element, i);
//			if (iF != -1)
//			{
//				for (int j = 0; j < n_functions; j++)
//				{
//					jF = get_function_global_number (k_element, j);
//
//					if (jF != -1)
//					{
//						// add G into respective places (by def_nodes)
//						eq_matrixes[k_system].add_to_entry (iF, jF, lambda * G->Elem (i, j));
//
//						// add bfbfbf for sigma into respective places
//						for (int k = 0; k < n_functions; k++)
//						{
//							sigma_integrate_coef[1] = k;
//							//eq_matrixes[k_system].add_to_entry (iF, jF, S->getElem (k) * Time_coef->getElem (0) * mesh.elements[k_element]->integrate (i, j, sigma_integrate_coef));
//							eq_matrixes[k_system].add_to_entry (iF, jF, S->getElem (k) * Time_coef->getElem (0) * DDD->getElem (i * n_functions * n_functions + j * n_functions + k));
//						}
//					}
//				}
//
//				// add f
//				for (int k = 0; k < n_functions; k++)
//				{
//					eq_matrixes[k_system].add_to_f_entry (iF, F->getElem (k) * mesh.elements[k_element]->integrate (i, k, f_integrate_coef));
//				}
//			}
//		}
//
//		// add prev time layers data
//		{
//			// go by amount of time layers that count
//			for (int t = 1; t < time_sampling; t++)
//			{
//				// get previous time layers solutions
//				for (int i = 0; i < n_functions; i++)
//				{
//					iF = get_function_global_number (k_element, i);
//					if (iF != -1)
//					{
//						U_prev->setElem (i, previous_time_layers_solutions[k_system][t - 1].getElem (iF));
//					}
//				}
//				// go by element's functions
//				for (int i = 0; i < n_functions; i++)
//				{
//					iF = get_function_global_number (k_element, i);
//					if (iF != -1)
//					{
//						for (int j = 0; j < n_functions; j++)
//						{
//							for (int k = 0; k < n_functions; k++)
//							{
//								sigma_integrate_coef[1] = k;
//								//eq_matrixes[k_system].add_to_f_entry (iF, S->getElem (k) * U_prev->getElem (j) * Time_coef->getElem (t) * mesh.elements[k_element]->integrate (i, j, sigma_integrate_coef));
//								eq_matrixes[k_system].add_to_f_entry (iF, S->getElem (k) * U_prev->getElem (j) * Time_coef->getElem (t) * DDD->getElem (i * n_functions * n_functions + j * n_functions + k));
//								//printf ("%.13lf\t%.13lf\n", DDD->getElem (i * n_functions * n_functions + j * n_functions + k), mesh.elements[k_element]->integrate (i, j, sigma_integrate_coef));
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//
//	delete G;
//	delete S;
//	delete F;
//	delete DDD;
//	delete U_prev;
//	delete Time_coef;
//}

//void Nonlinear_nonstationary_task::build_system (int k_system)
//{
//	Matrix * G;
//	MathVector * S;
//	MathVector * U_prev;
//	MathVector * F;
//	MathVector * Time_coef;
//	MathVector * DDD;
//
//	int n_functions;
//	int area;
//	double lambda;
//	int iF, jF;
//	int dim = mesh_pointer->get_dimentionality ();
//
//	n_functions = mesh_pointer->get_amount_non_zero_functions (0);
//	G = new Matrix (n_functions, n_functions);
//	S = new MathVector (eq_matrixes[k_system].Size());
//	F = new MathVector (eq_matrixes[k_system].Size ());
//	U_prev = new MathVector (n_functions);
//	DDD = new MathVector (n_functions * n_functions * n_functions);
//
//	int sigma_integrate_coef[] = { 3, 0 };
//	int f_integrate_coef[] = { 1 };
//
//	Time_coef = new MathVector (time_sampling);
//	get_time_coefficients (Time_coef);
//
//	// get sigma cooeficients
//	get_function_coefficients (k_system, 0, S);
//
//	// get F cooeficients
//	get_function_coefficients (k_system, 1, F);
//
//	eq_matrixes[k_system].Clear ();
//	for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
//	{
//		// get G matrix
//		mesh_pointer->get_G_local_matrix (k_element, G);
//		// get DDD matrix 
//		mesh_pointer->get_DDD_local_matrix (k_element, DDD);
//
//		// get average lambda on the element
//		area = mesh.elements[k_element]->get_area ();
//		lambda = 0.0;
//		for (int k_function = 0; k_function < n_functions; k_function++)
//		{
//			lambda += Lambda (k_system, k_function, area);
//		}
//		lambda /= (double)n_functions;
//
//		// go by element's functions
//		for (int i = 0; i < n_functions; i++)
//		{
//			iF = get_function_global_number (k_element, i);
//			if (iF != -1)
//			{
//				for (int j = 0; j < n_functions; j++)
//				{
//					jF = get_function_global_number (k_element, j);
//
//					if (jF != -1)
//					{
//						// add G into respective places (by def_nodes)
//						eq_matrixes[k_system].add_to_entry (iF, jF, lambda * G->Elem (i, j));
//
//						// add bfbfbf for sigma into respective places
//						for (int k = 0; k < n_functions; k++)
//						{
//							sigma_integrate_coef[1] = k;
//							//eq_matrixes[k_system].add_to_entry (iF, jF, S->getElem (k) * Time_coef->getElem (0) * mesh.elements[k_element]->integrate (i, j, sigma_integrate_coef));
//							eq_matrixes[k_system].add_to_entry (iF, jF, S->getElem (get_function_global_number (k_element, k)) * Time_coef->getElem (0) * DDD->getElem (i * n_functions * n_functions + j * n_functions + k));
//						}
//					}
//				}
//
//				// add f
//				for (int k = 0; k < n_functions; k++)
//				{
//					eq_matrixes[k_system].add_to_f_entry (iF, F->getElem (get_function_global_number (k_element, k)) * mesh.elements[k_element]->integrate (i, k, f_integrate_coef));
//				}
//			}
//		}
//
//		// add prev time layers data
//		{
//			// go by amount of time layers that count
//			for (int t = 1; t < time_sampling; t++)
//			{
//				// get previous time layers solutions
//				for (int i = 0; i < n_functions; i++)
//				{
//					iF = get_function_global_number (k_element, i);
//					if (iF != -1)
//					{
//						U_prev->setElem (i, previous_time_layers_solutions[k_system][t - 1].getElem (iF));
//					}
//				}
//				// go by element's functions
//				for (int i = 0; i < n_functions; i++)
//				{
//					iF = get_function_global_number (k_element, i);
//					if (iF != -1)
//					{
//						for (int j = 0; j < n_functions; j++)
//						{
//							for (int k = 0; k < n_functions; k++)
//							{
//								sigma_integrate_coef[1] = k;
//								//eq_matrixes[k_system].add_to_f_entry (iF, S->getElem (k) * U_prev->getElem (j) * Time_coef->getElem (t) * mesh.elements[k_element]->integrate (i, j, sigma_integrate_coef));
//								eq_matrixes[k_system].add_to_f_entry (iF, S->getElem (get_function_global_number (k_element, k)) * U_prev->getElem (j) * Time_coef->getElem (t) * DDD->getElem (i * n_functions * n_functions + j * n_functions + k));
//								//printf ("%.13lf\t%.13lf\n", DDD->getElem (i * n_functions * n_functions + j * n_functions + k), mesh.elements[k_element]->integrate (i, j, sigma_integrate_coef));
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//
//	delete G;
//	delete S;
//	delete F;
//	delete DDD;
//	delete U_prev;
//	delete Time_coef;
//}

void Nonlinear_nonstationary_task::build_system (int k_system)
{
	Matrix * M;
	MathVector * S;
	MathVector * U_prev;
	MathVector * F;
	MathVector * Time_coef;
	MathVector * DDD;

	int n_functions;
	int area;
	double lambda;
	int iF, jF;
	int dim = mesh_pointer->get_dimentionality ();

	n_functions = mesh_pointer->get_amount_non_zero_functions (0);
	M = new Matrix (n_functions);
	S = new MathVector (eq_matrixes[k_system].Size ());
	F = new MathVector (eq_matrixes[k_system].Size ());
	U_prev = new MathVector (n_functions);
	DDD = new MathVector (n_functions * n_functions * n_functions);

	int f_integrate_coef[] = { 1 };

	Time_coef = new MathVector (time_sampling);
	get_time_coefficients (Time_coef);

	// get sigma cooeficients
	get_function_coefficients (k_system, 0, S);

	// get F cooeficients
	get_function_coefficients (k_system, 1, F);

	// copy BM
	eq_matrixes[k_system] = BM[k_system];

	for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
	{
		// get DDD matrix 
		mesh_pointer->get_DDD_local_matrix (k_element, DDD);
		// get M matrix
		mesh_pointer->get_M_local_matrix (k_element, M);
		
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
						// add bfbfbf for sigma into respective places
						for (int k = 0; k < n_functions; k++)
						{
							eq_matrixes[k_system].add_to_entry (iF, jF, S->getElem (get_function_global_number (k_element, k)) * Time_coef->getElem (0) * DDD->getElem (i * n_functions * n_functions + j * n_functions + k));
						}
					}
				}

				// add f
				for (int k = 0; k < n_functions; k++)
				{
					//eq_matrixes[k_system].add_to_f_entry (iF, F->getElem (get_function_global_number (k_element, k)) * mesh.elements[k_element]->integrate (i, k, f_integrate_coef));

					eq_matrixes[k_system].add_to_f_entry (iF, F->getElem (get_function_global_number (k_element, k)) * M->Elem (i, k));
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
						U_prev->setElem (i, previous_time_layers_solutions[k_system][t - 1].getElem (iF));
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
							for (int k = 0; k < n_functions; k++)
							{
								eq_matrixes[k_system].add_to_f_entry (iF, S->getElem (get_function_global_number (k_element, k)) * U_prev->getElem (j) * Time_coef->getElem (t) * DDD->getElem (i * n_functions * n_functions + j * n_functions + k));
							}
						}
					}
				}
			}
		}
	}

	delete M;
	delete S;
	delete F;
	delete DDD;
	delete U_prev;
	delete Time_coef;
}

void Nonlinear_nonstationary_task::print_solutions (int knliter)
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
				sprintf (name, "Result//NLNSTT//s%i_t_%.5lf_nliter_%i.txt", k_system, local_time_stamps[i], knliter);
				file = fopen (name, "w");
				non_linear_layers_solutions[k_system].FPrint (file);
				fclose (file);
			}
		}
		delete c;
	}

}

void Nonlinear_nonstationary_task::print_solutions ()
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
				sprintf (name, "Result//NLNSTT//s%i_t_%.5lf.txt", k_system, local_time_stamps[i]);
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

				painter->draw_field (0, COLOR_SCALE_RED_YELLOW);
				painter->draw_contour_lines (k_system, 15);

				sprintf (name, "Pictures//NLNSTT//s%i_t_%.5lf.png", k_system, local_time_stamps[i]);

				painter->set_picture_name (pic_name);
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
					sprintf (name, "Result//NLNSTT//svert%i_t_%.5lf.txt", k_system, local_time_stamps[i]);
					save_slice (0, "Source Files//NLNSTT//slice_vert.txt", name, true);
				}
				{
					char name[64];
					sprintf (name, "Result//NLNSTT//shor%i_t_%.5lf.txt", k_system, local_time_stamps[i]);
					save_slice (0, "Source Files//NLNSTT//slice_hor.txt", name, true);
				}
			}
		}
	}
}

bool Nonlinear_nonstationary_task::solve_task (int solver_param[][5])
{

	{
		if (!prepared)
		{
			printf ("ERROR: missing preparation\n");
			return false;
		}
		// set base matrices
		BM = new compressed_matrix [n_systems];
		for (int k_system = 0; k_system < n_systems; k_system++)
		{
			build_BM (k_system);
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
					//relax (k_system, 1.0, prev_non_linear_layers_solutions[k_system]);
					w_rel = relax (k_system, prev_non_linear_layers_solutions[k_system]);
					//w_rel = relax (k_system, 0.0, 1.0, prev_non_linear_layers_solutions[k_system]);

					// non_linear_layers_solutions contains entire _new_ solution
					prev_non_linear_layers_solutions->Substract (non_linear_layers_solutions[k_system]);
					non_linear_dif[k_system] = prev_non_linear_layers_solutions->Norm () / non_linear_layers_solutions[k_system].Norm ();

					printf ("\t\tk_system: %i\tdif: %.3e\tw_rel: %.5lf\n", k_system, non_linear_dif[k_system], w_rel);
					sum_non_linear_dif += non_linear_dif[k_system];

					print_solutions (k_nl_iter);
				}
				if ((k_nl_iter > 0) && (sum_non_linear_dif < NLNSTT_NL_PRECIS))
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

	delete[] BM;
}

double Nonlinear_nonstationary_task::relax (int k_system, double w, const MathVector & prev_non_linear_layer_solution)
{
	eq_matrixes[k_system].get_solution (&(non_linear_layers_solutions[k_system])); // u k
	non_linear_layers_solutions[k_system].Linear_Combination (w, prev_non_linear_layer_solution, 1.0 - w);
	return 0.0;
}

double Nonlinear_nonstationary_task::relax (int k_system, double w_left, double w_right, const MathVector & prev_non_linear_layer_solution)
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

double Nonlinear_nonstationary_task::relax (int k_system, const MathVector & prev_non_linear_layer_solution)
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
		res_omp = new double[(par_num_threads - 1) * eq_matrixes[k_system].Size()];

	// left
	non_linear_layers_solutions[k_system].Linear_Combination (*sol, x1, prev_non_linear_layer_solution, 1.0 - x1);
	build_system (k_system);
	apply_boundary_conditions (k_system);
	eq_matrixes[k_system].mult_A_v (non_linear_layers_solutions[k_system], mv, res_omp, par_num_threads);
	mv->Substract (*(eq_matrixes[k_system].f));
	y1 = mv->Norm ();
	printf ("%lf %.3e\n", x1, y1);
	// right
	non_linear_layers_solutions[k_system].Linear_Combination (*sol, x2, prev_non_linear_layer_solution, 1.0 - x2);
	build_system (k_system);
	apply_boundary_conditions (k_system);
	eq_matrixes[k_system].mult_A_v (non_linear_layers_solutions[k_system], mv, res_omp, par_num_threads);
	mv->Substract (*(eq_matrixes[k_system].f));
	y2 = mv->Norm ();
	printf ("%lf %.3e\n", x2, y2);

	for (int i = 0; (i < 100) && (fabs (b - a) > 2e-2); i++)
	{
		printf ("%i", i);
		if (y1 >= y2)
		{
			a = x1;
			x1 = x2;
			y1 = y2;

			x2 = a + (b - a) / c;
			non_linear_layers_solutions[k_system].Linear_Combination (*sol, x2, prev_non_linear_layer_solution, 1.0 - x2);
						
			t1 = clock ();

			build_system (k_system);

			t2 = clock ();
			ms = (double)(t2 - t1);
			printf ("time spend building system: %lf s\n", ms / 1000);

			t1 = clock ();
			
			apply_boundary_conditions (k_system);

			t2 = clock ();
			ms = (double)(t2 - t1);
			printf ("time spend applying boundary cond: %lf s\n", ms / 1000);

			t1 = clock ();
			
			eq_matrixes[k_system].mult_A_v (non_linear_layers_solutions[k_system], mv, res_omp, par_num_threads);

			t2 = clock ();
			ms = (double)(t2 - t1);
			printf ("time spend multiplying: %lf s\n", ms / 1000);

			mv->Substract (*(eq_matrixes[k_system].f));
			y2 = mv->Norm ();
			printf ("%lf %.3e\n", x2, y2);
		}
		else
		{
			b = x2;
			x2 = x1;
			y2 = y1;

			x1 = b - (b - a) / c;
			non_linear_layers_solutions[k_system].Linear_Combination (*sol, x1, prev_non_linear_layer_solution, 1.0 - x1);
			t1 = clock ();

			build_system (k_system);

			t2 = clock ();
			ms = (double)(t2 - t1);
			printf ("time spend building system: %lf s\n", ms / 1000);

			t1 = clock ();

			apply_boundary_conditions (k_system);

			t2 = clock ();
			ms = (double)(t2 - t1);
			printf ("time spend applying boundary cond: %lf s\n", ms / 1000);

			t1 = clock ();

			eq_matrixes[k_system].mult_A_v (non_linear_layers_solutions[k_system], mv, res_omp, par_num_threads);

			t2 = clock ();
			ms = (double)(t2 - t1);
			printf ("time spend multiplying: %lf s\n", ms / 1000);

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

void Nonlinear_nonstationary_task::sigma_test ()
{
	// comment deletion of non_linear_layers_solutions in solve_task

	FILE * file = fopen ("Test//sigma.txt", "w");
	MathVector * Sigma = new MathVector (eq_matrixes[0].Size());
	Sigma->Zero ();

	int n_functions = mesh_pointer->get_amount_non_zero_functions (0);
	MathVector * S = new MathVector (n_functions);

	for (int i = 0; i < mesh_pointer->get_n_elements (); i++)
	{
		get_function_coefficients (i, 0, S);
		for (int k = 0; k < n_functions; k++)
		{
			Sigma->setElem (get_function_global_number (i, k), S->getElem (k));
		}
	}

	Sigma->FPrint (file);

	fclose (file);
	delete Sigma;
	delete S;
}
