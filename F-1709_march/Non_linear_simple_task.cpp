#include "Non_linear_simple_task.h"

double Non_linear_simple_task::get_function_value (int k_system, int k_function, double * coordinates, int area)
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

double Non_linear_simple_task::function_f (int k_system, double * coordinates, int area)
{
	double x = coordinates[0];
	double y = coordinates[1];
	double u = 0.0;
	get_non_linear_solution_in_point (k_system, coordinates, &u);
	//return (x + y) * u;
	//return (x*x + y*y) * u - 4.0;
	return -4.0;
}

double Non_linear_simple_task::function_sigma (int k_system, double * coordinates, int area)
{
	double u = 0.0;
	get_non_linear_solution_in_point (k_system, coordinates, &u);
	//return u;
	return 0.0;
}

double Non_linear_simple_task::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	double x = coordinates[0];
	double y = coordinates[1];
	//return x + y;
	return x*x + y*y;
}

void Non_linear_simple_task::build_system (int k_system)
{
	Matrix * G;
	MathVector * S;
	MathVector * F;
	int n_functions;
	int area;
	double lambda;
	int iF, jF;
	int dim = mesh_pointer->get_dimentionality ();

	n_functions = mesh_pointer->get_amount_non_zero_functions (0);
	G = new Matrix (n_functions, n_functions);
	S = new MathVector (n_functions);
	F = new MathVector (n_functions);

	int sigma_integrate_coef[] = {3, 0};
	int f_integrate_coef[] = {1};
	for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
	{
		// get G matrix
		mesh_pointer->get_G_local_matrix (k_element, G);

		// get sigma cooeficients
		get_function_coefficients (k_system, k_element, 0, S);

		// get F cooeficients
		get_function_coefficients (k_system, k_element, 1, F);

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
						eq_matrixes[k_system].add_to_entry (iF, jF, lambda * G->Elem (i, j));

						// add M for sigma into respective places
						for (int k = 0; k < n_functions; k++)
						{
							sigma_integrate_coef[1] = k;
							eq_matrixes[k_system].add_to_entry (iF, jF, S->getElem (k) * mesh.elements[k_element]->integrate (i, j, sigma_integrate_coef));
						}
					}
				}

				// add f
				for (int k = 0; k < n_functions; k++)
				{
					eq_matrixes[k_system].add_to_f_entry (iF, F->getElem (k) * mesh.elements[k_element]->integrate (i, k, f_integrate_coef));
				}
			}
		}
	}

	delete G;
	delete S;
	delete F;
}

void Non_linear_simple_task::print_solutions (int knliter)
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
				sprintf (name, "Result//NLST//s%i_nliter_%i_t_%.5lf.txt", k_system, knliter, local_time_stamps[i]);
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

				sprintf (name, "Pictures//NLST//s%i_nliter_%i_t_%.5lf.png", k_system, knliter, local_time_stamps[i]);

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
					sprintf (name, "Result//NLST//svert%i_nliter_%i_t_%.5lf.txt", k_system, knliter, local_time_stamps[i]);
					save_slice (0, "Source Files//NLST//slice_vert.txt", name, true);
				}
				{
					char name[64];
					sprintf (name, "Result//NLST//shor%i_nliter_%i_t_%.5lf.txt", k_system, knliter, local_time_stamps[i]);
					save_slice (0, "Source Files//NLST//slice_hor.txt", name, true);
				}
			}
		}
	
	}
}

bool Non_linear_simple_task::solve_task (int solver_param[][5])
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

	// 0 non_linear solution is a start solution
	for (int k_system = 0; k_system < n_systems; k_system++)
	{
		non_linear_layers_solutions[k_system].Copy ((previous_time_layers_solutions[k_system][0]));
	}

	int solver_iterations;
	bool found = false;
	for (int k_nl_iter = 0; k_nl_iter < MAX_NONLINEAR_ITER && !found; k_nl_iter++)
	{
		printf ("k_nl_iter: %i ", k_nl_iter);
		sum_non_linear_dif = 0.0;

		for (int k_system = 0; k_system < n_systems; k_system++)
		{
			// save prev non linear solution
			prev_non_linear_layers_solutions->Copy (non_linear_layers_solutions[k_system]);
			eq_matrixes[k_system].Clear ();
			// set starting point for the solver
			eq_matrixes[k_system].set_starting_point (&(non_linear_layers_solutions[k_system])); // set last layer solution as x0 for solver
			// build system
			build_system (k_system); 
			// apply boundary conditions
			apply_boundary_conditions (k_system); // apply boundary conditions
			// solve it
			solver_iterations = eq_matrixes[k_system].solve (solver_param[k_system]); 	

			// TEST get_nonlinear_solution
			{
				double coord[2];
				double v;
				for (int k_system = 0; k_system < n_systems; k_system++)
				{
					{
						char name[64];
						sprintf (name, "Result//NLST//s%i_prev_nonl_sol%i.txt", k_system, k_nl_iter);
						FILE * file = fopen (name, "w");
						for (int k_node = 0, k_node_N = mesh_pointer->get_n_nodes (); k_node < k_node_N; k_node++)
						{
							mesh_pointer->get_node_coordinates (k_node, coord);
							get_non_linear_solution_in_point (k_system, coord, &v);
							fprintf (file, "%.13lf\n", v);
						}
						fclose (file);
					}
				}
			}

			// TODO relax
			// for now just copy the solution from the matrix
			// both are needed in this one
			eq_matrixes[k_system].get_solution (&(non_linear_layers_solutions[k_system]));
			eq_matrixes[k_system].get_solution (&(previous_time_layers_solutions[k_system][0]));

			// non_linear_layers_solutions contains entire _new_ solution
			prev_non_linear_layers_solutions->Substract (non_linear_layers_solutions[k_system]);
			non_linear_dif[k_system] = prev_non_linear_layers_solutions->Norm () / non_linear_layers_solutions[k_system].Norm ();
			
			if (k_system != 0) printf ("              ");
			printf ("k_system: %i\tdif: %.3e\n", k_system, non_linear_dif[k_system]);
			sum_non_linear_dif += non_linear_dif[k_system];

			print_solutions (k_nl_iter);
		}
		if (sum_non_linear_dif < NLST_NL_PRECIS)
			found = true;
	}

	delete[] non_linear_dif;
	delete prev_non_linear_layers_solutions;
	delete[] non_linear_layers_solutions;
	non_linear_layers_solutions = NULL;
	return true;
}

double Non_linear_simple_task::relax (int k_system, MathVector * prev_solution)
{
	return 0.0;
}
