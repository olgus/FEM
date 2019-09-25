#include "Stephan_AMT.h"

template class Stephan_AMT_smoothing <Mesh_1D_L1>;
template class Stephan_AMT_coord_transf <Mesh_1D_L1>;

double autom_St_task::Laplace_function (double z)
{
	double h = 0.1;
	double zl, zr;
	bool fin = false;
	double r = 0.0;
	zl = 0.0;
	for (int i = 0; !fin; i++)
	{
		zr = (i + 1) * h;
		// check if hit boundary
		if (zr > z - PRECIS)
		{
			fin = true;
			zr = z;
		}
		// gauss 4
		{
			double weigths[4];
			weigths[0] = (18.0 + pow (30, 0.5)) / 36.0;
			weigths[1] = (18.0 + pow (30, 0.5)) / 36.0;
			weigths[2] = (18.0 - pow (30, 0.5)) / 36.0;
			weigths[3] = (18.0 - pow (30, 0.5)) / 36.0;

			double xi[4];
			xi[0] = -pow (3.0 / 7.0 - 2.0 * pow (6.0 / 5.0, 0.5) / 7.0, 0.5);
			xi[1] = pow (3.0 / 7.0 - 2.0 * pow (6.0 / 5.0, 0.5) / 7.0, 0.5);
			xi[2] = -pow (3.0 / 7.0 + 2.0 * pow (6.0 / 5.0, 0.5) / 7.0, 0.5);
			xi[3] = pow (3.0 / 7.0 + 2.0 * pow (6.0 / 5.0, 0.5) / 7.0, 0.5);

			double jac = (zr - zl) / 2.0;
			double r_segm = 0.0;
			double f, y;
			for (int k = 0; k < 4; k++)
			{
				y = jac * (1.0 + xi[k]) + zl;
				f = exp (-pow (y, 2.0));
				r_segm += weigths[k] * f;
			}
			r_segm *= jac;
			r += r_segm;
		}
		zl = zr;
	}

	r *= 2.0 / sqrt (PI);
	return r;
}

double autom_St_task::get_ph_ch_position (double t, double beta)
{
	return 2.0 * beta * sqrt (t);
}

double autom_St_task::get_solution (double t, double x, double beta)
{
	double st = 2.0 * sqrt (t);
	double pcb = st * beta;
	if ((x > pcb) || (st < PRECIS))
		return 0.0;
	return -1.0 + Laplace_function (x / st) / Laplace_function (beta);
}

template<class TypeMesh>
void Stephan_AMT<TypeMesh>::set_n_systems ()
{
	n_systems = 1;
	eq_matrixes = new compressed_matrix[n_systems];
	printf ("systems settled:\t%i\n", n_systems);
}

template<class TypeMesh>
double Stephan_AMT<TypeMesh>::function_starting_condition (int k_system, double * coordinates, int area)
{
	double t = time_layers[current_time_layer];
	double x = coordinates[0];
	double y = coordinates[1];
	double z = coordinates[0] * cos (angle) + coordinates[1] * sin (angle);
	//return 0.0;
	return autom_St_task::get_solution (t, z, beta);
}

template<class TypeMesh>
double Stephan_AMT<TypeMesh>::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	double r;
	// left vertical boundary cool, right hot
	switch (boundary)
	{
	case 0:
		r = T_cool;
		break;
	case 1:
		r = T_heat;
		break;

	}

	return r;
}

template<class TypeMesh>
double Stephan_AMT<TypeMesh>::function_f (int k_system, double * coordinates, int area)
{
	// no source
	return 0.0;
}

template<class TypeMesh>
void Stephan_AMT<TypeMesh>::print_solutions ()
{
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

		for (int k_system = 0; k_system < n_systems; k_system++)
		{
			for (size_t i = 0, i_end = local_time_stamps.size (); i < i_end; i++)
			{
				get_time_approx (local_time_stamps[i], c);

				printf ("\tsaving solutions into files, time layer:\t%.16lf\n", local_time_stamps[i]);
				//sprintf (name, "Result//s%i_t_%0.5lf.txt", k, local_time_stamps[i]);
				//file = fopen (name, "w");
				//previous_time_layers_solutions[k][0].FPrint (file);

				//double val;
				//for (int j = 0, j_end = previous_time_layers_solutions[k][0].getSize (); j < j_end; j++)
				//{
				//	val = 0;
				//	for (int t = 0; t < time_sampling; t++)
				//	{
				//		val += c->getElem (t) * previous_time_layers_solutions[k][t].getElem (j);
				//	}
				//	fprintf (file, "%.16lf\n", val);
				//}
				//fclose (file);

				if (painter != NULL)
				{
					char name[128];
					wchar_t wtemp[128];
					{
						painter->set_max_resolution (1000);
						painter->set_axis_names ("X", "Y");
						painter->draw_contour_lines (k_system, 20);
						painter->add_isoline (phys_T_phase_change);
						//painter->no_legend ();
						painter->draw_field (0, COLOR_SCALE_RED_YELLOW);

						sprintf (name, "Pictures//Stephan//s%i_t_%lf_iso.png", k_system, local_time_stamps[i]);
						mbstowcs (wtemp, name, strlen (name) + 1);
						std::wstring w_name = wtemp;
						painter->draw_to_file (w_name);
						painter->reset ();
					}
				}
			}
		}
		delete c;
	}
}

template<class TypeMesh>
double Stephan_AMT<TypeMesh>::get_ph_ch_x_position ()
{
	return 0.0;
}

template<class TypeMesh>
void Stephan_AMT<TypeMesh>::painter_pointer (Painter * p)
{
	painter = p;
}

template<class TypeMesh>
Stephan_AMT<TypeMesh>::Stephan_AMT ()
{
	non_linear = true;
	angle = 30;
	beta = 0.62;

	phase_change_og_position = 0;
}

template<class TypeMesh>
Stephan_AMT<TypeMesh>::Stephan_AMT (char * parameters_file_name)
{
	painter = NULL;
	non_linear = true;

	FILE * file = fopen (parameters_file_name, "r");
	fscanf (file, "%lf", &angle);
	angle = angle / 180 * PI; // radians

	fscanf (file, "%lf", &beta);
	fscanf (file, "%lf", &del_T_smoothing);
	fscanf (file, "%lf", &phys_T_phase_change);
	fscanf (file, "%lf", &T_heat);
	fscanf (file, "%lf", &T_cool);
	fscanf (file, "%lf", &L);
	fclose (file);
}

template<class TypeMesh>
Stephan_AMT<TypeMesh>::~Stephan_AMT ()
{
}

template<class TypeMesh>
void Stephan_AMT_smoothing<TypeMesh>::build_system (int k_system)
{
	Matrix * G;
	Matrix * M;
	MathVector * B;
	MathVector * Q;
	MathVector * MQ;
	MathVector * Time_coef;


	int n_points = mesh_pointer->amount_of_integration_points (0);
	double * weights = new double[n_points];
	double ** points = new double *[n_points];
	for (int k = 0; k < n_points; k++)
		points[k] = new double[2];

	int n_functions;
	int area;
	double lambda, gamma;
	int iF, jF;
	int dim = mesh_pointer->get_dimentionality ();

	Time_coef = new MathVector (time_sampling);
	get_time_coefficients (Time_coef);

	for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
	{
		n_functions = mesh_pointer->get_amount_non_zero_functions (k_element);
		// get G, M matrices 
		G = new Matrix (n_functions, n_functions);
		M = new Matrix (n_functions, n_functions);
		mesh_pointer->get_G_local_matrix (k_element, G);
		mesh_pointer->get_M_local_matrix (k_element, M);

		// get B-vector
		B = new MathVector (n_functions);
		get_local_B (k_system, k_element, B);

		// get average lambda/gamma on the element 
		area = mesh.elements[k_element]->get_area ();
		lambda = 0.0;
		gamma = 0.0;
		for (int k_function = 0; k_function < n_functions; k_function++)
		{
			lambda += Lambda (k_system, k_function, area);
			gamma += Gamma (k_system, k_function, area);
		}
		lambda /= (double)n_functions;
		gamma /= (double)n_functions;

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
						// add M into respective places (by def_nodes)
						eq_matrixes[k_system].add_to_entry (iF, jF, gamma * M->Elem (i, j) * Time_coef->getElem (0));
					}
				}
				// put B into respective places (by def_nodes)
				eq_matrixes[k_system].add_to_f_entry (iF, B->getElem (i));
			}
		}

		Q = new MathVector (n_functions);
		MQ = new MathVector (n_functions);
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
		if (L > PRECIS)
		{
			// add heat source to the left side of equation
			// L * time_coef * integral
			// integral = Dirac * psi i * psi j
			// go by element's functions
			// go by element's functions
			{
				// get integration points and everything
				double jac;
				mesh_pointer->integration_points (k_element, points, weights, &jac);
				double integral_value;
				double T;
				// make M, integrals L * bf_i * bf_j
				M->Zero ();
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
								integral_value = 0.0;
								for (int k = 0; k < n_points; k++)
								{
									// get solution, if T is non-linear
									{
										// get local functions' values of the element
										int func_amount = mesh_pointer->get_amount_non_zero_functions (k_element);
										MathVector * local_func = new MathVector (func_amount);
										mesh_pointer->get_local_function_values (k_element, points[k], local_func);
										// make vector of solutions in those nodes
										MathVector * solution = new MathVector (func_amount);
										int global_function;
										for (int i = 0; i < func_amount; i++)
										{
											global_function = get_function_global_number (k_element, i);
											if (global_function != -1)
											{
												solution->setElem (i, non_linear_layers_solutions[0].getElem (global_function));
											}
										}
										// multiply them by the solution's values in those nodes
										T = local_func->Scalar_Product (*solution);
										delete local_func;
										delete solution;
									}
									integral_value += weights[k] * mesh_pointer->get_basis_function_value (k_element, i, points[k]) * mesh_pointer->get_basis_function_value (k_element, j, points[k]) * Dirac_function (T, points[k]);
								}
								integral_value *= jac;
								M->setElem (i, j, integral_value);
							}
						}
					}
				}
				// M goes into A
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
								eq_matrixes[k_system].add_to_entry (iF, jF, M->Elem (i, j) * Time_coef->getElem (0) * L);
							}
						}
					}
				}
				// add heat source to the right side of equation
				// go by amount of time layers that count
				for (int t = 1; t < time_sampling; t++)
				{
					// get q's
					for (int i = 0; i < n_functions; i++)
					{
						iF = get_function_global_number (k_element, i);
						if (iF != -1)
						{
							Q->setElem (i, previous_time_layers_solutions[0][t - 1].getElem (iF));
						}
					}
					// multiply it by vector of T from respective layer
					M->MultiplyMatrixByVector (*Q, MQ);
					// add to f
					for (int i = 0; i < n_functions; i++)
					{
						iF = get_function_global_number (k_element, i);
						if (iF != -1)
						{
							eq_matrixes[k_system].add_to_f_entry (iF, MQ->getElem (i) * Time_coef->getElem (t) * L);
						}
					}
				}
			}
		}
		delete B;
		delete G;
		delete M;
		delete Q;
		delete MQ;
	}

	delete[] weights;
	for (int k = 0; k < n_points; k++)
		delete[] points[k];
	delete[] points;
	delete Time_coef;
}

template<class TypeMesh>
double Stephan_AMT_smoothing<TypeMesh>::Dirac_function (double T, double * coordinates)
{
	double r = 0.0;
	double t = time_layers[current_time_layer];
	double x = coordinates[0];
	double y = coordinates[1];
	r = exp (-pow (T - phys_T_phase_change, 2.0) / del_T_smoothing);
	return r;
}


template<class TypeMesh>
Stephan_AMT_smoothing<TypeMesh>::Stephan_AMT_smoothing ()
{
}

template<class TypeMesh>
Stephan_AMT_smoothing<TypeMesh>::Stephan_AMT_smoothing (char * parameters_file_name)
{
	Stephan_AMT::Stephan_AMT (parameters_file_name);
}

template<class TypeMesh>
Stephan_AMT_smoothing<TypeMesh>::Stephan_AMT_smoothing (const Stephan_AMT_smoothing & samts)
{
}

template<class TypeMesh>
Stephan_AMT_smoothing<TypeMesh>::~Stephan_AMT_smoothing ()
{
}

//template<class TypeMesh>
void Stephan_AMT_mesh_adapt::adapt_mesh ()
{
	// find specific point with phys_T_phase_change value
	double point[1];
	int k_element = -1;
	{
		double localMin = 1e+20;
		double localMax = -1e+20;
		double c[1];
		double v;
		int n[2];
		for (int k = 0, k_end = mesh_pointer->get_n_nodes (); k < k_end && k_element == -1; k++)
		{
			mesh_pointer->get_base_nodes (k, n);
			for (int i = 0; i < 2; i++)
			{
				mesh_pointer->get_node_coordinates (n[i], c);
				get_solution_in_point (0, c, &v);
				if (v > localMax)
					localMax = v;
				if (v < localMin)
					localMin = v;
			}
			if ((localMin - PRECIS < phys_T_phase_change) && (phys_T_phase_change < localMax + PRECIS))
			{
				k_element = k;
			}
		}
	}

	int new_L_heat_node;
	if (k_element != -1)
	{
		// find the point itself
		// COMMENT:
		// function uses previous_time_layers_solution
		get_isoline_section (0, -1, k_element, phys_T_phase_change, point, point);
		new_L_heat_node = -1;
		int n[2];
		// also check that it is not one of the already existing nodes
		{
			double c[1];
			mesh_pointer->get_base_nodes (k_element, n);
			// if it is, just keep it as L_heat_node
			for (int i = 0; i < 2 && new_L_heat_node == -1; i++)
			{
				mesh_pointer->get_node_coordinates (n[i], c);
				if (fabs (c[0] - point[0]) < PRECIS)
					new_L_heat_node = n[i];
			}
		}

		// remove previous L_heat_node
		if (!L_base_node)
		{
			//mesh.nodes.erase (std::remove (mesh.nodes.begin (), mesh.nodes.end (), mesh.nodes.begin () + L_heat_node));
		}

		L_heat_node = new_L_heat_node;
		if (new_L_heat_node == -1)
		{
			// add the node
			// find the place for the node
			Node_1D node;
			node.set_coordinates (point);
			mesh.nodes.insert (mesh.nodes.begin() + n[0] + 1, node);
			L_base_node = false;
		}
		else
		{
			L_base_node = true;
		}
	}
}

void Stephan_AMT_mesh_adapt::reset_areas ()
{
	int n[2];
	int cur_area = 1;
	for (int i = 0; i < mesh.get_n_elements (); i++)
	{
		// if n[0] < L_heat_node then 1, otherwise 2
		if (cur_area == 1)
		{
			mesh.get_base_nodes (i, n);
			if (n[0] > L_heat_node)
				cur_area = 2;
		}
		mesh.elements[i]->set_area (cur_area);
	}
}

bool Stephan_AMT_mesh_adapt::solve_task (int solver_param[][5])
{
	// TEST
	if (!prepared)
	{
		printf ("ERROR: missing preparation\n");
		return false;
	}

	int nonlinear_iter = 1;
	if (non_linear)
		nonlinear_iter = MAX_NONLINEAR_ITER;
	non_linear_layers_solutions = new MathVector[n_systems];
	for (int k_system = 0; k_system < n_systems; k_system++)
	{
		non_linear_layers_solutions[k_system].setSize (eq_matrixes[k_system].Size ()); // default starting solution point is 0
	}

	mesh_pointers_by_time_levels = new Mesh_Prototype* [n_time_layers];
	mesh_pointers_by_time_levels[0] = mesh_pointer;
	for (int i = 1; i < n_time_layers; i++)
		mesh_pointers_by_time_levels[i] = NULL;
	
	double * discr = new double[n_systems];
	double * prev_discr = new double[n_systems];
	bool stable = false; // indicator of reaching stable condition
	bool found = false;
	int solver_iterations;
	double sum_discr;
	double T_change;
	int unchanged_discr;

	// adapt current mesh to add a L_heat_node
	adapt_mesh ();

	for (current_time_layer = 1; current_time_layer < n_time_layers && !stable; current_time_layer++)
	{
		printf ("current time: %.7lf\n", time_layers[current_time_layer]);
		//printf ("%.7lf ", time_layers[current_time_layer]);
		// save previous time layers solution for non-linearity
		for (int k_system = 0; k_system < n_systems; k_system++)
		{
			non_linear_layers_solutions[k_system].Copy ((previous_time_layers_solutions[k_system][0]));
		}
		found = false;
		for (int i = 0; i < nonlinear_iter && !found; i++)
		{
			printf ("\tnonlinear iteration: %i\n", i);
			// iterate through each system of equations
			for (int k_system = 0; k_system < n_systems; k_system++)
			{
				eq_matrixes[k_system].Clear ();
				//eq_matrixes[k_system].set_starting_point (&(non_linear_layers_solutions[k_system])); // set last layer solution as x0 for solver
				eq_matrixes[k_system].set_starting_point (&(previous_time_layers_solutions[k_system][0])); // set last layer solution as x0 for solver

				build_system (k_system);  // reset system for current layer
				apply_boundary_conditions (k_system); // apply boundary conditions
				solver_iterations = eq_matrixes[k_system].solve (solver_param[k_system]); // solve it	
				discr[k_system] = relax (k_system, &(non_linear_layers_solutions[k_system]));
			}
			if (i > 0)
			{
				sum_discr = 0.0;
				unchanged_discr = 0;
				for (int k_system = 0; k_system < n_systems; k_system++)
				{
					sum_discr += discr[k_system];
					if (fabs (discr[k_system] - prev_discr[k_system]) < DISCR_CHANGE)
						unchanged_discr++;
				}
				if (sum_discr / (double)n_systems < NONLINEAR_DISCR)
					found = true;
				if (unchanged_discr >= n_systems)
					found = true;
				T_change = non_linear_layers_solutions[0].sum_sq_dif (previous_time_layers_solutions[0][0]);
				printf ("%i: %e %e %e\tT: %e %e\n", unchanged_discr, discr[0] - prev_discr[0], discr[1] - prev_discr[1], discr[2] - prev_discr[2], T_change, T_change / non_linear_layers_solutions[0].Norm ());
			}
			for (int k_system = 0; k_system < n_systems; k_system++)
			{
				prev_discr[k_system] = discr[k_system];
			}
		}

		// move meshes
		{
			// save old mesh by making new object for it
			Mesh_1D_L1 old_mesh;
			old_mesh.copy (mesh);
			mesh_pointers_by_time_levels[0] = &old_mesh;
			// move meshes up
			for (int i = n_time_layers - 1; i <= 0; i--)
			{
				mesh_pointers_by_time_levels[i] = mesh_pointers_by_time_levels[i - 1];
			}
			// adapt by adding the L_heat_node
			adapt_mesh ();
			// build elements on those nodes
			mesh.reset_elements ();
			reset_areas ();
			mesh_pointers_by_time_levels[0] = &mesh;

		}
		// move layers
		for (int k_system = 0; k_system < n_systems; k_system++)
			next_time_layer (k_system, &(non_linear_layers_solutions[k_system]));

		print_solutions ();
		print_extra_data ();
	}
	current_time_layer--;

	delete[] discr;
	delete[] prev_discr;
	delete[] non_linear_layers_solutions;
	non_linear_layers_solutions = NULL;
	return true;
}

void Stephan_AMT_mesh_adapt::set_starting_conditions ()
{
}

Stephan_AMT_mesh_adapt::Stephan_AMT_mesh_adapt ()
{
	mesh_pointers_by_time_levels = NULL;
}

Stephan_AMT_mesh_adapt::~Stephan_AMT_mesh_adapt ()
{
	if (mesh_pointers_by_time_levels != NULL)
	{
		// delete all meshes except 0, it gets deleted in parent destructor
		for (int i = 1; i < time_sampling; i++)
		{
			delete mesh_pointers_by_time_levels[i];
		}
	}
}
