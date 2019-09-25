#include "Melt_SF.h"

template class Task_Melt_SF <Triangular_Mesh>;
template class Task_Melt_SF <Triangular_Mesh_Hier>;
template class Task_Melt_SF <Mixed_Triangular_Mesh>;

template <class Type_Mesh>
void Task_Melt_SF<Type_Mesh>::set_n_systems ()
{
	n_systems = 3;
	eq_matrixes = new compressed_matrix[n_systems];
	printf ("systems settled:\t%i\n", n_systems);
	fprintf (log_file, "systems settled:\t%i\n", n_systems);
}

template<class TypeMesh>
double Task_Melt_SF<TypeMesh>::function_starting_condition (int k_system, double * coordinates, int area)
{
	double r = 0.0;
	switch (k_system)
	{
	case 0:
	{
		switch (area)
		{
		case 1:
			r = phys_T0_F;
			break;
		case 2:
			r = phys_T0_S;
			break;
		}
		break;
	}
	case 1:
	case 2:
	{
		r = 0.0;
		break;
	}
	}
	return r;
}

template<class TypeMesh>
void Task_Melt_SF<TypeMesh>::set_starting_conditions ()
{
	if (flag_calc_start)
	{
		for (int k_system = 0; k_system < n_systems; k_system++)
		{
			for (int i = 0; i < time_sampling; i++)
			{
				previous_time_layers_solutions[k_system][i].setSize (eq_matrixes[k_system].Size ()); // default starting solution point is 0
			}
		}

		// for T set only node functions 
		int area, iF;
		double coordinates[] = { 0.0, 0.0 };
		// go by elements
		for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
		{
			// set FCondition (area) for first 3 functions
			area = mesh_pointer->get_area (k_element);
			for (int i = 0; i < 3; i++)
			{
				iF = get_function_global_number (k_element, i);
				previous_time_layers_solutions[0][0].setElem (iF, function_starting_condition (0, coordinates, area));
			}
		}

		// make space for possible derivatives
		int dim = mesh_pointer->get_dimentionality ();
		derivative = new MathVector[dim]; // solutions are saved for each system
		for (int i = 0; i < dim; i++)
			derivative[i].setSize (eq_matrixes[0].Size ());
	}
	printf ("Melt_SF starting conditions set\n");
}

template<class TypeMesh>
void Task_Melt_SF<TypeMesh>::apply_first_boundary_conditions (int k_system)
{
	Task::apply_first_boundary_conditions (k_system);

	
	switch (k_system)
	{
		// for T add heating
	case 0:
	{
		if (time_layers[current_time_layer] > turn_on_source_time - 1e-10)
		{
			int boundary;
			int iF;
			MathVector * FCondition;
			int * functions;
			int f_amount;
			int k_element;
			int n1, n2;
			double coordinates[2];

			bool * replaced = new bool[N_functions];
			for (int i = 0; i < N_functions; i++)
			{
				replaced[i] = false;
			}

			// first condition by edges
			for (int i_edge = 0, i_end = mesh_pointer->edges->get_n_entries (); i_edge < i_end; i_edge++)
			{
				// get edges's nodes
				mesh_pointer->edges->get_edge_nodes (i_edge, &n1, &n2);
				// get boundary number for those nodes
				boundary = mesh_pointer->function_boundary_edge (n1, n2);

				if (boundary == 2)
				{
					// check if the edge lies within the heat source, if not, scratch that
					{
						mesh_pointer->get_node_coordinates (n1, coordinates);
						if (fabs (coordinates[0]) < L_source + SOURCE_PRECIS)
						{
							mesh_pointer->get_node_coordinates (n2, coordinates);
							if (fabs (coordinates[0]) < L_source + SOURCE_PRECIS)
							{
								// find element with that edge
								k_element = mesh_pointer->belonging_element (n1, n2);
								// get amount of values to put into f (from element)
								f_amount = mesh_pointer->get_amount_second_condition (k_element);
								FCondition = new MathVector (f_amount);

								functions = new int[f_amount];
								// get functions that are not 0 on the edge
								mesh_pointer->get_edge_functions (k_element, n1, n2, functions);

								// get first condition from element
								get_FCondition_edge (k_system, k_element, n1, n2, functions, 4, FCondition);

								// add it to f in places by functions
								for (int i = 0; i < f_amount; i++)
								{
									iF = get_function_global_number (k_element, functions[i]);
									if (iF != -1)
									{
										if (!replaced[iF])
										{
											replaced[iF] = true;
											eq_matrixes[k_system].clear_row (iF); // replace row
											eq_matrixes[k_system].set_entry (iF, iF, 1.0);
											eq_matrixes[k_system].set_f_entry (iF, FCondition->getElem (i));
										}
									}
								}
								delete FCondition;
								delete[] functions;
							}
						}
					}
				}
			}
			delete[] replaced;
		}
		break;
	}
	// for w and psi add 0 for all edges where T < Tpc
	case 1:
	case 2:
	{
		//go by elements
		// if T is smaller than T_phase_change in all def nodes, 0 all functions
		double coordinates[2];
		int counter;
		double T;
		int iF;
		int * nodes = new int[mesh_pointer->get_amount_of_def_nodes (0)];
		for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
		{
			counter = 0;
			// get def_nodes
			mesh_pointer->get_def_nodes (k_element, nodes);
			for (int i = 0; i < mesh_pointer->get_amount_of_def_nodes (0); i++)
			{
				// get coordinates
				mesh_pointer->get_node_coordinates (nodes[i], coordinates);
				// get T there
				get_solution_in_point (0, coordinates, &T);
				// if smaller, increase counter
				if (T < phys_T_phase_change - del_T_smoothing / 3.0)
					counter++;
			}
			if (counter == 3)

			//if (mesh_pointer->get_area (k_element) == 2)
			{
				// get functions
				int n_functions = mesh_pointer->get_amount_non_zero_functions (k_element);
				for (int i = 0; i < n_functions; i++)
				{
					iF = get_function_global_number (k_element, i);
					if (iF != -1)
					{
						eq_matrixes[k_system].clear_row (iF); // replace row
						eq_matrixes[k_system].set_entry (iF, iF, 1.0);
						eq_matrixes[k_system].set_f_entry (iF, 0.0);
					}
				}
			}
		}
		delete[] nodes;
		break;
	}
	}
}

template<class TypeMesh>
double Task_Melt_SF<TypeMesh>::function_FCondition (int k_system, double * coordinates, int area, int function)
{
	double r = 0.0;
	switch (k_system)
	{
	case 0:
	{
		//switch (boundary)
		//{
		//case 2: // don't panic, check file called conditions first
		//	r = (5.0 + 1.0) / (SCALE_TEMP);
		//	break;
		//case 3:
		//	r = phys_T0_S;
		//	break;
		//case 4:
		//	r = phys_T1;
		//	break;
		//}
		//break;
		switch (function)
		{
		case 1:
		{
			// heat
			r = phys_T0_F;
			break;
		}
		case 2:
		{
			// cool
			r = phys_T0_S;
			break;
		}
		case 4:
		{
			r = phys_T1;
			break;
		}
		}
		break;
	}
	case 1: // w
	{
		double step = 1e-2;
		double point[] = { 0.0, 0.0 };
		double c0[2], cn[2];
		mesh_pointer->get_0_boundaries (c0);
		mesh_pointer->get_N_boundaries (cn);
		switch (function)
		{
		case 0:
		{
			// anti-symmetry
			//return 0.0;
			point[0] = c0[0] + step;
			point[1] = coordinates[1];
			break;
		}
		case 1:
		{
			point[0] = cn[0] - step;
			point[1] = coordinates[1];
			break;
		}
		case 2:
		{
			point[0] = coordinates[0];
			point[1] = c0[1] + step;
			break;
		}
		case 3:
		{
			point[0] = coordinates[0];
			point[1] = cn[1] - step;
			break;
		}
		case 5:
		{
			// anti-symmetry
			return 0.0;
			break;
		}
		}

		double value = 0.0;
		get_non_linear_solution_in_point (2, point, &value);
		r = -2.0 * value / pow (step, 2.0);
		//r = 0.0;
		break;
	}
	case 2:
	{
		//switch (function)
		//{
		//case 5:
		//	r = 0.0;
		//	break;
		//}

		r = 0.0;
		break;
	}
	}
	return r;
}

//template<class TypeMesh>
//double Task_Melt_SF<TypeMesh>::function_FCondition (int k_system, double * coordinates, int area, int boundary)
//{
//	double x = coordinates[0];
//	double y = coordinates[1];
//	double t = time_layers[current_time_layer];
//
//	double r = 0.0;
//	switch (k_system)
//	{
//	case 0:
//	{
//		//r = x * x * y * y * t * t;
//		//r = x + y + t;
//		r = t + x * x + y * y;
//		break;
//	}
//	case 1:
//	{
//		//r = x * x + y * y + t;
//		//r = (x + y) * t;
//		r = t * (x + y);
//		break;
//	}
//	case 2:
//	{
//		//r = -x * x - y * y - t;
//		//r = t * (-x - y);
//		r = x + 2.0 * y + t;
//		break;
//	}
//	}
//	return r;
//}

template<class TypeMesh>
double Task_Melt_SF<TypeMesh>::function_f (int k_system, double * coordinates, int area)
{
	double x = coordinates[0];
	double y = coordinates[1];
	double t = time_layers[current_time_layer];

	double r = 0.0;
	switch (k_system)
	{
	case 0:
	{
		//r = 2.0 * x * x * y * y * t;
		////r *= (1.0 + x * x * y * y * t * t);
		//r *= (1.0 + exp (- pow (x * x * y * y * t * t - phys_T_phase_change, 2.0) / 0.01));
		//r += 4.0 * x * y * y * y * t * t;
		//r -= 4.0 * x * x * x * y * t * t;
		//r -= 2.0 * y * y * t * t;
		//r -= 2.0 * x * x * t * t;

		r = -3.0 + exp (-pow (t + x * x + y * y - phys_T_phase_change, 2.0) / del_T_smoothing);
		//r = -3.0;
		r += 4.0 * x;
		r -= 2.0 * y;
		//r = 1.0;
		break;
	}
	case 1:
	{
		//	r = - 4.0 + x * x + y * y + t;
			//r = t * (x + y);
		r = -x + y + t;
		break;
	}
	case 2:
	{
		//r = -8.0 * x * y + 3.0 - 2.0 * x * y * y * t * t;
		//r = - x - y - 1.0;
		r = -t * (x + y);

		break;
	}
	}
	return r;
}

template<class TypeMesh>
double Task_Melt_SF<TypeMesh>::get_lambda (double T)
{
	double r = 0.0;
	if (T < phys_T_phase_change - del_T_smoothing)
		r = lambda_solid;
	if (T > phys_T_phase_change + del_T_smoothing)
		r = lambda_fluid;
	if ((phys_T_phase_change - del_T_smoothing) < T && (T < phys_T_phase_change + del_T_smoothing))
		r = ((lambda_solid - lambda_fluid) / (2.0 * del_T_smoothing)) * T +
		+(del_T_smoothing * (lambda_fluid + lambda_solid) + phys_T_phase_change * (lambda_fluid - lambda_solid)) / (2.0 * del_T_smoothing);
	return r;
}

template<class TypeMesh>
void Task_Melt_SF<TypeMesh>::build_system (int k_system)
{
	// system is build fully for the whole domain
	// the only thing that is different if dirac
	// dirac is approximated and calculated over domain depending on T
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
		int area;

		Time_coef = new MathVector (time_sampling);
		get_time_coefficients (Time_coef);

		n_functions = mesh_pointer->get_amount_non_zero_functions (0);
		G = new Matrix (n_functions, n_functions);
		M = new Matrix (n_functions, n_functions);
		Psi = new MathVector (n_functions);
		DDFx = new Matrix (n_functions, n_functions);
		DDFy = new Matrix (n_functions, n_functions);
		Q = new MathVector (n_functions);
		MQ = new MathVector (n_functions);
		B = new MathVector (n_functions);
		double * coordinates = new double[dim];

		int n_points = mesh_pointer->amount_of_integration_points (0);
		double * weights = new double[n_points];
		double ** points = new double *[n_points];
		for (int k = 0; k < n_points; k++)
			points[k] = new double[2];

		double coef = 1.0 / (PR);
		//double coef = 1.0 / (PR * sqrt (GR));

		double lambda;
		double T_value;

		for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
		{
			n_functions_cur = mesh_pointer->get_amount_non_zero_functions (k_element);
			if (n_functions != n_functions_cur)
			{
				delete G;
				delete M;
				delete Q;
				delete MQ;
				delete Psi;
				delete DDFx;
				delete B;
				delete DDFy;

				n_functions = n_functions_cur;

				G = new Matrix (n_functions, n_functions);
				M = new Matrix (n_functions, n_functions);
				Psi = new MathVector (n_functions);
				DDFx = new Matrix (n_functions, n_functions);
				DDFy = new Matrix (n_functions, n_functions);
				Q = new MathVector (n_functions);
				MQ = new MathVector (n_functions);
				B = new MathVector (n_functions);
			}

			// get G, M matrices 
			mesh_pointer->get_G_local_matrix (k_element, G);
			mesh_pointer->get_M_local_matrix (k_element, M);

			area = mesh_pointer->get_area (k_element);
			// lambda as average
			lambda = 0.0;
			double sum = 0;
			for (int i = 0, i_end = mesh_pointer->get_amount_of_def_nodes (k_element); i < i_end; i++)
			{
				mesh_pointer->get_node_coordinates (mesh_pointer->get_node_number (k_element, i), coordinates);
				get_non_linear_solution_in_point (0, coordinates, &T_value);
				lambda = get_lambda (T_value);
				sum += lambda;
			}
			lambda = sum / (double)mesh_pointer->get_amount_of_def_nodes (k_element);

			// get psi values
			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (k_element, i);
				if (iF != -1)
				{
					Psi->setElem (i, non_linear_layers_solutions[2].getElem (iF));
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
							eq_matrixes[k_system].add_to_entry (iF, jF, G->Elem (i, j) * coef * lambda / lambda_fluid);
							// add M (by t) into respective places (by def_nodes)
							eq_matrixes[k_system].add_to_entry (iF, jF, M->Elem (i, j) * Time_coef->getElem (0));
						}
					}
				}
			}

			// add derivatives
			for (int k = 0; k < n_functions; k++)
			{
				mesh.elements[k_element]->get_DDF (k, 1, 0, DDFy);
				mesh.elements[k_element]->get_DDF (k, 0, 1, DDFx);

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
			if (Ste < 1e+9)
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
									eq_matrixes[k_system].add_to_entry (iF, jF, M->Elem (i, j) * Time_coef->getElem (0) / Ste);
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
								eq_matrixes[k_system].add_to_f_entry (iF, MQ->getElem (i) * Time_coef->getElem (t) / Ste);
							}
						}
					}
				}
			}
		}

		delete Time_coef;
		delete G;
		delete M;
		delete Q;
		delete MQ;
		delete Psi;
		delete DDFx;
		delete B;
		delete DDFy;
		delete[] weights;
		for (int k = 0; k < n_points; k++)
			delete[] points[k];
		delete[] points;
		delete[] coordinates;
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
		MathVector * TD;

		int n_functions;
		int n_functions_cur;
		int iF, jF;
		int dim = mesh_pointer->get_dimentionality ();
		int area;

		Time_coef = new MathVector (time_sampling);
		get_time_coefficients (Time_coef);

		n_functions = mesh_pointer->get_amount_non_zero_functions (0);
		G = new Matrix (n_functions, n_functions);
		M = new Matrix (n_functions, n_functions);
		Psi = new MathVector (n_functions);
		DDFx = new Matrix (n_functions, n_functions);
		DDFy = new Matrix (n_functions, n_functions);
		Q = new MathVector (n_functions);
		MQ = new MathVector (n_functions);
		T = new MathVector (n_functions);
		TD = new MathVector (n_functions);
		D = new Matrix (n_functions, n_functions);
		B = new MathVector (n_functions);

		double * coordinates = new double[dim];
		double coef = GR;
		double inversion_coef;
		double T_value;
		//double coef = 1.0 / sqrt (GR);
		for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
		{
			n_functions_cur = mesh_pointer->get_amount_non_zero_functions (k_element);
			if (n_functions != n_functions_cur)
			{
				delete G;
				delete M;
				delete Q;
				delete MQ;
				delete Psi;
				delete DDFx;
				delete DDFy;
				delete T;
				delete TD;
				delete D;
				delete B;

				n_functions = n_functions_cur;

				G = new Matrix (n_functions, n_functions);
				M = new Matrix (n_functions, n_functions);
				Psi = new MathVector (n_functions);
				DDFx = new Matrix (n_functions, n_functions);
				DDFy = new Matrix (n_functions, n_functions);
				Q = new MathVector (n_functions);
				MQ = new MathVector (n_functions);
				T = new MathVector (n_functions);
				TD = new MathVector (n_functions);
				D = new Matrix (n_functions, n_functions);
				B = new MathVector (n_functions);
			}

			area = mesh_pointer->get_area (k_element);
			// get G, M matrices 
			mesh_pointer->get_G_local_matrix (k_element, G);
			mesh_pointer->get_M_local_matrix (k_element, M);

			// get psi values
			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (k_element, i);
				if (iF != -1)
				{
					Psi->setElem (i, non_linear_layers_solutions[2].getElem (iF));
				}
			}

			// get T values
			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (k_element, i);
				if (iF != -1)
				{
					T->setElem (i, non_linear_layers_solutions[0].getElem (iF));
				}
			}

			// get dT/dx
			mesh.elements[k_element]->get_D (0, D);
			D->MultiplyMatrixByVector (*T, TD);

			// average on the element
			if (density_inversion)
			{
				double sum = 0;
				for (int i = 0, i_end = mesh_pointer->get_amount_of_def_nodes (k_element); i < i_end; i++)
				{
					mesh_pointer->get_node_coordinates (mesh_pointer->get_node_number (k_element, i), coordinates);
					get_non_linear_solution_in_point (0, coordinates, &T_value);
					inversion_coef = density_inversion_coef (T_value);
					sum += inversion_coef;
				}
				inversion_coef = sum / (double)mesh_pointer->get_amount_of_def_nodes (k_element);
			}
			else
				inversion_coef = 1.0;

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
							//eq_matrixes[k_system].add_to_entry (iF, jF, coef * G->Elem (i, j));
							eq_matrixes[k_system].add_to_entry (iF, jF, G->Elem (i, j));
							// add M (by t) into respective places (by def_nodes)
							eq_matrixes[k_system].add_to_entry (iF, jF, M->Elem (i, j) * Time_coef->getElem (0));
						}
					}

					// add dt/dx
					//if (!density_inversion)
						eq_matrixes[k_system].add_to_f_entry (iF, inversion_coef * coef * TD->getElem (i));
					//if (density_inversion)
					//{
					//	// sum by derivatives
					//	for (int j = 0; j < n_functions; j++)
					//	{
					//		eq_matrixes[k_system].add_to_f_entry (iF, D->Elem (i, j) * GR * inversion_coef / SCALE_TEMP);
					//	}
					//}
					//eq_matrixes[k_system].add_to_f_entry (iF, TD->getElem (i));
				}
			}

			// add derivatives
			for (int k = 0; k < n_functions; k++)
			{
				mesh.elements[k_element]->get_DDF (k, 1, 0, DDFy);
				mesh.elements[k_element]->get_DDF (k, 0, 1, DDFx);

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
		delete[] coordinates;
		delete G;
		delete M;
		delete Q;
		delete MQ;
		delete Psi;
		delete DDFx;
		delete DDFy;
		delete T;
		delete TD;
		delete D;
		delete B;
		break;
	}
	case 2:
	{
		Matrix * G;
		Matrix * M;
		MathVector * W;
		MathVector * B;
		MathVector * WM;
		int n_functions;
		int n_functions_cur;
		int iF, jF;
		int area;
		int dim = mesh_pointer->get_dimentionality ();

		n_functions = mesh_pointer->get_amount_non_zero_functions (0);
		G = new Matrix (n_functions, n_functions);
		M = new Matrix (n_functions, n_functions);
		W = new MathVector (n_functions);
		WM = new MathVector (n_functions);
		B = new MathVector (n_functions);

		for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
		{
			n_functions_cur = mesh_pointer->get_amount_non_zero_functions (k_element);
			if (n_functions != n_functions_cur)
			{
				delete G;
				delete M;
				delete W;
				delete WM;
				delete B;

				n_functions = n_functions_cur;

				G = new Matrix (n_functions, n_functions);
				M = new Matrix (n_functions, n_functions);
				W = new MathVector (n_functions);
				WM = new MathVector (n_functions);
				B = new MathVector (n_functions);
			}

			area = mesh_pointer->get_area (k_element);
			// get G, M matrices 
			mesh_pointer->get_G_local_matrix (k_element, G);
			mesh_pointer->get_M_local_matrix (k_element, M);

			// get W
			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (k_element, i);
				if (iF != -1)
				{
					W->setElem (i, non_linear_layers_solutions[1].getElem (iF));
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

					// put W into respective places (by def_nodes)
					eq_matrixes[k_system].add_to_f_entry (iF, WM->getElem (i));
				}
			}
		}
		delete G;
		delete M;
		delete W;
		delete WM;
		delete B;
		break;
	}
	}
}

//template<class TypeMesh>
//bool Task_Melt_SF<TypeMesh>::solve_task (int method, int d_type, int depth)
//{
//	if (!prepared)
//	{
//		printf ("ERROR: missing preparation\n");
//		return false;
//	}
//
//	non_linear_layers_solutions = new MathVector[n_systems];
//	for (int k_system = 0; k_system < n_systems; k_system++)
//	{
//		non_linear_layers_solutions[k_system].setSize (eq_matrixes[k_system].Size ()); // default starting solution point is 0
//	}
//
//	double * discr = new double[n_systems];
//	double * prev_discr = new double[n_systems];
//	bool stable = false; // indicator of reaching stable condition
//	bool found = false;
//	int solver_iterations;
//	double sum_discr;
//	int unchanged_discr;
//	for (current_time_layer = 1; current_time_layer < n_time_layers && !stable; current_time_layer++)
//	{
//		printf ("current time: %.7lf\n", time_layers[current_time_layer]);
//		// save previous time layers solution for non-linearity
//		for (int k_system = 0; k_system < n_systems; k_system++)
//		{
//			non_linear_layers_solutions[k_system].Copy ((previous_time_layers_solutions[k_system][0]));
//		}
//		found = false;
//		for (int i = 0; i < MAX_NONLINEAR_ITER_MELT && !found; i++)
//		{
//			// iterate through each system of equations
//			for (int k_system = 0; k_system < n_systems; k_system++)
//			{
//				eq_matrixes[k_system].Clear ();
//				eq_matrixes[k_system].set_starting_point (&(non_linear_layers_solutions[k_system])); // set last layer solution as x0 for solver
//				// fprintf starting vector
//				//char name[256];
//				//sprintf (name, "Result//Melt_SF//Matrix//s%i_sv.txt", k_system);
//				//non_linear_layers_solutions[k_system].FPrint (name);
//
//				build_system (k_system);  // reset system for current layer
//				apply_boundary_conditions (k_system); // apply boundary conditions
//				solver_iterations = eq_matrixes[k_system].solve_LOS (0); // solve it	
//				discr[k_system] = relax (k_system, &(non_linear_layers_solutions[k_system]));
//			}
//			if (i > 0)
//			{
//				sum_discr = 0.0;
//				unchanged_discr = 0;
//				for (int k_system = 0; k_system < n_systems; k_system++)
//				{
//					sum_discr += discr[k_system];
//					if (fabs (discr[k_system] - prev_discr[k_system]) < DISCR_CHANGE)
//						unchanged_discr++;
//				}
//				if (sum_discr / 3.0 < NONLINEAR_DISCR)
//					found = true;
//				if (unchanged_discr == 3)
//					found = true;
//
//			}
//			for (int k_system = 0; k_system < n_systems; k_system++)
//			{
//				prev_discr[i] = discr[i];
//			}
//		}
//		// move layers
//		for (int k_system = 0; k_system < n_systems; k_system++)
//			next_time_layer (k_system, &(non_linear_layers_solutions[k_system]));
//		//print_solutions ();
//		//print_extra_data ();
//	}
//	current_time_layer--;
//
//	delete[] discr;
//	delete[] prev_discr;
//	delete[] non_linear_layers_solutions;
//	non_linear_layers_solutions = NULL;
//	return true;
//}
//
//template<class TypeMesh>
//bool Task_Melt_SF<TypeMesh>::solve_task (int solver_param[][5])
//{
//	// TEST
//	if (!prepared)
//	{
//		printf ("ERROR: missing preparation\n");
//		return false;
//	}
//
//	non_linear_layers_solutions = new MathVector[n_systems];
//	for (int k_system = 0; k_system < n_systems; k_system++)
//	{
//		non_linear_layers_solutions[k_system].setSize (eq_matrixes[k_system].Size ()); // default starting solution point is 0
//	}
//
//	double * discr = new double[n_systems];
//	double * prev_discr = new double[n_systems];
//	bool stable = false; // indicator of reaching stable condition
//	bool found = false;
//	int solver_iterations;
//	double sum_discr;
//	double T_change;
//	int unchanged_discr;
//	FILE * file_T = fopen ("Result//Melt_SF//T_change.txt", "w");
//	for (current_time_layer = 1; current_time_layer < n_time_layers && !stable; current_time_layer++)
//	{
//		printf ("current time: %.7lf\n", time_layers[current_time_layer]);
//		//printf ("%.7lf ", time_layers[current_time_layer]);
//		// save previous time layers solution for non-linearity
//		for (int k_system = 0; k_system < n_systems; k_system++)
//		{
//			non_linear_layers_solutions[k_system].Copy ((previous_time_layers_solutions[k_system][0]));
//		}
//		found = false;
//		for (int i = 0; i < MAX_NONLINEAR_ITER_MELT && !found; i++)
//		{
//			printf ("\tnonlinear iteration: %i\n", i);
//			// iterate through each system of equations
//			for (int k_system = 0; k_system < n_systems; k_system++)
//			{
//				eq_matrixes[k_system].Clear ();
//				//eq_matrixes[k_system].set_starting_point (&(non_linear_layers_solutions[k_system])); // set last layer solution as x0 for solver
//				eq_matrixes[k_system].set_starting_point (&(previous_time_layers_solutions[k_system][0])); // set last layer solution as x0 for solver
//																									 // fprintf starting vector
//				//char name[256];
//				//sprintf (name, "Result//Melt_SF//Matrix//s%i_sv.txt", k_system);
//				//non_linear_layers_solutions[k_system].FPrint (name);
//
//				build_system (k_system);  // reset system for current layer
//				apply_boundary_conditions (k_system); // apply boundary conditions
//				solver_iterations = eq_matrixes[k_system].solve (solver_param[k_system]); // solve it	
//				discr[k_system] = relax (k_system, &(non_linear_layers_solutions[k_system]));
//			}
//			if (i > 0)
//			{
//				sum_discr = 0.0;
//				unchanged_discr = 0;
//				for (int k_system = 0; k_system < n_systems; k_system++)
//				{
//					sum_discr += discr[k_system];
//					if (fabs (discr[k_system] - prev_discr[k_system]) < DISCR_CHANGE_MELT)
//						unchanged_discr++;
//				}
//				if (sum_discr / 3.0 < NONLINEAR_DISCR_MELT)
//					found = true;
//				if (unchanged_discr >= n_systems)
//					found = true;
//				T_change = non_linear_layers_solutions[0].sum_sq_dif (previous_time_layers_solutions[0][0]);
//				printf ("%i: %e %e %e\tT: %e %e\n", unchanged_discr, discr[0] - prev_discr[0], discr[1] - prev_discr[1], discr[2] - prev_discr[2], T_change, T_change / non_linear_layers_solutions[0].Norm ());
//				fprintf (file_T, "%lf %e %e\n", time_layers[current_time_layer], T_change, T_change / non_linear_layers_solutions[0].Norm ());
//			}
//			for (int k_system = 0; k_system < n_systems; k_system++)
//			{
//				prev_discr[k_system] = discr[k_system];
//			}
//		}
//		// move layers
//		for (int k_system = 0; k_system < n_systems; k_system++)
//			next_time_layer (k_system, &(non_linear_layers_solutions[k_system]));
//		//calc_gradient (0);
//		//calc_gradient (2);
//		//print_solutions ();
//		//print_extra_data ();
//	}
//	current_time_layer--;
//	fclose (file_T);
//
//	//FILE * file_source_power = fopen ("Result//Melt_SF//source_power.txt", "w");
//	//fprintf (file_source_power, "%lf\n", source_power);
//	//source_power *= time_layers[current_time_layer] / n_time_layers;
//	//source_power *= lambda_fluid;
//	//source_power *= SCALE_TEMP / (L_source * 2.0 * (SCALE_SIZE / 1000) * time_layers[current_time_layer] * SCALE_TIME);
//	//fprintf (file_source_power, "%lf\n", source_power);
//	//printf ("Source power: %lf\n", source_power);
//	//fclose (file_source_power);
//
//	delete[] discr;
//	delete[] prev_discr;
//	delete[] non_linear_layers_solutions;
//	non_linear_layers_solutions = NULL;
//	return true;
//}

template<class TypeMesh>
void Task_Melt_SF<TypeMesh>::get_scales ()
{
	FILE * file = fopen ("Source Files//Melt_SF//Scales.txt", "r");
	fscanf (file, "%lf", &SCALE_TEMP);
	fscanf (file, "%lf", &SCALE_TIME);
	fscanf (file, "%lf", &SCALE_SIZE);
	fscanf (file, "%lf", &SCALE_SPEED);
	fclose (file);
}

template<class TypeMesh>
void Task_Melt_SF<TypeMesh>::get_scales (char * scales_file_name)
{
	FILE * file = fopen (scales_file_name, "r");
	fscanf (file, "%lf", &SCALE_TEMP);
	fscanf (file, "%lf", &SCALE_TIME);
	fscanf (file, "%lf", &SCALE_SIZE);
	fscanf (file, "%lf", &SCALE_SPEED);
	fclose (file);
}

template<class TypeMesh>
void Task_Melt_SF<TypeMesh>::use_density_inversion ()
{
	density_inversion = true;
}

template<class TypeMesh>
void Task_Melt_SF<TypeMesh>::set_starting_conditions_from_task (Task_pointer * task)
{
	flag_calc_start = false;
	
	// assume the same mesh 

	for (int k_system = 0; k_system < n_systems; k_system++)
	{
		for (int i = 0; i < time_sampling; i++)
		{
			previous_time_layers_solutions[k_system][i].setSize (eq_matrixes[k_system].Size ()); // default starting solution point is 0
		}
		task->eq_matrixes[k_system].get_solution (&previous_time_layers_solutions[k_system][0]);
	}


	// make space for possible derivatives
	int dim = mesh_pointer->get_dimentionality ();
	derivative = new MathVector[dim]; // solutions are saved for each system
	for (int i = 0; i < dim; i++)
		derivative[i].setSize (eq_matrixes[0].Size ());
}

template<class TypeMesh>
void Task_Melt_SF<TypeMesh>::reset_parameters (double _SCALE_TEMP, double _GR, double _Ste, double _phys_T_phase_change, double _phys_T0_F)
{
	GR = _GR;
	Ste = _Ste;
	phys_T_phase_change = _phys_T_phase_change;
	phys_T0_F = _phys_T0_F;
	SCALE_TEMP = _SCALE_TEMP;
}

template<class TypeMesh>
double Task_Melt_SF<TypeMesh>::get_front_position (double x)
{
	double p[2];
	mesh_pointer->get_N_boundaries (p);
	double step = 0.1;
	p[0] = x;
	double y = p[1];
	double value = 0.0;

	for (int i = 0; i < 1000 && fabs (value - phys_T_phase_change) > 1e-12; i++)
	{
		// go down by step
		p[1] = y - step;
		get_solution_in_point (0, p, &value);
		// if value is bigger than T_phys_change
		if (value + 1e-15 > phys_T_phase_change)
		{
			// step = step / 2
			step = step / 2.0;
		}
		else
		{
			// else new y
			y = p[1];
		}
	}
	

	return y;
}

template<class TypeMesh>
void Task_Melt_SF<TypeMesh>::print_solutions ()
{
	// save pics for T, w, psi 
	// over time_layers
	// colored for now
	// for T add Tpc isoline specifically
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
		for (int k_system = 0; k_system < 3; k_system++)
		{
			// fprint solution
			for (size_t i = 0, i_end = local_time_stamps.size (); i < i_end; i++)
			{
				get_time_approx (local_time_stamps[i], c);
				printf ("\tsaving solutions into files, time layer:\t%.16lf\n", local_time_stamps[i]);
				sprintf (name, "Result//Melt_SF//Solution//s%i_t_%.5lf.txt", k_system, local_time_stamps[i]);
				file = fopen (name, "w");
				double val;
				for (int j = 0, j_end = previous_time_layers_solutions[k_system][0].getSize (); j < j_end; j++)
				{
					val = 0;
					for (int t = 0; t < time_sampling; t++)
					{
						val += c->getElem (t) * previous_time_layers_solutions[k_system][t].getElem (j);
					}
					fprintf (file, "%.16lf\n", val);
				}
				fclose (file);

				// fprint matrix
				//sprintf (name, "Result//Melt_SF//Matrix//s%i_t_%.5lf.txt", k_system, local_time_stamps[i]);
				//eq_matrixes[k_system].fprint (name);

			}

			for (int a = 1; a < 2/* && k_system == 0*/; a++)
			{
				char pic_name[128];
				wchar_t wtemp[128];

				// draw solution
				for (size_t i = 0, i_end = local_time_stamps.size (); i < i_end; i++)
				{
					painter->set_max_resolution (1000);
					painter->set_axis_names ("X", "Y");
					//painter->add_area_dividers ();
					//painter->draw_field (0, COLOR_SCALE_BLACK_WHITE);
					painter->draw_contour_lines (k_system, 15);
					//painter->draw_contour_lines (k_system, 10);
					painter->set_legend_full_iso ();
					//painter->no_legend ();
					painter->set_iso_line_width (1.0);
					painter->set_axes_scale (10, 10);

					painter->set_area (a);	
					switch (k_system)
					{
					case 0:
						sprintf (pic_name, "t = %.5lf", local_time_stamps[i]);
						//set_symmetry (1);
						painter->add_isoline (phys_T_phase_change);
						if (a == 1)
						{
							painter->add_isoline ((4.029 - TEMP_T0) / (SCALE_TEMP));
							painter->set_min_max (0.0, 1.0);
							painter->draw_contour_lines (k_system, 20);
							//painter->no_legend ();
							painter->draw_field (0, COLOR_SCALE_RED_YELLOW);
							//double s0[] = {-0.4, 0.0};
							//double sN[] = {0.4, 1.0};
							//painter->select_section (s0, sN);
							sprintf (name, "Pictures//Melt_SF//s%i_t_%.5lf.png", k_system, local_time_stamps[i]);
						}
						else
						{
							//painter->set_min_max (phys_T0_S, phys_T1);
							painter->set_min_max (0.0, 0.4);
							painter->set_scale (1.0, 2.5);
							painter->set_axes_scale (17, 6);
							painter->set_iso_line_width (2.0);
							//painter->draw_contour_lines (k_system, 10);
							painter->draw_contour_lines (k_system, 12, COLOR_SCALE_RAINBOW);
							sprintf (name, "Pictures//Melt_SF//s%i_top_t_%.5lf.png", k_system, local_time_stamps[i]);
						}
						break;
					case 1:
						painter->draw_contour_lines (k_system, 50, COLOR_SCALE_RAINBOW);
						painter->set_min_max (-1500.0, 1500.0);
						//painter->no_legend ();
						sprintf (pic_name, "t = %.5lf", local_time_stamps[i]);
						//set_symmetry (2);
						sprintf (name, "Pictures//Melt_SF//s%i_t_%.5lf.png", k_system, local_time_stamps[i]);
						break;
					case 2:
						painter->draw_contour_lines (k_system, 40, COLOR_SCALE_RAINBOW);
						painter->set_min_max (-20.0, 20.0);
						//painter->no_legend ();
						sprintf (pic_name, "t = %.5lf", local_time_stamps[i]);
						//set_symmetry (2);
						sprintf (name, "Pictures//Melt_SF//s%i_t_%.5lf.png", k_system, local_time_stamps[i]);
						break;
					}
					painter->set_picture_name (pic_name);
					mbstowcs (wtemp, name, strlen (name) + 1);
					std::wstring w_name = wtemp;

					if (k_system == 0 || a == 1)
						painter->draw_to_file (w_name);
					painter->reset ();

					//if ((k_system == 0))
					//{
					//	// gradient pic 
					//	//painter->add_area_dividers ();
					//	painter->set_max_resolution (800);
					//	painter->set_min_max (0.0, 1.0);
					//	painter->draw_contour_lines (k_system, 10);
					//	painter->set_legend_full_iso ();
					//	painter->set_iso_line_width (1.0);

					//	painter->set_axes_scale (8, 10);
					//	painter->draw_antigradient_field (0);
					//	//painter->set_vector_net (20, 20);
					//	painter->set_vector_net (30, 30);
					//	double c0[] = {-0.5, 0.0};
					//	double cN[] = { 0.5, 1.0};
					//	painter->select_section (c0, cN);
					//	//painter->no_legend ();
					//	//painter->add_isoline (phys_T_phase_change);
					//	painter->add_isoline ((4.029 - TEMP_T0) / (SCALE_TEMP));
					//	set_symmetry (1);
					//	painter->set_picture_name (" ");
					//	sprintf (name, "Pictures//Melt_SF//grad_T_t_%.5lf.png", local_time_stamps[i]);
					//	mbstowcs (wtemp, name, strlen (name) + 1);
					//	w_name = wtemp;
					//	painter->draw_to_file (w_name);
					//	painter->reset ();
					//}
				}
			}


			//// print solution
			//for (size_t i = 0, i_end = local_time_stamps.size (); i < i_end; i++)
			//{
			//	char name[128];
			//	//sprintf (name, "Result//Melt_SF//s%i_t_%.5lf.txt", k_system, local_time_stamps[i]);
			//	sprintf (name, "Result//Plumes//s%i_t_%.5lf.txt", k_system, local_time_stamps[i]);

			//	file = fopen (name, "w");
			//	previous_time_layers_solutions[k_system][0].FPrint (file);
			//	fclose (file);
			//}
		}

		delete c;
	}
}

template<class TypeMesh>
double Task_Melt_SF<TypeMesh>::Dirac_function (double T, double * coordinates)
{
	double r = 0.0;
	double t = time_layers[current_time_layer];
	double x = coordinates[0];
	double y = coordinates[1];
	//if (fabs (T - phys_T_phase_change) < del_T_smoothing + 1e-15)
	//{
	//r = 1.0 / (2.0 * del_T_smoothing);
	//	r *= (1.0 - pow (T - phys_T_phase_change, 2.0) / pow (del_T_smoothing, 2.0));
	//}
	r = exp (-pow (T - phys_T_phase_change, 2.0) / del_T_smoothing);
	//r = 0.0;
	//r = T * T;
	return r;
}

template<class TypeMesh>
double Task_Melt_SF<TypeMesh>::density_inversion_coef (double T1)
{
	double T = T1 * SCALE_TEMP + TEMP_T0;
	double T0 = 4.029325;
	double q = 1.894816;
	if (T - T0 > 0)
		return q * pow (T - T0, q - 1.0);
	return -q * pow (T0 - T, q - 1.0);

	//double T = T1 * SCALE_TEMP + TEMP_T0;
	//double T0 = 4.029325;
	//double q = 1.894816;
	//if (T - T0 > 0)
	//	return pow (T - T0, q - 1.0);
	//return pow (T0 - T, q - 1.0);
}

template<class TypeMesh>
void Task_Melt_SF<TypeMesh>::next_time_layer (int k_system, MathVector * solution)
{
	for (int j = time_sampling - 1; j > 0; j--)
	{
		previous_time_layers_solutions[k_system][j].Copy ((previous_time_layers_solutions[k_system][j - 1]));
	}

	previous_time_layers_solutions[k_system][0].Copy (*solution);
}

template<class TypeMesh>
void Task_Melt_SF<TypeMesh>::print_extra_data ()
{
	// save pics for T, w, psi 
	// over time_layers
	// colored for now
	// for T add Tpc isoline specifically
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
		char file_name[128];
		char pic_name[128];
		char name[128];
		wchar_t wtemp[128];

		for (size_t i = 0, i_end = local_time_stamps.size (); i < i_end; i++)
		{
			double weights[] = { 1.0, 1.0, 0.0 };

			//// local heat flux on the boundary
			//weights[0] = 1.0;
			//weights[1] = -1.0;					
			//set_symmetry (1);
			//sprintf (file_name, "Result//Melt_SF//hft_t_%.5lf.txt", local_time_stamps[i]);
			//save_slice_first_derivative_weighted (0, 1, weights, "Source Files//Melt_SF//bound_slice.txt", file_name, true);

			//// T vertical over the heat source
			//weights[0] = SCALE_SIZE;
			//weights[1] = SCALE_TEMP;
			//weights[2] = TEMP_T0;
			//set_symmetry (0);
			////set_symmetry (1);
			//sprintf (file_name, "Result//Melt_SF//T_source_axis_%.5lf.txt", local_time_stamps[i]);
			//save_slice_weighted (0, weights, "Source Files//Melt_SF//Slice_T_vert.txt", file_name, true);
			//weights[2] = 0.0;

			//// T vertical over the heat source
			//weights[0] = SCALE_SIZE;
			//weights[1] = SCALE_TEMP;
			//set_symmetry (1);
			//sprintf (file_name, "Result//Melt_SF//T_source_axis_dim_%.5lf.txt", local_time_stamps[i]);
			//save_slice_weighted (0, weights, "Source Files//Melt_SF//Slice_T_vert.txt", file_name, true);

			// T horizontal
			//if ((fabs (local_time_stamps[i] - 0.026) < 1e-7) || (fabs (local_time_stamps[i] - 0.05) < 1e-7)
			//	|| (fabs (local_time_stamps[i] - 0.1) < 1e-7)
			//	)
			//if ((local_time_stamps[i] > 0.04 - 1e-10) && (local_time_stamps[i] < 0.5 + 1e-10))
			//{
			//	weights[0] = SCALE_SIZE;
			//	weights[1] = SCALE_TEMP;
			//	set_symmetry (1);

			//	sprintf (file_name, "Result//Melt_SF//T_horizontal1_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_weighted (0, weights, "Source Files//Melt_SF//slice_Vx1.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//T_horizontal2_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_weighted (0, weights, "Source Files//Melt_SF//slice_Vx2.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//T_horizontal3_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_weighted (0, weights, "Source Files//Melt_SF//slice_Vx3.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//T_horizontal4_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_weighted (0, weights, "Source Files//Melt_SF//slice_Vx4.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//T_horizontal5_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_weighted (0, weights, "Source Files//Melt_SF//slice_Vx5.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//T_horizontal6_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_weighted (0, weights, "Source Files//Melt_SF//slice_Vx6.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//T_horizontal7_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_weighted (0, weights, "Source Files//Melt_SF//slice_Vx7.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//T_horizontal8_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_weighted (0, weights, "Source Files//Melt_SF//slice_Vx8.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//T_horizontal9_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_weighted (0, weights, "Source Files//Melt_SF//slice_Vx9.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//T_horizontal10_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_weighted (0, weights, "Source Files//Melt_SF//slice_Vx10.txt", file_name, true);
			//}
			// T gradient pic
			//painter->add_area_dividers ();
			//painter->set_max_resolution (600);
			//painter->draw_contour_lines (0, 10);
			//painter->set_iso_line_width (1.0);
			//painter->set_axes_scale (6, 3);
			//painter->set_min_max (phys_T0_S, phys_T1);
			//painter->draw_antigradient_field (0);
			//painter->set_vector_net (20, 10);
			//set_symmetry (1);
			//double s0[] = {-1.2, 0.0};
			//double sN[] = {1.2, 1.2};
			//painter->select_section (s0, sN);
			//painter->set_picture_name (" ");
			//sprintf (name, "Pictures//Melt_SF//grad_T_t_%.5lf.png", local_time_stamps[i]);
			//mbstowcs (wtemp, name, strlen (name) + 1);
			//std::wstring w_name = wtemp;
			//painter->draw_to_file (w_name);
			//painter->reset ();

			//// local heat flux on the bottom
			//weights[0] = 1.0;
			//weights[1] = -1.0;
			//set_symmetry (1);
			//sprintf (file_name, "Result//Melt_SF//hfb_t_%.5lf.txt", local_time_stamps[i]);
			//save_slice_first_derivative_weighted (0, 1, weights, "Source Files//Melt_SF//bound_slice_bottom.txt", file_name, true);

			////// local heat flux on the bottom weighted
			//weights[0] = SCALE_SIZE;
			//weights[1] = -lambda_fluid * SCALE_TEMP / SCALE_SIZE;
			//set_symmetry (1);
			//sprintf (file_name, "Result//Melt_SF//hfb_dim_t_%.5lf.txt", local_time_stamps[i]);
			//save_slice_first_derivative_weighted (0, 1, weights, "Source Files//Melt_SF//bound_slice_bottom.txt", file_name, true);

			//// Vy = dPsi/dx
			//weights[0] = SCALE_SIZE;
			//weights[1] = -SCALE_SPEED;
			//set_symmetry (0);
			////set_symmetry (1);
			//sprintf (file_name, "Result//Melt_SF//speed_y_t_%.5lf.txt", local_time_stamps[i]);
			//save_slice_first_derivative_weighted (2, 0, weights, "Source Files//Melt_SF//slice_Vy.txt", file_name, true);

			//// Vx = -dPsi/dy
			////if ((fabs (local_time_stamps[i] - 0.004) < 1e-7) || (fabs (local_time_stamps[i] - 0.005) < 1e-7))
			//if ((local_time_stamps[i] > 0.04 - 1e-10) && (local_time_stamps[i] < 0.5 + 1e-10))
			//{
			//	weights[0] = SCALE_SIZE;
			//	weights[1] = -SCALE_SPEED;
			//	set_symmetry (1);

			//	sprintf (file_name, "Result//Melt_SF//speed_x_t_horizontal1_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_first_derivative_weighted (2, 1, weights, "Source Files//Melt_SF//slice_Vy1.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//speed_x_t_horizontal2_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_first_derivative_weighted (2, 1, weights, "Source Files//Melt_SF//slice_Vy2.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//speed_x_t_horizontal3_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_first_derivative_weighted (2, 1, weights, "Source Files//Melt_SF//slice_Vy3.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//speed_x_t_horizontal4_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_first_derivative_weighted (2, 1, weights, "Source Files//Melt_SF//slice_Vy4.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//speed_x_t_horizontal5_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_first_derivative_weighted (2, 1, weights, "Source Files//Melt_SF//slice_Vy5.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//speed_x_t_horizontal6_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_first_derivative_weighted (2, 1, weights, "Source Files//Melt_SF//slice_Vy6.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//speed_x_t_horizontal7_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_first_derivative_weighted (2, 1, weights, "Source Files//Melt_SF//slice_Vy7.txt", file_name, true);
			//}
			//// Vy = dPsi/dx
			//if ((local_time_stamps[i] > 0.04 - 1e-10) && (local_time_stamps[i] < 0.5 + 1e-10))
			//{
			//	weights[0] = 1.0;
			//	weights[1] = -1.0;
			//	set_symmetry (1);
			//	sprintf (file_name, "Result//Melt_SF//speed_y_t_horizontal1_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_first_derivative_weighted (2, 0, weights, "Source Files//Melt_SF//slice_Vy1.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//speed_y_t_horizontal2_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_first_derivative_weighted (2, 0, weights, "Source Files//Melt_SF//slice_Vy2.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//speed_y_t_horizontal3_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_first_derivative_weighted (2, 0, weights, "Source Files//Melt_SF//slice_Vy3.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//speed_y_t_horizontal4_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_first_derivative_weighted (2, 0, weights, "Source Files//Melt_SF//slice_Vy4.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//speed_y_t_horizontal5_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_first_derivative_weighted (2, 0, weights, "Source Files//Melt_SF//slice_Vy5.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//speed_y_t_horizontal6_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_first_derivative_weighted (2, 0, weights, "Source Files//Melt_SF//slice_Vy6.txt", file_name, true);
			//	sprintf (file_name, "Result//Melt_SF//speed_y_t_horizontal7_%.5lf.txt", local_time_stamps[i]);
			//	save_slice_first_derivative_weighted (2, 0, weights, "Source Files//Melt_SF//slice_Vy7.txt", file_name, true);

			//	//weights[0] = SCALE_SIZE;
			//	//weights[1] = -SCALE_SPEED;
			//	//set_symmetry (1);
			//	//sprintf (file_name, "Result//Melt_SF//speed_y_t_horizontal_dim1_%.5lf.txt", local_time_stamps[i]);
			//	//save_slice_first_derivative_weighted (2, 0, weights, "Source Files//Melt_SF//slice_Vy1.txt", file_name, true);
			//	//sprintf (file_name, "Result//Melt_SF//speed_y_t_horizontal_dim2_%.5lf.txt", local_time_stamps[i]);
			//	//save_slice_first_derivative_weighted (2, 0, weights, "Source Files//Melt_SF//slice_Vy2.txt", file_name, true);
			//	//sprintf (file_name, "Result//Melt_SF//speed_y_t_horizontal_dim3_%.5lf.txt", local_time_stamps[i]);
			//	//save_slice_first_derivative_weighted (2, 0, weights, "Source Files//Melt_SF//slice_Vy3.txt", file_name, true);
			//	//sprintf (file_name, "Result//Melt_SF//speed_y_t_horizontal_dim4_%.5lf.txt", local_time_stamps[i]);
			//	//save_slice_first_derivative_weighted (2, 0, weights, "Source Files//Melt_SF//slice_Vy4.txt", file_name, true);
			//	//sprintf (file_name, "Result//Melt_SF//speed_y_t_horizontal_dim5_%.5lf.txt", local_time_stamps[i]);
			//	//save_slice_first_derivative_weighted (2, 0, weights, "Source Files//Melt_SF//slice_Vy5.txt", file_name, true);
			//	//sprintf (file_name, "Result//Melt_SF//speed_y_t_horizontal_dim6_%.5lf.txt", local_time_stamps[i]);
			//	//save_slice_first_derivative_weighted (2, 0, weights, "Source Files//Melt_SF//slice_Vy6.txt", file_name, true);
			//	//sprintf (file_name, "Result//Melt_SF//speed_y_t_horizontal_dim7_%.5lf.txt", local_time_stamps[i]);
			//	//save_slice_first_derivative_weighted (2, 0, weights, "Source Files//Melt_SF//slice_Vy7.txt", file_name, true);
			//}

			//// gradient norm
			//weights[0] = 1.0;
			//weights[1] = -1.0;
			//set_symmetry (1);
			//sprintf (file_name, "Result//Melt_SF//grad_norm_%.5lf.txt", local_time_stamps[i]);
			//save_slice_first_derivative (2, 10, "Source Files//Melt_SF//bound_slice_bottom.txt", file_name, true);
			//
			//set_symmetry (1);
			//weights[0] = SCALE_SIZE;
			//weights[1] = -lambda_fluid * SCALE_TEMP / SCALE_SIZE;
			//sprintf (file_name, "Result//Melt_SF//grad_norm_dim_%.5lf.txt", local_time_stamps[i]);
			//save_slice_first_derivative_weighted (2, 10, weights, "Source Files//Melt_SF//bound_slice_bottom.txt", file_name, true);

			//// integral value
			//double x = 0.0;
			//double coordinates[2];
			//double value;
			//double sum_bottom = 0.0, sum_top = 0.0;
			//int counter_top = 0, counter_bottom = 0;
			//double coordinatesN[2];
			//mesh_pointer->get_N_boundaries (coordinatesN);

			//FILE * heat_file_bottom, *heat_file_bottom_dim, *heat_file_top, *heat_file_top_dim;
			//if (flag_append)
			//{
			//	heat_file_bottom = fopen ("Result//Melt_SF//integral_heat_flux_bottom.txt", "w");
			//	heat_file_bottom_dim = fopen ("Result//Melt_SF//integral_heat_flux_bottom_dim.txt", "w");
			//	heat_file_top = fopen ("Result//Melt_SF//integral_heat_flux_top.txt", "w");
			//	heat_file_top_dim = fopen ("Result//Melt_SF//integral_heat_flux_top_dim.txt", "w");
			//	flag_append = false;
			//}
			//else
			//{
			//	heat_file_bottom = fopen ("Result//Melt_SF//integral_heat_flux_bottom.txt", "a+");
			//	heat_file_bottom_dim = fopen ("Result//Melt_SF//integral_heat_flux_bottom_dim.txt", "a+");
			//	heat_file_top = fopen ("Result//Melt_SF//integral_heat_flux_top.txt", "a+");
			//	heat_file_top_dim = fopen ("Result//Melt_SF//integral_heat_flux_top_dim.txt", "a+");
			//}
			//for (int i = 0; x < coordinatesN[0] + 1e-7; i++)
			//{
			//	x = i * 0.001;
			//	coordinates[0] = x;

			//	coordinates[1] = 0.0;
			//	if (get_derivative (0, 1, coordinates, &value))
			//	{
			//		sum_bottom += value;
			//		counter_bottom++;
			//	}

			//	coordinates[1] = phase_change_og_position;
			//	if (get_derivative (0, 1, coordinates, &value))
			//	{
			//		sum_top += value;
			//		counter_top++;
			//	}
			//}

			//if (counter_bottom == 0)
			//	counter_bottom = 1;
			//sum_bottom = -sum_bottom * 2.0 * coordinatesN[0] / (double)counter_bottom;
			//fprintf (heat_file_bottom, "%.5lf\t%e\n", local_time_stamps[i], sum_bottom);
			//if (counter_top == 0)
			//	counter_top = 1;
			//sum_top = -sum_top * 2.0 * coordinatesN[0] / (double)counter_top;
			//fprintf (heat_file_top, "%.5lf\t%e\n", local_time_stamps[i], sum_top);

			//sum_bottom *= (L_source * 2.0 * SCALE_SIZE) * lambda_fluid * SCALE_TEMP / SCALE_SIZE;
			//fprintf (heat_file_bottom_dim, "%.5lf\t%e\n", local_time_stamps[i] * SCALE_TIME, sum_bottom);
			//
			//double c0[2], cN[2];
			//mesh_pointer->get_0_boundaries (c0);
			//mesh_pointer->get_N_boundaries (cN);
			//sum_top *= ((cN[0] - c0[0]) * 2.0 * SCALE_SIZE) * lambda_fluid * SCALE_TEMP / SCALE_SIZE;
			//fprintf (heat_file_top_dim, "%.5lf\t%e\n", local_time_stamps[i] * SCALE_TIME, sum_top);
			//fclose (heat_file_bottom);
			//fclose (heat_file_bottom_dim);			
			//fclose (heat_file_top);
			//fclose (heat_file_top_dim);

			// melt rate
			//melt_rate_counter++;
			//if (melt_rate_counter >= 10)
			//	melt_rate_counter = 0;
			//if (melt_rate_counter == 0)
			//{
			//	FILE * melt_rate_file = fopen ("Result//Melt_SF//melt_rate.txt", "a+");
			//	fprintf (melt_rate_file, "%.5lf\t%.7lf\n", local_time_stamps[i], get_melt_volume ());
			//	fclose (melt_rate_file);
			//}
			//// phase change position
			//sprintf (file_name, "Result//Melt_SF//phase_change_pos_%.5lf.txt", local_time_stamps[i]);
			//save_slice (0, phys_T_phase_change, 0, file_name, true);

			//// thermoconductive heat flux
			//double Qt, Qc;
			//get_Qs (0.25, &Qt, &Qc);
			//FILE * Qtc_file = fopen ("Result//Melt_SF//Qtc_025.txt", "a+");
			//fprintf (Qtc_file, "%.5lf\t%e\t%e\n", local_time_stamps[i], Qt, Qc);
			//fclose (Qtc_file);

			//get_Qs (0.5, &Qt, &Qc);
			//Qtc_file = fopen ("Result//Melt_SF//Qtc_05.txt", "a+");
			//fprintf (Qtc_file, "%.5lf\t%e\t%e\n", local_time_stamps[i], Qt, Qc);
			//fclose (Qtc_file);

			//get_Qs (0.75, &Qt, &Qc);
			//Qtc_file = fopen ("Result//Melt_SF//Qtc_075.txt", "a+");
			//fprintf (Qtc_file, "%.5lf\t%e\t%e\n", local_time_stamps[i], Qt, Qc);
			//fclose (Qtc_file);

			//get_Qs (0.95, &Qt, &Qc);
			//Qtc_file = fopen ("Result//Melt_SF//Qtc_095.txt", "a+");
			//fprintf (Qtc_file, "%.5lf\t%e\t%e\n", local_time_stamps[i], Qt, Qc);
			//fclose (Qtc_file);
		}
	}
}

template<class TypeMesh>
double Task_Melt_SF<TypeMesh>::get_melt_volume ()
{
	// two integrals
	double integral_fluid = 0.0;
	double integral_solid = 0.0;
	double * coordinates = new double[mesh_pointer->get_dimentionality ()];
	int counter_f;
	int counter_s;
	double value;

	// go by elements
	for (int i = 0, i_end = mesh_pointer->get_n_elements (); i < i_end; i++)
	{
		{
			// check how much of it is fluid
			// and check how much of it is solid
			counter_f = counter_s = 0;
			for (int k = 0, k_end = mesh_pointer->get_amount_of_def_nodes (i); k < k_end; k++)
			{
				mesh_pointer->get_node_coordinates (mesh_pointer->get_node_number (i, k), coordinates);
				get_solution_in_point (0, coordinates, &value);
				if (value > phys_T_phase_change + 1e-5)
					counter_f++;
				if (value < phys_T_phase_change - 1e-5)
					counter_s++;
			}
			// add the area to the fluid_integral
			integral_fluid += mesh_pointer->get_geometrical_area (i) * (double)counter_f / (double)mesh_pointer->get_amount_of_def_nodes (i);
			integral_solid += mesh_pointer->get_geometrical_area (i) * (double)counter_s / (double)mesh_pointer->get_amount_of_def_nodes (i);
		}
	}

	delete[] coordinates;
	// return difference between them 
	return integral_solid / integral_fluid;
}

template<class TypeMesh>
void Task_Melt_SF<TypeMesh>::get_Qs (double y, double * Qt, double * Qc)
{
	*Qt = *Qc = 0.0;
	double c0[2], cN[2];
	mesh_pointer->get_0_boundaries (c0);
	mesh_pointer->get_N_boundaries (cN);

	double h = (cN[0] - c0[0]) / (double)INTEGRAL_POINTS;
	double point[2];
	point[1] = y;
	double value;
	double der_value;
	for (int i = 0; i < INTEGRAL_POINTS + 1; i++)
	{
		point[0] = c0[0] + h * i;

		// Q = integrate -dT(x, y)/dy dx 
		if (get_derivative (0, 1, point, &der_value))
		{
			*Qt += -der_value;
		}

		if (get_solution_in_point (0, 0, point, &value))
		{
			// Q = integrate -dpsi(x,0.5)/dx T(x,0.5)dx
			if (get_derivative (2, 0, point, &der_value))
			{
				*Qc += -der_value * value;
			}
		}
	}
	*Qt *= h;
	*Qc *= h;
}

template<class TypeMesh>
Task_Melt_SF<TypeMesh>::Task_Melt_SF ()
{
	non_linear = true;

	FILE * file = fopen ("Source Files//Melt_SF//Parameters.txt", "r");
	fscanf (file, "%lf", &PR);
	fscanf (file, "%lf", &GR);
	fscanf (file, "%lf", &phase_change_og_position);
	fscanf (file, "%lf", &L_source);
	fscanf (file, "%lf", &phys_T1);
	fscanf (file, "%lf", &phys_T_phase_change);
	fscanf (file, "%lf", &del_T_smoothing);
	fscanf (file, "%lf", &Ste);
	fscanf (file, "%lf", &phys_T0_F);
	fscanf (file, "%lf", &phys_T0_S);
	fscanf (file, "%lf", &lambda_fluid);
	fscanf (file, "%lf", &lambda_solid);
	fclose (file);

	SCALE_TEMP = SCALE_TIME = SCALE_SIZE = SCALE_SPEED = 1.0;
	source_power = 0.0;
	flag_append = true;
	flag_calc_start = true;

	density_inversion = false;
	melt_rate_counter = 0;
}

template<class TypeMesh>
Task_Melt_SF<TypeMesh>::Task_Melt_SF (char * parameters_file_name)
{
	non_linear = true;
	
	FILE * file = fopen (parameters_file_name, "r");
	fscanf (file, "%lf", &PR);
	fscanf (file, "%lf", &GR);
	fscanf (file, "%lf", &phase_change_og_position);
	fscanf (file, "%lf", &L_source);
	fscanf (file, "%lf", &phys_T1);
	fscanf (file, "%lf", &phys_T_phase_change);
	fscanf (file, "%lf", &del_T_smoothing);
	fscanf (file, "%lf", &Ste);
	fscanf (file, "%lf", &phys_T0_F);
	fscanf (file, "%lf", &phys_T0_S);
	fscanf (file, "%lf", &lambda_fluid);
	fscanf (file, "%lf", &lambda_solid);
	fscanf (file, "%lf", &TEMP_T0);
	fclose (file);

	SCALE_TEMP = SCALE_TIME = SCALE_SIZE = SCALE_SPEED = 1.0;
	source_power = 0.0;
	flag_append = true;
	flag_calc_start = true;

	density_inversion = false;
	melt_rate_counter = 0;

	turn_on_source_time = -1.0;
	bottom_temp = 0.7500;
}

template<class TypeMesh>
Task_Melt_SF<TypeMesh>::Task_Melt_SF (const Task_Melt_SF & nst)
{
}

template<class TypeMesh>
Task_Melt_SF<TypeMesh>::~Task_Melt_SF ()
{
}
