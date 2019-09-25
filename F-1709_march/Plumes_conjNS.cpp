#include "Plumes_conjNs.h"

template class Plumes_Prototype <Triangular_Mesh>;
template class Plumes_Prototype <Triangular_Mesh_Hier>;
template class Plumes_Prototype <Mixed_Triangular_Mesh>;

template <class Type_Mesh>
void Plumes_Prototype<Type_Mesh>::set_n_systems ()
{
	n_systems = 3;
	eq_matrixes = new compressed_matrix[n_systems];
	printf ("systems settled:\t%i\n", n_systems);
	fprintf (log_file, "systems settled:\t%i\n", n_systems);
}

template <class Type_Mesh>
void Plumes_Prototype<Type_Mesh>::build_system (int k_system)
{
	switch (k_system)
	{
	case 0:
	{
		Matrix * G;
		Matrix * M;
		Matrix * DDFx;
		Matrix * DDFy;
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

				n_functions = n_functions_cur;

				G = new Matrix (n_functions, n_functions);
				M = new Matrix (n_functions, n_functions);
				Psi = new MathVector (n_functions);
				DDFx = new Matrix (n_functions, n_functions);
				DDFy = new Matrix (n_functions, n_functions);
				Q = new MathVector (n_functions);
				MQ = new MathVector (n_functions);
			}

			// get G, M matrices 
			mesh_pointer->get_G_local_matrix (k_element, G);
			mesh_pointer->get_M_local_matrix (k_element, M);
			
			area = mesh_pointer->get_area (k_element);
			if (area == 1) // fluid
			{
				double coef = 1.0 / (PR * sqrt (GR));
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
								eq_matrixes[k_system].add_to_entry (iF, jF, G->Elem (i, j) * coef);
								// add M (by t) into respective places (by def_nodes)
								eq_matrixes[k_system].add_to_entry (iF, jF, M->Elem (i, j) * Time_coef->getElem (0));
							}
						}
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
			else
			{
				double coef = 1.0 / PR;
				coef *= therm_solid / therm_fluid;

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
								eq_matrixes[k_system].add_to_entry (iF, jF, G->Elem (i, j) * coef);
								// add M (by t) into respective places (by def_nodes)
								eq_matrixes[k_system].add_to_entry (iF, jF, M->Elem (i, j) * Time_coef->getElem (0));
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

		}

		delete Time_coef;
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

		double coef;
		coef = 1.0 / sqrt (GR);

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
			}

			area = mesh_pointer->get_area (k_element);
			if (area == 1) // fluid
			{
				// get G, M matrices 
				mesh_pointer->get_G_local_matrix (k_element, G);
				mesh_pointer->get_M_local_matrix (k_element, M);

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
				D->MultiplyMatrixByVector (*T, TD);

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
								eq_matrixes[k_system].add_to_entry (iF, jF, G->Elem (i, j) * coef);
								// add M (by t) into respective places (by def_nodes)
								eq_matrixes[k_system].add_to_entry (iF, jF, M->Elem (i, j) * Time_coef->getElem (0));
							}
						}

						// add dt/dx
						eq_matrixes[k_system].add_to_f_entry (iF, TD->getElem (i));
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
		}

		for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
		{
			n_functions_cur = mesh_pointer->get_amount_non_zero_functions (k_element);
			if (n_functions != n_functions_cur)
			{			
				n_functions = n_functions_cur;
			}

			area = mesh_pointer->get_area (k_element);
			if (area == 2) // solid
			{
				// "first condition" 0
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
		delete Time_coef;
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
		break;
	}
	case 2:
	{
		Matrix * G;
		Matrix * M;
		MathVector * W;
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

		for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
		{
			n_functions_cur = mesh_pointer->get_amount_non_zero_functions (k_element);
			if (n_functions != n_functions_cur)
			{
				delete G;
				delete M;
				delete W;
				delete WM;

				n_functions = n_functions_cur;

				G = new Matrix (n_functions, n_functions);
				M = new Matrix (n_functions, n_functions);
				W = new MathVector (n_functions);
				WM = new MathVector (n_functions);
			}
			
			area = mesh_pointer->get_area (k_element);
			if (area == 1) // fluid
			{
				// get G, M matrices 
				mesh_pointer->get_G_local_matrix (k_element, G);
				mesh_pointer->get_M_local_matrix (k_element, M);

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

						// put W into respective places (by def_nodes)
						eq_matrixes[k_system].add_to_f_entry (iF, WM->getElem (i));
					}
				}
			}
		}

		for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
		{
			n_functions_cur = mesh_pointer->get_amount_non_zero_functions (k_element);
			if (n_functions != n_functions_cur)
			{
				n_functions = n_functions_cur;
			}

			area = mesh_pointer->get_area (k_element);
			if (area == 2) // solid
			{
				// "first condition" 0
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
		delete G;
		delete M;
		delete W;
		delete WM;
		break;

	}
	}
}

template <class Type_Mesh>
void Plumes_Prototype<Type_Mesh>::apply_first_boundary_conditions (int k_system)
{
	Task::apply_first_boundary_conditions (k_system);

	// for T add temperature on the bottom
	if (k_system == 0)
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
					if (fabs (coordinates[0]) < L_source + ZERO_Plumes)
					{
						mesh_pointer->get_node_coordinates (n2, coordinates);
						if (fabs (coordinates[0]) < L_source + ZERO_Plumes)
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
							get_FCondition_edge (k_system, k_element, n1, n2, functions, boundary, FCondition);

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
	}
}

template <class Type_Mesh>
void Plumes_Prototype<Type_Mesh>::print_solutions ()
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

		char name[64];
		for (int k = 0; k < n_systems; k++)
		{
			for (size_t i = 0, i_end = local_time_stamps.size (); i < i_end; i++)
			{
				get_time_approx (local_time_stamps[i], c);

				printf ("\tsaving solutions into files, time layer:\t%.16lf\n", local_time_stamps[i]);
				sprintf (name, "Result//By time layers//s%i_t_%0.5lf.txt", k, local_time_stamps[i]);
				file = fopen (name, "w");
				double val;
				for (int j = 0, j_end = previous_time_layers_solutions[k][0].getSize (); j < j_end; j++)
				{
					val = 0;
					//for (int t = 0; t < time_sampling; t++)
					//{
					//	val += c->getElem (t) * previous_time_layers_solutions[k][t].getElem (j);
					//}
					val = previous_time_layers_solutions[k][0].getElem (j);
					fprintf (file, "%.16lf\n", val);
				}
				fclose (file);

				if (false)
				{
					//// add divider into painter
					//std::vector <std::pair<int, int>> divider;
					//mesh.get_area_dividers (&divider);
					//double coordinates_1[2];
					//double coordinates_2[2];
					//for (size_t i = 0, i_end = divider.size (); i < i_end; i++)
					//{
					//	mesh_pointer->get_node_coordinates (divider[i].first, coordinates_1);
					//	mesh_pointer->get_node_coordinates (divider[i].second, coordinates_2);
					//	painter->addLine (coordinates_1[0], coordinates_1[1], coordinates_2[0], coordinates_2[1]);
					//}

					wchar_t pic_name[128];
					int n_materials = 1;
					if (k == 0)
					{
						//painter->set_min_max (phys_T0, phys_T1);
						n_materials = 2;
					}

					for (int i_material = 0; i_material < n_materials; i_material++)
					{
						//painter->draw_subsection (i_material + 1);
						{
							swprintf (pic_name, 128, L"Pictures//s%i%i_t_%lf.png", k, i_material, local_time_stamps[i]);
							std::wstring pic_name_field (pic_name);
							//painter->draw_isovalues (k, pic_name_field);
						}
						{
							swprintf (pic_name, 128, L"Pictures//s%i%i_t_%lf_iso.png", k, i_material, local_time_stamps[i]);
							std::wstring pic_name_iso (pic_name);
							//painter->draw_isovalues (k, pic_name_iso);
						}
						{
							swprintf (pic_name, 128, L"Pictures//s%i%i_t_%lf_field.png", k, i_material, local_time_stamps[i]);
							std::wstring pic_name_field (pic_name);
							//painter->draw_field (k, pic_name_field);
						}
					}
				}
			}


		}
	}
}

template <class Type_Mesh>
double Plumes_Prototype<Type_Mesh>::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	// theoretically zero flow should be applied automatically through Theta in get_Theta_edge
	double r = 0.0;

	switch (k_system)
	{
	case 0:
	{
		switch (boundary)
		{
		case 2:
			r = phys_T1;
			break;
		case 3:
			r = phys_T0;
			break;
		}
		break;
	}
	case 1:
	{
		double point[2];
		point[0] = coordinates[0];
		point[1] = coordinates[1];

		if (slip == 1)
		{
			r = 0.0;
		}
		else
		{
			switch (boundary)
			{
			case 0:
				point[0] += P_step;
				break;
			case 1:
				point[0] -= P_step;
				break;
			case 2:
				point[1] += P_step;
				break;
			case 3:
				point[1] -= P_step;
				break;
			}
			break;
			get_solution_in_point (2, point, &r);
			r *= -2.0 / pow (P_step, 2.0);
		}
		break;
	}
	case 2:
	{
		r = 0.0;
		break;
	}
	}
	return r;
}

template<class TypeMesh>
double Plumes_Prototype<TypeMesh>::function_starting_condition (int k_system, double * coordinates, int area)
{
	double r = 0.0;
	switch (k_system)
	{
	case 0:
	{
		r = phys_T0;
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
void Plumes_Prototype<TypeMesh>::solution_step (int k_system, int method, int d_type, int depth)
{
	double discr = 1.0; // discrepancy

						// implicit simple iteration
	int solver_iterations;
	// previous non-linear solution
	MathVector * prev_solution = new MathVector (eq_matrixes[k_system].Size ());

	eq_matrixes[k_system].Clear ();
	// set last time layer's solution as starting point 
	eq_matrixes[k_system].set_starting_point (&(previous_time_layers_solutions[k_system][0]));
	build_system (k_system);  // reset system for current layer
	apply_boundary_conditions (k_system); //apply boundary conditions

	// solve
	switch (method)
	{
	case 0:
		solver_iterations = eq_matrixes[k_system].solve_GMRES (d_type, depth); // solve it	
		break;
	case 1:
		solver_iterations = eq_matrixes[k_system].solve_LOS (d_type); // solve it	
		break;
	default:
		solver_iterations = eq_matrixes[k_system].solve_LOS (d_type); // solve it	
		break;
	}

	// get 0-step solution
	eq_matrixes[k_system].get_solution (prev_solution);

	printf ("system: %i, time layer:\t%lf\n", k_system, time_layers[current_time_layer]);
	fprintf (log_file, "system: %i, time layer: %lf\n", k_system, time_layers[current_time_layer]);

	int i = 0;
	for (i; i < IMPLICIT_SIMPLE_ITERATION_MAX_ITER && discr > PL_DISCRPANSCY; i++)
	{
		eq_matrixes[k_system].set_starting_point (prev_solution); // set last non-linear step solution as x0 for solver

		switch (method)
		{
		case 0:
			solver_iterations = eq_matrixes[k_system].solve_GMRES (d_type, depth); // solve it	
			break;
		case 1:
			solver_iterations = eq_matrixes[k_system].solve_LOS (d_type); // solve it	
			break;
		default:
			solver_iterations = eq_matrixes[k_system].solve_LOS (d_type); // solve it	
			break;
		}
		fprintf (log_file, "non_linear iter: %i, solver: %i", i, solver_iterations);

		// relax with prev_solution and x
		discr = relax (0, prev_solution);

		// save current solution as prev one
		eq_matrixes[k_system].get_solution (prev_solution);
	}

	printf ("\tnon-linear iterations: %i, discr: %e\n", i, discr);
	next_time_layer (k_system); // save solution, move up time layer

	delete prev_solution;
}

template <class Type_Mesh>
Plumes_Prototype<Type_Mesh>::Plumes_Prototype ()
{
	FILE * file = fopen ("Source Files//Plumes//Parameters.txt", "r");
	fscanf (file, "%lf", &PR);
	fscanf (file, "%lf", &GR);
	fscanf (file, "%lf", &L_source);
	fscanf (file, "%lf", &phys_T0);
	fscanf (file, "%lf", &phys_T1);
	fscanf (file, "%lf", &P_step);
	fscanf (file, "%i", &slip);	
	fscanf (file, "%lf", &therm_fluid);
	fscanf (file, "%lf", &therm_solid);
	fclose (file);
}

template <class Type_Mesh>
Plumes_Prototype<Type_Mesh>::Plumes_Prototype (const Plumes_Prototype<Type_Mesh> & nst)
{
}

template <class Type_Mesh>
Plumes_Prototype<Type_Mesh>::~Plumes_Prototype ()
{
}
