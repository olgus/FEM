#include "Time_Dependent_Task.h"

template class Time_Dependent_Task <Mesh<Node<Point_Prototype>, Element>>;
template class Time_Dependent_Task <Mesh<Node<Point_2D>, Element>>;
template class Time_Dependent_Task <Mesh<Node<Point_3D>, Element>>;
template class Time_Dependent_Task <Mesh<Node<Point_2D_Polar>, Element>>;
template class Time_Dependent_Task <Mesh<Node_2D, Element>>;
template class Time_Dependent_Task <Mesh<Node_2D_Polar, Element>>;
template class Time_Dependent_Task <Mesh<Node_1D, Element>>;
template class Time_Dependent_Task <Mesh<Node_3D, Element>>;

template class Time_Dependent_Task <Mesh_1D_L1>;
template class Time_Dependent_Task <Mesh_1D_Hier>;
template class Time_Dependent_Task <Mesh_1D_Hermitian>;
template class Time_Dependent_Task <Triangular_Mesh>;
template class Time_Dependent_Task <Mixed_Triangular_Mesh>;
template class Time_Dependent_Task <Triangular_Polar_Mesh>;
template class Time_Dependent_Task <Triangular_Mesh_Hier>;
template class Time_Dependent_Task <Rectangular_Mesh>;
template class Time_Dependent_Task <Rectangular_Mesh_3>;
template class Time_Dependent_Task <Rectangular_Mesh_Spline>;
template class Time_Dependent_Task <Prismatic_Mesh>;
template class Time_Dependent_Task <Cubic_Mesh>;

template<class Type_Mesh>
void Time_Dependent_Task<Type_Mesh>::map_time (char * file_name)
{
	FILE * file = fopen (file_name, "r");

	int lvl; // nest lvl
	fscanf (file, "%i", &lvl);
	int n;
	fscanf (file, "%i", &n);

	// get sections' boundaries values
	double * boundaries = new double[n + 1];
	for (int i = 0; i < n + 1; i++)
	{
		fscanf (file, "%lf", &boundaries[i]);
	}

	// get amount of subsections
	int * n_sections = new int[n];
	for (int i = 0; i < n; i++)
	{
		fscanf (file, "%i", &n_sections[i]);
		n_sections[i] *= (int)round (pow (2.0, lvl));
	}

	// get coefficients of subsections
	double * coef = new double[n];
	for (int i = 0; i < n; i++)
	{
		fscanf (file, "%lf", &coef[i]);
	}

	// get directions of subsections
	int direc;
	for (int i = 0; i < n; i++)
	{
		fscanf (file, "%i", &direc);
		if (direc == -1)
			coef[i] = 1.0 / coef[i];
	}
	fclose (file);

	int counter = 1;
	double l0; // first lenght
	double lcur; // length of current segment
	double GPS; // geometric sum
	double L; // length of a section
	double val;

	val = boundaries[0];
	time_stamps.push_back (val);
	for (int i = 0; i < n; i++)
	{
		// get GPS
		if (fabs (coef[i] - 1.0) < ZERO_mesh_1D)
		{
			GPS = n_sections[i];
		}
		else
		{
			GPS = (1.0 - pow (coef[i], n_sections[i])) / (1.0 - coef[i]); // TEST
		}
		// get L
		L = (boundaries[i + 1] - boundaries[i]);
		// get l0
		l0 = L / GPS;
		for (int j = 0; j < n_sections[i]; j++)
		{
			lcur = l0 * pow (coef[i], j);
			val += lcur;
			time_stamps.push_back (val);
		}
	}
	delete[] coef;
	delete[] n_sections;
	delete[] boundaries;
}

template<class Type_Mesh>
Time_Dependent_Task<Type_Mesh>::Time_Dependent_Task ()
{
	start_cond = true;

	time_sampling = 2;
}

template<class Type_Mesh>
Time_Dependent_Task<Type_Mesh>::Time_Dependent_Task (const Time_Dependent_Task & tdt)
{
}

template<class Type_Mesh>
Time_Dependent_Task<Type_Mesh>::~Time_Dependent_Task ()
{
}

template<class Type_Mesh>
void Time_Dependent_Task<Type_Mesh>::set_scheme_order (int Scheme_order)
{
	time_sampling = Scheme_order + 1;
}

template<class Type_Mesh>
void Time_Dependent_Task<Type_Mesh>::get_time_coefficients (MathVector * c)
{
	int scheme = current_time_layer < (time_sampling - 1) ? current_time_layer : (time_sampling - 1);

	c->Zero ();

	switch (scheme)
	{
	case 1:
	{
		double dt = time_layers[current_time_layer] - time_layers[current_time_layer - 1];
		c->setElem (0, 1.0 / dt);
		c->setElem (1, 1.0 / dt);
		break;
	}
	case 2:
	{
		double t0 = time_layers[current_time_layer];
		double t1 = time_layers[current_time_layer - 1];
		double t2 = time_layers[current_time_layer - 2];

		double dt = t0 - t2;
		double dt1 = t1 - t2;
		double dt0 = t0 - t1;
		
		c->setElem (0, (dt + dt0) / (dt * dt0));
		c->setElem (1, dt / (dt0 * dt1));
		c->setElem (2, -dt0 / (dt * dt1));		
		break;
	}
	case 3:
	{
		double t0 = time_layers[current_time_layer];
		double t1 = time_layers[current_time_layer - 1];
		double t2 = time_layers[current_time_layer - 2];
		double t3 = time_layers[current_time_layer - 3];

		double dt1 = t0 - t1;
		double dt2 = t0 - t2;
		double dt3 = t0 - t3;
		double dt12 = t1 - t2;
		double dt13 = t1 - t3;
		double dt23 = t2 - t3;

		c->setElem (0, 1.0 / dt3 + 1.0 / dt2 + 1.0 / dt1);
		c->setElem (1, (dt2 * dt3) / (dt13 * dt12 * dt1));
		c->setElem (2, -(dt3 * dt1) / (dt23 * dt12 * dt2));
		c->setElem (3, (dt2 * dt1) / (dt23 * dt13 * dt3));
		break;
	}
	default:
	{
		break;
	}
	}
}

template<class Type_Mesh>
void Time_Dependent_Task<Type_Mesh>::get_time_approx (double t, MathVector * c)
{
	int scheme = current_time_layer < (time_sampling - 1) ? current_time_layer : (time_sampling - 1);

	c->Zero ();

	switch (scheme)
	{
	case 1:
	{
		double dt = time_layers[current_time_layer] - time_layers[current_time_layer - 1];
		c->setElem (0, (t - time_layers[current_time_layer - 1]) / dt);
		c->setElem (1, (time_layers[current_time_layer] - t) / dt);
		break;
	}
	case 2:
	{
		double t0 = time_layers[current_time_layer];
		double t1 = time_layers[current_time_layer - 1];
		double t2 = time_layers[current_time_layer - 2];

		double dt = t0 - t2;
		double dt1 = t1 - t2;
		double dt0 = t0 - t1;

		c->setElem (0, (t - t1) * (t - t2) / (dt * dt0));
		c->setElem (1, -(t - t2) * (t - t0) / (dt0 * dt1));
		c->setElem (2, (t - t0) * (t - t1) / (dt * dt1));
		break;
	}
	case 3:
	{
		double t0 = time_layers[current_time_layer];
		double t1 = time_layers[current_time_layer - 1];
		double t2 = time_layers[current_time_layer - 2];
		double t3 = time_layers[current_time_layer - 3];

		double dt1 = t0 - t1;
		double dt2 = t0 - t2;
		double dt3 = t0 - t3;
		double dt12 = t1 - t2;
		double dt13 = t1 - t3;
		double dt23 = t2 - t3;

		c->setElem (0, (t - t3) * (t - t2) * (t - t1) / (dt1 * dt2 * dt3));
		c->setElem (1, -(t - t3) * (t - t2) * (t - t0) / (dt13 * dt12 * dt1));
		c->setElem (2, (t - t3) * (t - t1) * (t - t0) / (dt23 * dt12 * dt2));
		c->setElem (3, -(t - t2) * (t - t1) * (t - t0) / (dt23 * dt13 * dt3));
		break;
	}
	default:
	{
		break;
	}
	}
}

template<class Type_Mesh>
void Time_Dependent_Task<Type_Mesh>::set_time_layers (char * file_name)
{
	FILE * file = fopen (file_name, "r");
	
	fscanf (file, "%i", &time_sampling);
	time_sampling++;

	previous_time_layers_solutions = new MathVector *[n_systems]; // solutions are saved for each system
	for (int k_systems = 0; k_systems < n_systems; k_systems++)
	{
		previous_time_layers_solutions[k_systems] = new MathVector[time_sampling];
	}

	int lvl; // nest lvl
	fscanf (file, "%i", &lvl);
	int n;
	fscanf (file, "%i", &n);

	// get sections' boundaries values
	double * boundaries = new double[n + 1];
	for (int i = 0; i < n + 1; i++)
	{
		fscanf (file, "%lf", &boundaries[i]);
	}

	n_time_layers = 0;
	// get amount of subsections
	int * n_sections = new int[n];
	for (int i = 0; i < n; i++)
	{
		fscanf (file, "%i", &n_sections[i]);
		n_sections[i] *= (int)round (pow (2.0, lvl));
		n_time_layers += n_sections[i];
	}

	// get coefficients of subsections
	double * coef = new double[n];
	for (int i = 0; i < n; i++)
	{
		fscanf (file, "%lf", &coef[i]);
	}

	// get directions of subsections
	int direc;
	for (int i = 0; i < n; i++)
	{
		fscanf (file, "%i", &direc);
		if (direc == -1)
			coef[i] = 1.0 / coef[i];
	}
	fclose (file);

	n_time_layers++;
	// map layers according to the "net"
	time_layers = new double[n_time_layers];

	time_layers[0] = boundaries[0]; 
	int counter = 1;
	double l0; // first lenght
	double lcur; // length of current segment
	double GPS; // geometric sum
	double L; // length of a section

	for (int i = 0; i < n; i++)
	{
		// get GPS
		if (fabs (coef[i] - 1.0) < ZERO_mesh_1D)
		{
			GPS = n_sections[i];
		}
		else
		{
			GPS = (1.0 - pow (coef[i], n_sections[i])) / (1.0 - coef[i]); // TEST
		}
		// get L
		L = (boundaries[i + 1] - boundaries[i]);
		// get l0
		l0 = L / GPS;

		for (int j = 0; j < n_sections[i]; j++)
		{
			lcur = l0 * pow (coef[i], j);
			time_layers[counter] = time_layers[counter - 1] + lcur;
			counter++;
		}
	}
	delete[] coef;
	delete[] n_sections;
	delete[] boundaries;

	printf ("time layers:\t%i\ntime sampling:\t%i \n", n_time_layers, time_sampling);
	fprintf (log_file, "time layers:\t%i\ntime sampling:\t%i \n", n_time_layers, time_sampling);
}

template<class Type_Mesh>
void Time_Dependent_Task<Type_Mesh>::build_system (int k_system)
{
	Matrix * G;
	Matrix * M;
	MathVector * B;
	MathVector * Q;
	MathVector * MQ;
	MathVector * Time_coef;
	
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
		area = mesh. elements[k_element]->get_area ();
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
		delete B;
		delete G;
		delete M;
		delete Q;
		delete MQ;
	}

	delete Time_coef;
}

template<class Type_Mesh>
double Time_Dependent_Task<Type_Mesh>::function_f (int k_system, double * coordinates, int area)
{
	double t = time_layers[current_time_layer];
	double x = coordinates[0];
	double y = coordinates[1];

	//return 1.0;
	//return 2.0 * t + 1.0;
	return 1.0 + 2.0 * t + 3.0 * t * t;
	//return -sin (t);
}

template<class Type_Mesh>
double Time_Dependent_Task<Type_Mesh>::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	double t = time_layers[current_time_layer];
	double x = coordinates[0];
	double y = coordinates[1];
	
	//return t;
	//return  t * t + t;
	return  t * t * t + t * t + t + x + y;
	//return cos (t);
}

template<class Type_Mesh>
double Time_Dependent_Task<Type_Mesh>::function_starting_condition (int k_system, double * coordinates, int area)
{
	double t = time_layers[current_time_layer];
	double x = coordinates[0];
	double y = coordinates[1];

	//return t;
	return  t * t * t + t * t + t + x + y;
}

template<class Type_Mesh>
void Time_Dependent_Task<Type_Mesh>::set_starting_conditions ()
{
	for (int k_system = 0; k_system < n_systems; k_system++)
	{
		for (int i = 0; i < time_sampling; i++)
		{
			previous_time_layers_solutions[k_system][i].setSize (eq_matrixes[k_system].Size ()); // default starting solution point is 0
		}
	}

	Matrix * M;
	MathVector * B;
	int n_functions;
	int n_functions_cur;
	int iF, jF;
	int dim = mesh_pointer->get_dimentionality ();

	n_functions = mesh_pointer->get_amount_non_zero_functions (0);
	M = new Matrix (n_functions, n_functions);
	B = new MathVector (n_functions);

	for (int k_system = 0; k_system < n_systems; k_system++)
	{
		for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
		{
			n_functions_cur = mesh_pointer->get_amount_non_zero_functions (k_element);
			if (n_functions != n_functions_cur)
			{
				delete B;
				delete M;

				n_functions = n_functions_cur;

				M = new Matrix (n_functions, n_functions);
				B = new MathVector (n_functions);
			}

			// get M matrix 
			mesh_pointer->get_M_local_matrix (k_element, M);

			// get B-vector
			get_local_SCB (k_system, k_element, B);

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
							eq_matrixes[k_system].add_to_entry (iF, jF, M->Elem (i, j));
						}
					}

					// put B into respective places (by def_nodes)
					eq_matrixes[k_system].add_to_f_entry (iF, B->getElem (i));
				}
			}
		}
	
		eq_matrixes[k_system].solve_LOS (1);
		eq_matrixes[k_system].get_solution (&(previous_time_layers_solutions[k_system][0]));
		//char name[128];
		//sprintf (name, "Result//Melt_SF//test_s%i.txt", k_system);
		//FILE * file = fopen (name, "w");
		//previous_time_layers_solutions[k_system][0].FPrint (file);
		//fclose (file);
	}

	//double c[2];
	//for (int k_system = 0; k_system < n_systems; k_system++)
	//{
	//	for (int k = 0; k < mesh_pointer->get_n_nodes (); k++)
	//	{
	//		mesh_pointer->get_node_coordinates (k, c);
	//		previous_time_layers_solutions[k_system][0].setElem (k, function_starting_condition (k_system, c, 1));
	//	}
	//}
	printf ("\tstarting conditions set\n");
	delete B;
	delete M;
}

template<class Type_Mesh>
void Time_Dependent_Task<Type_Mesh>::get_local_SCB (int k_system, int k_element, MathVector * B)
{
	int n_funcs = B->getSize ();
	int area = mesh_pointer->get_area (k_element);

	double jac = 0.0; // jacobian
	int dim = mesh_pointer->get_dimentionality ();
	int n_integr_points = mesh.elements[k_element]->amount_of_integration_points ();
	double * weigths = new double[n_integr_points];
	double ** points = new double *[n_integr_points];
	for (int i = 0; i < n_integr_points; i++)
	{
		points[i] = new double[dim];
	}

	// get integration points for the element
	mesh.elements[k_element]->integration_points (points, weigths, &jac);
	
	double val;
	// go by element's functions
	for (int i = 0; i < n_funcs; i++)
	{
		val = 0.0;
		// integrate them * basis_functions on the element
		for (int j = 0; j < n_integr_points; j++)
		{
			// sum them multiplying by weigths
			val += mesh.elements[k_element]->get_basis_function_value (i, points[j]) * function_starting_condition (k_system, points[j], area) * weigths[j];
		}
		// multiply by jac
		val *= jac;
		B->setElem (i, val);
	}
}

template<class Type_Mesh>
void Time_Dependent_Task<Type_Mesh>::prepare (char * mesh_file_name, char * bound_file_name, char * time_file_name, char * stamps_file_name, bool save)
{
	Task<Type_Mesh>::prepare (mesh_file_name, bound_file_name, time_file_name, save);

	use_time_mapping = false;

	// map time
	map_time (stamps_file_name);
}

template<class Type_Mesh>
void Time_Dependent_Task<Type_Mesh>::prepare (char * mesh_file_name, char * bound_file_name[], char * time_file_name, char * stamps_file_name, bool save)
{
	Task<Type_Mesh>::prepare (mesh_file_name, bound_file_name, time_file_name, save);

	use_time_mapping = false;

	// map time
	map_time (stamps_file_name);
}

template<class Type_Mesh>
void Time_Dependent_Task<Type_Mesh>::prepare (char * nodes_file_name, char * elements_file_name, char * bound_file_name, char * time_file_name, char * stamps_file_name)
{
	Task<Type_Mesh>::prepare (nodes_file_name, elements_file_name, bound_file_name, time_file_name);

	use_time_mapping = false;

	// map time
	map_time (stamps_file_name);
}

template<class Type_Mesh>
void Time_Dependent_Task<Type_Mesh>::prepare (char * nodes_file_name, char * elements_file_name, char * bound_file_name[], char * time_file_name, char * stamps_file_name)
{
	Task<Type_Mesh>::prepare (nodes_file_name, elements_file_name, bound_file_name, time_file_name);

	use_time_mapping = false;

	// map time
	map_time (stamps_file_name);
}

template<class Type_Mesh>
bool Time_Dependent_Task<Type_Mesh>::start_from_selected_time_layer (double time_value, char * solutions[])
{
	start_cond = false;

	// find time_layer in time_layers 
	int pos = n_time_layers;

	int first = 0, last = n_time_layers - 1, mid;
	while (first < last)
	{
		mid = (first + last) / 2;
		if (time_value < time_layers[mid] + 1e-7) // <=
			last = mid;
		else
			first = mid + 1;
	}
	if (last < n_time_layers && time_layers[last])
	{
		pos = last;
	}

	if (pos == n_time_layers)
		return false;
	else
	{
		{
			int new_n_time_layers = n_time_layers - pos;
			double * new_time_layers = new double[new_n_time_layers];
			for (int i = 0; i < new_n_time_layers; i++)
			{
				new_time_layers[i] = time_layers[i + pos];
			}
			delete[] time_layers;
			time_layers = new double[new_n_time_layers];
			for (int i = 0; i < new_n_time_layers; i++)
			{
				time_layers[i] = new_time_layers[i];
			}
			n_time_layers = new_n_time_layers;

			delete[] new_time_layers;
		}

		std::vector < double> new_time_stamps;
		for (int i = 0; i < (int)time_stamps.size (); i++)
		{
			if (time_stamps[i] > time_value)
				new_time_stamps.push_back (time_stamps[i]);
		}
		time_stamps.clear ();
		time_stamps.insert (time_stamps.begin (), new_time_stamps.begin (), new_time_stamps.end ());

		// read solutions into k_time_layer = 0
		for (int i = 0; i < n_systems; i++)
		{
			read_solution (solutions[i], i, 0);
		}
		return true;
	}
}
