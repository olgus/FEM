#include "Spline_1D_Smooth.h"

void Spline_1D_Smooth::get_local_A_matrix (int k_element, Matrix * A_matrix)
{
	// set A_matrix to zero
	A_matrix->Zero ();
	// allocate some memory
	int matr_dim = mesh.elements[k_element]->get_amount_non_zero_functions ();
	Matrix * A = new Matrix (matr_dim, matr_dim);

	int dim = mesh.get_dimentionality ();
	double * coordinates = new double[dim];

	int m[] = { 1 };
	// go through solution points
	for (int k = 0; k < n_sol_points; k++)
	{
		for (int j = 0; j < dim; j++)
		{
			coordinates[j] = solution_points->getElem (k);
		}

		// check if it's in the current element 
		if (mesh.elements[k_element]->point_inside (mesh, coordinates))
		{
			// if it is, get matrix A[ij] = Psi[i] (solution_points[k]) * Psi[j] (solution_points[k])
			for (int i = 0; i < matr_dim; i++)
			{
				for (int j = 0; j < matr_dim; j++)
				{
					A->setElem (i, j, mesh.elements[k_element]->get_function_value (i, j, m, coordinates));
				}
			}
			// multiply by weights[k]
			// A->MultiplyByNumber (weigths->getElem (k));
			// add it to A_matrix
			A_matrix->Add (*A);
		}
	}

	delete A;
	delete[] coordinates;
}

void Spline_1D_Smooth::get_local_B (int k_element, MathVector * B_vector)
{
	// set B_vector to zero
	B_vector->Zero ();
	// allocate some memory
	int matr_dim = mesh.elements[k_element]->get_amount_non_zero_functions ();
	MathVector * B = new MathVector (matr_dim);

	int dim = mesh.get_dimentionality ();
	double * coordinates = new double[dim];
	double val;
	// go through solution points
	for (int k = 0; k < n_sol_points; k++)
	{
		for (int j = 0; j < dim; j++)
		{
			coordinates[j] = solution_points->getElem (k);
		}

		// check if it's in the current element 
		if (mesh.elements[k_element]->point_inside (mesh, coordinates))
		{
			// if it is, get Psi[i] (solution_points[k]) 
			for (int i = 0; i < matr_dim; i++)
			{
				val = mesh.elements[k_element]->get_basis_function_value (i, coordinates);
				B->setElem (i, val);
			}
			// multiply it by solution_values[k]
			B->MultiplyByNumber (solution_values->getElem (k));
			// multiply it by weights[k]
			//B->MultiplyByNumber (weigths->getElem (k));
			// add it to B_vector
			B_vector->Add (*B);
		}
	}

	delete B;
	delete[] coordinates;
}

double Spline_1D_Smooth::function_f (int k_system, double * coordinates, int area)
{
	return 0.0;
}

void Spline_1D_Smooth::build_system (int k_system)
{
	Matrix * A;
	Matrix * G;
	MathVector * B;
	int * def_nodes;
	int n_entries;
	int n_def_nodes;
	int iB, jB;

	// code for test purposes
	FILE * fileG = fopen ("Result Files Extra//G.txt", "w");

	for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
	{
		n_entries = mesh.elements[k_element]->get_amount_non_zero_functions ();
		n_def_nodes = mesh.elements[k_element]->get_amount_of_def_nodes ();
		// get A, G, D matrices and B vector
		A = new Matrix (n_entries, n_entries);
		G = new Matrix (n_entries, n_entries);
		B = new MathVector (n_entries);

		get_local_A_matrix (k_element, A);
		mesh.elements[k_element]->get_G_local_matrix (G);
		get_local_B (k_element, B);

		if (k_element == 0)
		{
			G->FPrint (fileG);
		}

		def_nodes = new int[n_def_nodes];
		mesh.elements[k_element]->get_def_nodes (def_nodes);
		// go by element's functions
		for (int i = 0; i < n_entries; i++)
		{
			iB = get_function_global_number (k_element, i);
			for (int j = 0; j < n_entries; j++)
			{
				jB = get_function_global_number (k_element, j);

				eq_matrixes[k_system].add_to_entry (iB, jB, A->Elem (i, j));
				eq_matrixes[k_system].add_to_entry (iB, jB, alpha * G->Elem (i, j));
			}
			eq_matrixes[k_system].add_to_f_entry (iB, B->getElem (i));
		}
		delete[] def_nodes;
		delete A;
		delete G;
		delete B;
	}

	fclose (fileG);
}

void Spline_1D_Smooth::apply_boundary_conditions (int k_system)
{
	// no conditions
}

Spline_1D_Smooth::Spline_1D_Smooth ()
{
	solution_values = NULL;
	solution_points = NULL;
}

Spline_1D_Smooth::Spline_1D_Smooth (const Spline_1D_Smooth & spline)
{
}

Spline_1D_Smooth::~Spline_1D_Smooth ()
{
	if (solution_values != NULL)
		delete solution_values;	
	if (solution_points != NULL)
		delete solution_points;
}

void Spline_1D_Smooth::prepare (char * file_name, double * extra)
{
	double x, y;
	std::vector <std::pair<double, double>> points_data;

	FILE * f = fopen (file_name, "r");
	while (!feof (f))
	{
		fscanf (f, "%lf %lf", &x, &y);
		points_data.push_back (std::pair <double, double> (x, y));
	}
	fclose (f);

	prepare (points_data, extra);
}

void Spline_1D_Smooth::prepare (std::vector<std::pair<double, double>> points_data, double * extra)
{
	int n = (int)(pow ((double)points_data.size (), 1.0 / 3.0));
	//n = (int)points_data.size ();
	//int n = 2;
	double c0 = 1e+10;
	double cN = -1e+10;
	
	for (int i = 0; i < (int)points_data.size (); i++)
	{
		if (points_data[i].first > cN)
			cN = points_data[i].first;
		if (points_data[i].first < c0)
			c0 = points_data[i].first;
	}
	//c0 -= 0.001 * (cN - c0);
	//cN += 0.001 * (cN - c0);
	prepare (c0, cN, n, extra[0], extra[1], points_data);
}

void Spline_1D_Smooth::prepare (char * file_name, double * extra, double n_points_sqrt)
{
	double x, y;
	std::vector <std::pair<double, double>> points_data;

	FILE * f = fopen (file_name, "r");
	while (!feof (f))
	{
		fscanf (f, "%lf %lf", &x, &y);
		points_data.push_back (std::pair <double, double> (x, y));
	}
	fclose (f);

	prepare (points_data, extra, n_points_sqrt);
}

void Spline_1D_Smooth::prepare (std::vector<std::pair<double, double>> points_data, double * extra, double n_points_sqrt)
{
	int n = (int)(pow ((double)points_data.size (), n_points_sqrt));
	//int n = 2;
	double c0 = 1e+10;
	double cN = -1e+10;

	for (int i = 0; i < (int)points_data.size (); i++)
	{
		if (points_data[i].first > cN)
			cN = points_data[i].first;
		if (points_data[i].first < c0)
			c0 = points_data[i].first;
	}
	//c0 -= 0.001 * (cN - c0);
	//cN += 0.001 * (cN - c0);
	prepare (c0, cN, n, extra[0], extra[1], points_data);
}

bool Spline_1D_Smooth::get_solution_in_point (double * c, double * value)
{
	return Task::get_solution_in_point (0, c, value);
}

void Spline_1D_Smooth::solve_task ()
{
	Task::solve_task (2, 0, 10);
}

void Spline_1D_Smooth::get_boundaries (double * c0, double * cN)
{
	mesh.get_0_boundaries (c0);
	mesh.get_N_boundaries (cN);
}

void Spline_1D_Smooth::prepare (double c0, double cN, int n_x, double Alpha, double Beta, std::vector<std::pair<double, double>> points_data)
{
	log_file = fopen ("Result Files Extra//log.txt", "w");

	// build a mesh
	mesh_pointer = &mesh;
	prepared = mesh.build_Mesh (c0, cN, n_x, &N_functions);
	
	solution_points = new MathVector ((int) points_data.size ());
	solution_values = new MathVector ((int) points_data.size ());

	alpha = Alpha;
	beta = Beta;

	// save data
	for (size_t i = 0, i_end = points_data.size (); i < i_end; i++)
	{
		solution_points->setElem ((int)i, points_data[i].first);
		solution_values->setElem ((int)i, points_data[i].second);
	}
	n_sol_points = solution_points->getSize ();

	// set amount of equations
	set_n_systems ();
	// set time layers
	set_time_layers ("");
	// prepares elements
	for (int i = 0, i_end = mesh_pointer->get_n_elements (); i < i_end; i++)
	{
		mesh_pointer->prepare_element (i);
	}

	use_time_mapping = true;

	// build matrices' portraits
	for (int k_system = 0; k_system < n_systems; k_system++)
		if (!build_portrait (k_system))
		{
			printf ("ERROR: portrait\n");
			prepared = false;
		}

	// if portrait building was successful, set starting conditions
	if (prepared)
		set_starting_conditions ();
}

void Spline_1D_Smooth::prepare (double c0, double cN, int n_x, double Alpha, double Beta, char * file_name)
{
	std::vector<std::pair<double, double>> points_data;
	FILE * file = fopen (file_name, "r");
	double x, f;
	while (!feof (file))
	{
		fscanf (file, "%lf %lf", &x, &f);
		points_data.push_back (std::pair <double, double> (x, f));
	}
	fclose (file);

	prepare (c0, cN, n_x, Alpha, Beta, points_data);
}
