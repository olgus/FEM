#include "Spline.h"

void Spline::get_local_A_matrix (int k_element, Matrix * A_matrix)
{
	// set A_matrix to zero
	A_matrix->Zero ();
	// allocate some memory
	int matr_dim = mesh. elements[k_element]->get_amount_non_zero_functions ();
	Matrix * A = new Matrix (matr_dim, matr_dim);

	int dim = mesh.get_dimentionality ();
	double * coordinates = new double[dim];
	double val;

	int m[] = { 1 };
	// go through solution points
	for (int k = 0; k < n_sol_points; k++)
	{
		for (int j = 0; j < dim; j++)
		{
			coordinates[j] = solution_points->Elem (k, j);
		}

		// check if it's in the current element 
		if (mesh. elements[k_element]->point_inside (mesh, coordinates))
		{
			// if it is, get matrix A[ij] = Psi[i] (solution_points[k]) * Psi[j] (solution_points[k])
			for (int i = 0; i < matr_dim; i++)
			{
				for (int j = 0; j < matr_dim; j++)
				{
					val = mesh. elements[k_element]->get_function_value (i, j, m, coordinates);
					A->setElem (i, j, val);
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

void Spline::get_local_B (int k_element, MathVector * B_vector)
{
	// set B_vector to zero
	B_vector->Zero ();
	// allocate some memory
	int matr_dim = mesh. elements[k_element]->get_amount_non_zero_functions ();
	MathVector * B = new MathVector (matr_dim);

	int dim = mesh.get_dimentionality ();
	double * coordinates = new double[dim];
	double val;
	// go through solution points
	for (int k = 0; k < n_sol_points; k++)
	{
		for (int j = 0; j < dim; j++)
		{
			coordinates[j] = solution_points->Elem (k, j);
		}

		// check if it's in the current element 
		if (mesh. elements[k_element]->point_inside (mesh, coordinates))
		{
			// if it is, get Psi[i] (solution_points[k]) 
			for (int i = 0; i < matr_dim; i++)
			{
				val = mesh. elements[k_element]->get_basis_function_value (i, coordinates);
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

double Spline::function_f (int k_system, double * coordinates, int area)
{
	return 0.0;
}

void Spline::build_system (int k_system)
{
	Matrix * A;
	Matrix * G;
	//Matrix * D;
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
		//D = new Matrix (n_entries, n_entries);
		B = new MathVector (n_entries);

		get_local_A_matrix (k_element, A);
		mesh.elements[k_element]->get_G_local_matrix (G);
		//mesh.elements[k_element]->get_D_local_matrix (D);
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
				//eq_matrixes[k_system].add_to_entry (iB, jB, beta * D->Elem (i, j));
			}
			eq_matrixes[k_system].add_to_f_entry (iB, B->getElem (i));
		}
		delete[] def_nodes;
		delete A;
		delete G;
		//delete D;
		delete B;
	}


	fclose (fileG);

}

void Spline::apply_boundary_conditions (int k_system)
{
	// no conditions
}

void Spline::print_solutions ()
{
	FILE * file_1 = fopen ("Result//spline_F.txt", "w");
	FILE * file_2 = fopen ("Result//spline_dF_dx.txt", "w");
	FILE * file_3 = fopen ("Result//spline_dF_dy.txt", "w");
	FILE * file_4 = fopen ("Result//spline_dF2_dxdy.txt", "w");
	FILE * file_5 = fopen ("Result//spline_all.txt", "w");

	// 4 files
	// f
	// df/dx
	// df/dy
	// df^2/(dx*dy)
	MathVector * v = new MathVector (eq_matrixes[0].Size());
	eq_matrixes[0].get_solution (v);

	for (int i = 0, i_end = mesh_pointer->get_n_nodes (); i < i_end; i++)
	{
		fprintf (file_1, "%.16lf\n", v->getElem (i * 4));
		fprintf (file_2, "%.16lf\n", v->getElem (i * 4 + 1));
		fprintf (file_3, "%.16lf\n", v->getElem (i * 4 + 2));
		fprintf (file_4, "%.16lf\n", v->getElem (i * 4 + 3));
	}

	for (int i = 0, i_end = eq_matrixes[0].Size (); i < i_end; i++)
	{
		fprintf (file_5, "%.16lf\n", v->getElem (i));
	}
	fclose (file_1);
	fclose (file_2);
	fclose (file_3);
	fclose (file_4);
	fclose (file_5);
	delete v;
}

void Spline::get_conditions (int k_system, char * file_name)
{
	// no conditions
}
 
bool Spline::get_solution_in_point (int k_system, double * coordinates, double * value)
{
	bool found = false;
	int p_element;
	// find element that contains the point 
	for (int i = 0, i_end = mesh_pointer->get_n_elements (); i < i_end; i++)
	{
		if (mesh_pointer->point_inside (i, coordinates))
		{
			p_element = i;
			found = true;
		}
	}
	if (!found)
	{
		*value = 0.0;
		return false;
	}

	// get local functions' values of the element
	int node_amount = mesh_pointer->get_amount_of_def_nodes (p_element);
	int func_amount = mesh_pointer->get_amount_non_zero_functions (p_element);
	MathVector * local_func = new MathVector (func_amount);
	mesh_pointer->get_local_function_values (p_element, coordinates, local_func);

	// make vector of solutions in those nodes
	MathVector * solution = new MathVector (func_amount);
	int * def_nodes = new int[node_amount];
	mesh_pointer->get_def_nodes (p_element, def_nodes);
	for (int i = 0; i < node_amount; i++)
	{
		for (int k = 0; k < 4; k++)
		{
			solution->setElem (i * 4 + k, previous_time_layers_solutions[k_system][0].getElem (def_nodes[i] * 4 + k));
		}
	}

	// multiply them by the solution's values in those nodes
	*value = local_func->Scalar_Product (*solution);

	delete local_func;
	delete solution;
	delete[] def_nodes;

	return true;
}

bool Spline::get_first_derivative (int k_system, int k_var, double * coordinates, double * value)
{
	bool found = false;
	int p_element;
	// find element that contains the point 
	for (int i = 0, i_end = mesh_pointer->get_n_elements (); i < i_end; i++)
	{
		if (mesh_pointer->point_inside (i, coordinates))
		{
			p_element = i;
			found = true;
		}
	}
	if (!found)
	{
		*value = 0.0;
		return false;
	}

	// get local functions' values of the element
	int node_amount = mesh_pointer->get_amount_of_def_nodes (p_element);
	int func_amount = mesh_pointer->get_amount_non_zero_functions (p_element);
	MathVector * local_func = new MathVector (func_amount);
	// get values of k_var derivative 

	for (int i = 0; i < func_amount; i++)
	{
		local_func->setElem (i, mesh.elements[p_element]->get_basis_function_derivative (i, k_var, coordinates));
	}

	// make vector of solutions in those nodes
	MathVector * solution = new MathVector (func_amount);
	int * def_nodes = new int[node_amount];
	mesh_pointer->get_def_nodes (p_element, def_nodes);
	for (int i = 0; i < node_amount; i++)
	{
		for (int k = 0; k < 4; k++)
		{
			solution->setElem (i * 4 + k, previous_time_layers_solutions[k_system][0].getElem (def_nodes[i] * 4 + k));
		}
	}
	// multiply them by the solution's values in those nodes
	*value = local_func->Scalar_Product (*solution);
	delete local_func;
	delete solution;
	delete[] def_nodes;

	return true;
}

void Spline::prepare (char * mesh_file_name, char * bound_file_name, bool save)
{
	Task::prepare (mesh_file_name, bound_file_name, save);

	// get points to smooth
	FILE * file = fopen ("Result Files Extra//points_amount.txt", "r");
	fscanf (file, "%i", &n_sol_points);
	printf ("Spline for %i points\n", n_sol_points);
	fclose (file);

	// get alpha and beta
	file = fopen ("Source Files//spline.txt", "r");
	fscanf (file, "%lf %lf", &alpha, &beta);
	fclose (file);

	// allocate memory
	int dim = mesh_pointer->get_dimentionality ();
	solution_values = new MathVector (n_sol_points);
	solution_points = new Matrix (n_sol_points, dim);
	weigths = new MathVector (n_sol_points);

	// get spline points from file
	double val;
	file = fopen ("Source Files//spline_points.txt", "r");
	for (int i = 0; i < n_sol_points; i++)
	{
		// solution point's coordinates
		for (int j = 0; j < dim; j++)
		{
			fscanf (file, "%lf", &val);
			solution_points->setElem (i, j, val);
		}
		// solution's value
		fscanf (file, "%lf", &val);
		solution_values->setElem (i, val);
	}
	fclose (file);

	double * w = new double[n_sol_points];
	// get weigths
	//file = fopen ("Spline Points//weights.txt", "r");
	//if (!feof (file)) // if they are in the file, read them
	//{
	//	for (int i = 0; i < n_sol_points; i++)
	//	{
	//		fscanf (file, "%lf", &val);
	//		w[i] = val;
	//	}
	//}
	//else // or set them as 1.0 (default value)
	{
		for (int i = 0; i < n_sol_points; i++)
		{
			w[i] = 1.0;
		}
	}
	weigths->Initialize (w);
	//fclose (file);
	delete[] w;
}

void Spline::prepare (Mesh_Prototype * parent_mesh, char * spline_for_mesh_file)
{
	printf ("\tpreparing time layers and elements\n");

	log_file = fopen ("Result Files Extra//log.txt", "w");
	// get c0 and cN from parent_mesh
	double c0[2];
	double cN[2];
	parent_mesh->get_0_boundaries (c0);
	parent_mesh->get_N_boundaries (cN);

	// TEST
	alpha = beta = 1e-3;

	// get points
	FILE * file = fopen (spline_for_mesh_file, "r");
	fscanf (file, "%i", &n_sol_points);
	double min = 1e+15;
	double max = -1e+15;
	double value;
	solution_values = new MathVector (n_sol_points);
	solution_points = new Matrix (n_sol_points, 2);
	weigths = new MathVector (n_sol_points);

	for (int i = 0; i < n_sol_points; i++)
	{
		for (int k = 0; k < 2; k++)
		{
			fscanf (file, "%lf", &value);
			solution_points->setElem (i, k, value);
		}

		// solution's value
		fscanf (file, "%lf", &value);
		solution_values->setElem (i, value);
		if (value < min)
			min = value;
		if (value > max)
			max = value;
	}

	// amount of nodes on sides = sqrt (n_points)
	int N_axis[2];
	N_axis[0] = N_axis[1] = (int)(sqrt (n_sol_points)) + 1;

	// set weights
	double * w = new double[n_sol_points];
	// 1.0 for most nodes, 2.0 for minimal values
	double div_value = (max - min) * 0.25 + min;
	for (int i = 0; i < n_sol_points; i++)
	{
		if (solution_values->getElem (i) < div_value + 1e-10)
			w[i] = 2.0;
		else
			w[i] = 1.0;

	}
	weigths->Initialize (w);

	// build mesh
	mesh_pointer = &mesh;
	prepared = mesh.build_Mesh (c0, cN, N_axis, &N_functions);
	// prepares elements
	for (int i = 0, i_end = mesh_pointer->get_n_elements (); i < i_end; i++)
	{
		mesh_pointer->prepare_element (i);
	}	
	// set amount of equations
	set_n_systems ();
	// set time layers
	set_time_layers ("");
	delete[] w;
}

void Spline::prepare (double * c0, double * cN, double lvl, char * spline_for_mesh_file)
{
	//int n_points_file, n_gen_points;
	printf ("\tpreparing time layers and elements\n");

	log_file = fopen ("Result Files Extra//log.txt", "w");

	// TEST
	alpha = 1e-2;
	beta = 1e-7;

	// get points
	FILE * file = fopen (spline_for_mesh_file, "r");
	fscanf (file, "%i", &n_sol_points);
	//n_gen_points = n_points_file * n_points_file;
	//n_sol_points = n_points_file + n_gen_points;

	double min = 1e+15;
	double max = -1e+15;
	double value;
	solution_values = new MathVector (n_sol_points);
	solution_points = new Matrix (n_sol_points, 2);
	weigths = new MathVector (n_sol_points);

	for (int i = 0; i < n_sol_points; i++)
	{
		for (int k = 0; k < 2; k++)
		{
			fscanf (file, "%lf", &value);
			solution_points->setElem (i, k, value);
		}

		// solution's value
		fscanf (file, "%lf", &value);
		value /=  pow (2.0, lvl);
		solution_values->setElem (i, value);
		if (value < min)
			min = value;
		if (value > max)
			max = value;
	}

	//// generate points
	//{
	//	double hx = (cN[0] - c0[0]) / (double) n_points_file;
	//	double hy = (cN[1] - c0[1]) / (double) n_points_file;
	//	double x, y;
	//	int counter = 0;
	//	for (int i = 0; i < n_points_file; i++)
	//	{
	//		x = c0[0] * i + hx / 2.0;
	//		for (int j = 0; j < n_points_file; j++)
	//		{
	//			y = c0[1] * j + hy / 2.0;
	//			solution_points->setElem (counter, 0, x);
	//			solution_points->setElem (counter, 1, y);
	//			solution_values->setElem (counter, def_value);
	//			counter++;
	//		}
	//	}
	//}

	// amount of nodes on sides = sqrt (n_points)
	int N_axis[2];
	N_axis[0] = N_axis[1] = (int)(/*sqrt*/ (n_sol_points));

	// set weights
	double * w = new double[n_sol_points];
	// 1.0 for most nodes, 2.0 for minimal values
	double div_value = (max - min) * 0.25 + min;
	for (int i = 0; i < n_sol_points; i++)
	{
		if (solution_values->getElem (i) < div_value + 1e-10)
			w[i] = 10.0;
		else
			w[i] = 1.0;

	}
	//for (int i = n_points_file; i < n_sol_points; i++)
	//	w[i] = 1.0;

	weigths->Initialize (w);

	// build mesh
	mesh_pointer = &mesh;
	prepared = mesh.build_Mesh (c0, cN, N_axis, &N_functions);
	// prepares elements
	for (int i = 0, i_end = mesh_pointer->get_n_elements (); i < i_end; i++)
	{
		mesh_pointer->prepare_element (i);
	}
	// set amount of equations
	set_n_systems ();
	// set time layers
	set_time_layers ("");
	
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
	delete[] w;
}

void Spline::save_slice_first_derivative (int k_system, char * file_sour, char * file_dest, bool keep_points)
{
	int dim_node = mesh_pointer->get_dimentionality ();
	double * coordinates = new double[dim_node];
	double * coordinatesN = new double[dim_node];
	double value; // value of the solution
	bool exists; // flag of existence of solution's value 
	int n_amount = 0; // amount of points written in file 
	int k_var; // variable to iterate through
	double h; // iteration step
	int iter_amount; // amount of points to iterate through
	double vstart, vend; // starting and ending values of variable 

	FILE * points = fopen (file_sour, "r");
	fscanf (points, "%i %lf", &k_var, &h); // get number of the variable to iterate through and step for it
	if (k_var >= 0 && k_var < dim_node) // if the variable number is valide
	{
		// get other variables' values (for slice)
		for (int i = 0; i < dim_node; i++)
		{
			if (i != k_var)
			{
				fscanf (points, "%lf", &coordinates[i]);
			}
		}
		// get starting value
		mesh_pointer->get_0_boundaries (coordinatesN);
		vstart = coordinatesN[k_var];
		// get ending value
		mesh_pointer->get_N_boundaries (coordinatesN);
		vend = coordinatesN[k_var];
		// get amount of points to iterate through
		iter_amount = (int)((vend - vstart) / h);

		// open file for result input
		FILE * result_points = fopen (file_dest, "w");
		for (int i = 0; i < iter_amount + 1; i++)
		{
			coordinates[k_var] = vstart + i * h;
			//calculate the point
			value = 0.0;
			// get value in the point
			exists = get_first_derivative (k_system, k_var, coordinates, &value);
			// if value has to be saved with the point
			if (keep_points)
			{
				// then print in out only of it exists
				if (exists)
				{
					for (int i = 0; i < dim_node; i++)
					{
						fprintf (result_points, "%.16lf ", coordinates[i]);
					}
					fprintf (result_points, "%.16lf\n", value);
					n_amount++;
				}
			}
			else // otherwise print it anyway, it'll be zero
			{
				fprintf (result_points, "%.16lf\n", value);
				n_amount++;
			}
		}
		fclose (result_points);
	}

	fclose (points);

	delete[] coordinates;
	delete[] coordinatesN;

	printf ("\trequested solution is saved\n");
}

void Spline::save_slice_mixed_derivative (int k_system, char * file_sour, char * file_dest, bool keep_points)
{
	int dim_node = mesh_pointer->get_dimentionality ();
	double * coordinates = new double[dim_node];
	double * coordinatesN = new double[dim_node];

	double value; // value of the solution
	bool exists; // flag of existence of solution's value 
	int n_amount = 0; // amount of points written in file 
	int k_var; // variable to iterate through
	double h; // iteration step
	int iter_amount; // amount of points to iterate through
	double vstart, vend; // starting and ending values of variable 

	FILE * points = fopen (file_sour, "r");
	fscanf (points, "%i %lf", &k_var, &h); // get number of the variable to iterate through and step for it
	if (k_var >= 0 && k_var < dim_node) // if the variable number is valide
	{
		// get other variables' values (for slice)
		for (int i = 0; i < dim_node; i++)
		{
			if (i != k_var)
			{
				fscanf (points, "%lf", &coordinates[i]);
			}
		}
		// get starting value
		mesh_pointer->get_0_boundaries (coordinatesN);
		vstart = coordinatesN[k_var];
		// get ending value
		mesh_pointer->get_N_boundaries (coordinatesN);
		vend = coordinatesN[k_var];
		// get amount of points to iterate through
		iter_amount = (int)((vend - vstart) / h);

		// open file for result input
		FILE * result_points = fopen (file_dest, "w");
		for (int i = 0; i < iter_amount + 1; i++)
		{
			coordinates[k_var] = vstart + i * h;
			//calculate the point
			value = 0.0;
			// get value in the point
			exists = get_first_derivative (k_system, 10, coordinates, &value);
			// if value has to be saved with the point
			if (keep_points)
			{
				// then print in out only of it exists
				if (exists)
				{
					for (int i = 0; i < dim_node; i++)
					{
						fprintf (result_points, "%.16lf ", coordinates[i]);
					}
					fprintf (result_points, "%.16lf\n", value);
					n_amount++;
				}
			}
			else // otherwise print it anyway, it'll be zero
			{
				fprintf (result_points, "%.16lf\n", value);
				n_amount++;
			}
		}
		fclose (result_points);
	}

	fclose (points);

	delete[] coordinates;
	delete[] coordinatesN;

	printf ("\trequested solution is saved\n");
}

void Spline::default_value (double value)
{
	def_value = value;
}

void Spline::set_alpha_beta (double Alpha, double Beta)
{
	alpha = Alpha;
	beta = Beta;
}

Spline::Spline ()
{
	def_value = 1.0;
	n_sol_points = 0;
	alpha = beta = 1e-5;
	solution_values = NULL;
	solution_points = NULL;
	weigths = NULL;
}

Spline::Spline (const Spline & spline)
{
	n_sol_points = spline.n_sol_points;
	alpha = spline.alpha;
	beta = spline.beta;
			
	if (solution_values != NULL)
		delete solution_values;
	if (solution_points != NULL)
		delete solution_points;
	if (weigths != NULL)
		delete weigths;

	if (spline.solution_values != NULL)
	{
		solution_values = new MathVector (spline.solution_values->getSize ());
		solution_values->Copy (*(spline.solution_values));
	}
	if (spline.solution_points != NULL)
	{
		solution_points = new Matrix (spline.solution_points->Size0 (), spline.solution_points->Size1 ());
		solution_points->Copy (*(spline.solution_points));
	}
	if (spline.weigths != NULL)
	{
		weigths = new MathVector (spline.weigths->getSize ());
		weigths->Copy (*(spline.weigths));
	}
}

Spline::~Spline ()
{
	if (solution_values != NULL)
		delete solution_values;
	if (solution_points != NULL)
		delete solution_points;
	if (weigths != NULL)
		delete weigths;
}
