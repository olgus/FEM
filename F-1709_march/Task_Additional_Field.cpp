#include "Task_Additional_Field.h"

void Task_2D_Additional_field::get_local_B (int k_system, int k_element, MathVector * B)
{
	int n_funcs = mesh_pointer->get_amount_non_zero_functions (k_element);

	double jac = 0.0; // jacobian
	int dim = mesh_pointer->get_dimentionality ();
	int n_integr_points = mesh.elements[k_element]->edge_amount_of_integration_points ();
	double * weigths = new double[n_integr_points];
	double ** points = new double *[n_integr_points];
	for (int i = 0; i < n_integr_points; i++)
	{
		points[i] = new double[dim];
	}

	// get integration points from the element
	mesh.elements[k_element]->integration_points (points, weigths, &jac);

	// make matrix 
	Matrix * M = new Matrix (n_funcs, n_funcs);
	MathVector * Q = new MathVector (n_funcs);
	mesh.get_M_local_matrix (k_element, M);

	double val;
	double main_value;
	for (int i = 0; i < n_funcs; i++)
	{
		// get vector of integrals f * u0
		val = 0.0;
		// get theta values * basis_functions in those points
		for (int k = 0; k < n_integr_points; k++)
		{
			task_main_field->get_solution_in_point (k_insertion, points[k], &main_value);
			// sum them multiplying by weigths
			val += mesh.elements[k_element]->get_basis_function_value (i, points[k]) * main_value * weigths[k];
		}
		// multiply by jac
		val *= jac;
		B->setElem (i, val);
	}
	// solve it
	M->solve (Q, *B);

	// multiply by G
	Matrix * G = new Matrix (n_funcs, n_funcs);
	mesh.get_G_local_matrix (k_element, G);
	G->MultiplyByNumber (lambda_val[1] - lambda_val[mesh_pointer->get_area (k_element)]);
	//printf ("%lf %i %lf\n", lambda_val[1], mesh_pointer->get_area (k_element), lambda_val[mesh_pointer->get_area (k_element)]);
	G->MultiplyMatrixByVector (*Q, B);

	delete Q;
	delete G;
	delete M;
	delete[] weigths;
	for (int i = 0; i < n_integr_points; i++)
	{
		delete[] points[i];
	}
	delete[] points;

}

double Task_2D_Additional_field::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	// to set starting vector to 0
	return 0.0;
}

Task_2D_Additional_field::Task_2D_Additional_field ()
{
}

Task_2D_Additional_field::~Task_2D_Additional_field ()
{
}

void Task_2D_Additional_field::prepare (int K_insertion, Task_pointer * Task_main_field, char * file_name_nodes, char * file_name_triangles, std::map <int, double> Lambda_val)
{
	k_insertion = K_insertion;
	task_main_field = Task_main_field;
	log_file = fopen ("Result Files Extra//log.txt", "w");


	int N_axis[] = { 20, 20 };
	// make uniform mesh, all areas = 1
	mesh.build_Mesh (file_name_nodes, file_name_triangles, &N_functions);
	mesh_pointer = &mesh;
	// set amount of equations
	set_n_systems ();
	// set time layers
	set_time_layers ("");
	// prepares elements
	for (int i = 0, i_end = mesh_pointer->get_n_elements (); i < i_end; i++)
	{
		mesh_pointer->prepare_element (i);
	}
	// build matrices' portraits
	for (int k_system = 0; k_system < n_systems; k_system++)
		if (!build_portrait (k_system))
		{
			printf ("ERROR: portrait\n");
			prepared = false;
		}
	set_starting_conditions ();

	// gamma 0
	gamma_val.clear ();
	for (int i = 0; i < k_insertion + 1; i++)
		gamma_val.insert (std::pair<int, double> (i + 1, 0));
	// copy lambda
	lambda_val.clear ();
	for (int i = 0; i < k_insertion + 1; i++)
		lambda_val.insert (std::pair<int, double> (i + 1, Lambda_val[i + 1]));
	// first boundary condition
	// set conditions
	conditions = new Condition *[n_systems];
	for (int i = 0; i < n_systems; i++)
	{
		conditions[i] = new Condition[mesh_pointer->get_dimentionality () * 2];
		for (int k = 0; k < 4; k++)
		{
			conditions[i][k].set_Type (1);
			conditions[i][k].set_function (1);
		}
	}

	prepared = true;
}