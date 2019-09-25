#include "Task_Mesh_Interpolation.h"

void Task_2D_Mesh_Interpolation::apply_boundary_conditions (int k_system)
{
	// no boundary conditions
}

double Task_2D_Mesh_Interpolation::function_f (int k_system, double * coordinates, int area)
{
	double r = 0.0;
	task->get_solution_in_point (task_k_system, coordinates, &r);
	return r;
}

Task_2D_Mesh_Interpolation::Task_2D_Mesh_Interpolation ()
{
	task = NULL;

	// lambda zero
	lambda_val.clear ();
	lambda_val.insert (std::pair<int, double> (1, 0));
}

Task_2D_Mesh_Interpolation::~Task_2D_Mesh_Interpolation ()
{
}

void Task_2D_Mesh_Interpolation::prepare (Task_pointer * Task, int K_system, double * c0, double * cN, int * N_axis)
{
	log_file = fopen ("Result Files Extra//log.txt", "w");

	task = Task;
	task_k_system = K_system;

	// make uniform mesh, all areas = 1
	N_functions = mesh.build_mesh (c0, cN, N_axis);
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
	prepared = true;
}
