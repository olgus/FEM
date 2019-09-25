#include "Spline_1D_Interpolative.h"

void Spline_1D_Interpolative::renumerate_functions ()
{
	int functions[4];
	int n_elements = (int)mesh.elements.size ();
	for (size_t i = 0; i < n_elements; i++)
	{
		functions[0] = (int)(n_elements + 1 + i);
		functions[1] = (int)i;
		functions[2] = (int)(n_elements + 1 + i + 1);
		functions[3] = (int)(i + 1);

		mesh.elements[i]->set_global_functions (functions);
	}

	N_functions = 2 * (n_elements + 1);
}

void Spline_1D_Interpolative::build_system (int k_system)
{
	// only for first N_functions / 2
	int der_functions = N_functions / 2;
	double a, h0, h1;

	for (int i = 1; i < der_functions - 1; i++)
	{
		// i - 1
		h0 = mesh.elements[i - 1]->get_geometrical_area ();
		a = 2.0 / h0;
		eq_matrixes[k_system].set_entry (i, i - 1, a);

		// i
		h1 = mesh.elements[i]->get_geometrical_area ();
		a = 4.0 * (1.0 / h0 + 1.0 / h1);
		eq_matrixes[k_system].set_entry (i, i, a);

		// i + 1
		a = 2.0 / h1;
		eq_matrixes[k_system].set_entry (i, i + 1, a);

		// f entry
		a = 6.0 * (-function_values[i - 1] / (h0 * h0) +
			function_values[i] * (1.0 / (h0 * h0) - 1.0 / (h1 * h1)) + function_values[i + 1] / (h1 * h1));
		
		eq_matrixes[k_system].set_f_entry (i, a);
	}
}

void Spline_1D_Interpolative::apply_boundary_conditions (int k_system)
{
	apply_first_boundary_conditions (k_system);
}

void Spline_1D_Interpolative::apply_first_boundary_conditions (int k_system)
{
	// every function after first N_functions / 2
	int der_functions = N_functions / 2;
	for (int i = der_functions; i < N_functions; i++)
	{
		eq_matrixes[k_system].set_entry (i, i, 1.0);
		eq_matrixes[k_system].set_f_entry (i, function_values[i - der_functions]);
	}
	// also add for functions 0 and N_functions / 2 - 1
	int i = 0;
	double h1 = mesh.elements[0]->get_geometrical_area ();
	double h2 = mesh.elements[1]->get_geometrical_area ();
	double a = (- function_values[0] * (2.0 * h1 + h2) / (h1 * (h1 + h2))
		+ function_values[1] * (h1 + h2) / (h1 * h2)
		- function_values[2] * (h1) / (h2 * (h1 + h2)));

	eq_matrixes[k_system].set_entry (i, i, 1.0);
	eq_matrixes[k_system].set_f_entry (i, a);

	i = der_functions - 1;
	int n_elements = (int)mesh.elements.size () - 1;

	h2 = mesh.elements[n_elements - 1]->get_geometrical_area ();
	h1 = mesh.elements[n_elements]->get_geometrical_area ();
	a = (function_values[i - 2] * ( h1 / (h2 * (h1 + h2)))
		- function_values[i - 1] * (h1 + h2) / (h1 * h2)
		+ function_values[i] * (2.0 * h1 + h2) / (h1 * (h1 + h2)));
	eq_matrixes[k_system].set_entry (i, i, 1.0);
	eq_matrixes[k_system].set_f_entry (i, a);

}

void Spline_1D_Interpolative::prepare (std::vector<std::pair<double, double>> points_data, double * extra)
{
	log_file = fopen ("Result Files Extra//log.txt", "w");
	// build a mesh
	mesh_pointer = &mesh;

	std::vector <double> coordinates;
	for (size_t i = 0, i_end = points_data.size (); i < i_end; i++)
	{
		coordinates.push_back (points_data[i].first);
		function_values.push_back (points_data[i].second);
	}
	prepared = mesh.build_Mesh (coordinates, &N_functions);

	// set amount of equations
	set_n_systems ();
	// set time layers
	set_time_layers ("");
	// prepares elements
	for (int i = 0, i_end = mesh_pointer->get_n_elements (); i < i_end; i++)
	{
		mesh_pointer->prepare_element (i);
	}

	// no conditions for spline

	use_time_mapping = true;

	// renumerate for better numbers
	renumerate_functions ();

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

void Spline_1D_Interpolative::prepare (char * file_name)
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

	double extra[1];
	prepare (points_data, extra);
}

bool Spline_1D_Interpolative::get_solution_in_point (double * c, double * value)
{
	return Task::get_solution_in_point (0, c, value);
}

void Spline_1D_Interpolative::solve_task ()
{
	Task::solve_task (2, 0, 1);
}

void Spline_1D_Interpolative::get_boundaries (double * c0, double * cN)
{
	mesh.get_0_boundaries (c0);
	mesh.get_N_boundaries (cN);
}

Spline_1D_Interpolative::Spline_1D_Interpolative ()
{
}

Spline_1D_Interpolative::Spline_1D_Interpolative (const Spline_1D_Interpolative & spline)
{
}

Spline_1D_Interpolative::~Spline_1D_Interpolative ()
{
}
