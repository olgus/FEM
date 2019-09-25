#include "Task_2D_point_source.h"

double Task_2D_point_source::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	//return coordinates[0] + coordinates[1];
	return 0.0;
}

double Task_2D_point_source::function_f (int k_system, double * coordinates, int area)
{
	//return coordinates[0] + coordinates[1];
	return 0.0;
}

double Task_2D_point_source::function_Theta (int k_system, double * coordinates, int area, int boundary)
{
	double x = coordinates[0];
	double y = coordinates[1];
	double res = 0.0;
	return res;
}

void Task_2D_point_source::apply_boundary_conditions (int k_system)
{
	apply_first_boundary_conditions (k_system);
	double coordinates[2];
	int iF;

	// if r if bigger than 1000, zero it
	double r;
	int counter;
	int nodes[3];
	for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
	{
		counter = 0;
		// get def_nodes
		mesh_pointer->get_def_nodes (k_element, nodes);
		for (int i = 0; i < mesh_pointer->get_amount_of_def_nodes (0); i++)
		{
			// get coordinates
			mesh_pointer->get_node_coordinates (nodes[i], coordinates);
			// get r there
			r = sqrt (coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1]);
			// if smaller, increase counter
			if (r > 1000.0)
				counter++;
		}
		if (counter == 3)
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
	// add sources
	MathVector * B;
	int n_functions;
	bool counted;

	for (int i = 0; i < n_sources; i++)
	{
		coordinates[0] = sources[i].x;
		coordinates[1] = sources[i].y;

		counted = false;
		for (int k_element = 0, k_end = mesh.get_n_elements (); k_element < k_end && !counted; k_element++)
		{
			// find element that has it
			if (mesh.point_inside (k_element, coordinates))
			{
				n_functions = mesh.get_amount_non_zero_functions (k_element);
				B = new MathVector (n_functions);
				// get basis functions values in point of source
				for (int j = 0; j < n_functions; j++)
				{
					B->setElem (j, mesh. elements[k_element]->get_basis_function_value (j, coordinates));
				}
				// make vector of basis_functions_values * power
				B->MultiplyByNumber (sources[i].power);
				// add it respectively into F
				for (int j = 0; j < n_functions; j++)
				{
					iF = get_function_global_number (k_element, j);
					if (iF != -1)
					{
						eq_matrixes[k_system].add_to_f_entry (iF, B->getElem (j));
					}
				}
				delete B;
				counted = true;
			}
		}
	}
}

void Task_2D_point_source::prepare (char * mesh_file_name, char * bound_file_name, bool save)
{
	Task::prepare (mesh_file_name, bound_file_name, save);
	// get sources
	FILE * file = fopen ("Source Files\\point_sources.txt", "r");
	fscanf (file, "%i", &n_sources);
	if (n_sources > 0)
	{
		sources = new point_source_2D[n_sources];
		for (int i = 0; i < n_sources; i++)
		{
			fscanf (file, "%lf %lf %lf", &sources[i].x, &sources[i].y, &sources[i].power);
		}
	}
	fclose (file);
}

void Task_2D_point_source::prepare (char * file_name_nodes, char * file_name_elements, char * bound_file_name)
{
	Task::prepare (file_name_nodes, file_name_elements, bound_file_name);
	// get sources
	FILE * file = fopen ("Source Files\\point_sources.txt", "r");
	fscanf (file, "%i", &n_sources);
	if (n_sources > 0)
	{
		sources = new point_source_2D[n_sources];
		for (int i = 0; i < n_sources; i++)
		{
			fscanf (file, "%lf %lf %lf", &sources[i].x, &sources[i].y, &sources[i].power);
		}
	}
	fclose (file);
}

void Task_2D_point_source::prepare (char * file_name_nodes, char * file_name_elements, char * bound_file_name, char * file_sources)
{
	Task::prepare (file_name_nodes, file_name_elements, bound_file_name);
	// get sources
	FILE * file = fopen (file_sources, "r");
	fscanf (file, "%i", &n_sources);
	if (n_sources > 0)
	{
		sources = new point_source_2D[n_sources];
		for (int i = 0; i < n_sources; i++)
		{
			fscanf (file, "%lf %lf %lf", &sources[i].x, &sources[i].y, &sources[i].power);
		}
	}
	fclose (file);
}

Task_2D_point_source::Task_2D_point_source ()
{
	n_sources = 0;
	sources = NULL;

	gamma_val.clear ();
	lambda_val.clear ();
	for (int i = 0; i < 10; i++)
	{
		gamma_val.insert (std::pair<int, double> (i, 0.0));
		lambda_val.insert (std::pair<int, double> (i, 1.0));
	}
}

Task_2D_point_source::~Task_2D_point_source ()
{
	if (sources != NULL)
		delete[] sources;
}
