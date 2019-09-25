#include "Task_3D_point_source.h"

double Task_3D_point_source::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	return 0.0;
}

void Task_3D_point_source::get_local_B (int k_system, int k_element, MathVector * B)
{
	B->Zero ();
}

void Task_3D_point_source::apply_boundary_conditions (int k_system)
{
	apply_first_boundary_conditions (k_system);

	// add sources
	double coordinates[3];
	MathVector * B;
	int n_functions;
	int iF;
	bool counted;

	for (int i = 0; i < n_sources; i++)
	{
		coordinates[0] = sources[i].x;
		coordinates[1] = sources[i].y;
		coordinates[2] = sources[i].z;

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

void Task_3D_point_source::prepare (char * mesh_file_name, char * bound_file_name, bool save)
{
	Task::prepare (mesh_file_name, bound_file_name, save);
	// get sources
	FILE * file = fopen ("Source Files\\point_sources.txt", "r");
	fscanf (file, "%i", &n_sources);
	if (n_sources > 0)
	{
		sources = new point_source[n_sources];
		for (int i = 0; i < n_sources; i++)
		{
			fscanf (file, "%lf %lf %lf %lf", &sources[i].x, &sources[i].y, &sources[i].z, &sources[i].power);
		}
	}
	fclose (file);
}

void Task_3D_point_source::prepare (char * file_name_nodes, char * file_name_elements, char * bound_file_name)
{
	Task::prepare (file_name_nodes, file_name_elements, bound_file_name);
	// get sources
	FILE * file = fopen ("Source Files\\point_sources.txt", "r");
	fscanf (file, "%i", &n_sources);
	if (n_sources > 0)
	{
		sources = new point_source[n_sources];
		for (int i = 0; i < n_sources; i++)
		{
			fscanf (file, "%lf %lf %lf %lf", &sources[i].x, &sources[i].y, &sources[i].z, &sources[i].power);
		}
	}
	fclose (file);
}

double Task_3D_point_source::Gamma (int k_system, int k_function, int area)
{
	return 0.0;
}

Task_3D_point_source::Task_3D_point_source ()
{
	n_sources = 0;
	sources = NULL;
}

Task_3D_point_source::~Task_3D_point_source ()
{
	if (sources != NULL)
		delete[] sources;
}
