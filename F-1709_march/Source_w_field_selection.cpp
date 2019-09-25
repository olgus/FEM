#include "Source_w_field_selection.h"

double Task_2D_source_w_field_selection::function_f (int k_system, double * coordinates, int area)
{
	return 0.0;
}

void Task_2D_source_w_field_selection::apply_boundary_conditions (int k_system)
{
}

double Task_2D_source_w_field_selection::get_main_field (double * coordinates)
{
	double v = 0.0;
	double r;
	for (int i = 0; i < n_sources; i++)
	{
		r = sqrt (pow (coordinates[0] - sources[i].x, 2.0) + pow (coordinates[1] - sources[i].y, 2.0));
		if (fabs (r) < 1e-10)
			v += 1e+10;
		else
		{
			v += 1.0 - (sources[i].power / (2.0 * PI * Lambda (0, 0, lambda_val[1]))) * log (r / 1e-5);
		}
	}
	//printf ("%lf\t", lambda_val[1]);
	return v;
}

Task_2D_source_w_field_selection::Task_2D_source_w_field_selection ()
{
	n_sources = 0;
	sources = NULL;
	mesh_pointer = NULL;
	tasks_additional_field = NULL;
}

Task_2D_source_w_field_selection::~Task_2D_source_w_field_selection ()
{
	if (sources != NULL)
		delete[] sources;
	if (tasks_additional_field != NULL)
		delete[] tasks_additional_field;
}

void Task_2D_source_w_field_selection::prepare (int N_insertions, char * file_name_sources, char * file_name_lambdas, char * file_name_insertion_nodes[], char * file_name_insertion_triangles[])
{
	n_insertions = N_insertions;

	// get sources
	FILE * file = fopen (file_name_sources, "r");
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
	
	// read lambdas
	lambda_val.clear ();
	file = fopen (file_name_lambdas, "r");

	int area;
	double value;
	while (!feof (file))
	{
		fscanf (file, "%i %lf", &area, &value);
		lambda_val.insert (std::pair<int, double> (area, value));
	}

	fclose (file);

	tasks_additional_field = new Task_2D_Additional_field[n_insertions];
	Task_pointer * cur = this;
	for (int i = 0; i < n_insertions; i++)
	{
		tasks_additional_field[i].prepare (i + 2, cur, file_name_insertion_nodes[i], file_name_insertion_triangles[i], lambda_val);
	}
}

void Task_2D_source_w_field_selection::solve_task ()
{
	Task_pointer * cur = this;
	// for each insertion
	for (int i = 0; i < n_insertions; i++)
	{
		tasks_additional_field[i].solve_task (2, 0, 0);
	}
}

bool Task_2D_source_w_field_selection::get_solution_in_point (int k_system, double * coordinates, double * value)
{
	*value = get_main_field (coordinates);
	if (k_system == 0)
	{
		//*value = 0.0;
		for (int i = 0; i < n_insertions; i++)
		{
			double ins_value = 0.0;
			tasks_additional_field[i].get_solution_in_point (0, coordinates, &ins_value);
			*value += ins_value;
		}
	}
	else
	{
		for (int i = 0; i < k_system - 2; i++)
		{
			double ins_value = 0.0;
			tasks_additional_field[i].get_solution_in_point (0, coordinates, &ins_value);
			*value += ins_value;
		}
	}

	return true;
}

int Task_2D_source_w_field_selection::get_mesh_dim ()
{
	return 2;
}
