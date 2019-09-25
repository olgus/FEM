#include "Test_1D_Task.h"

template class Task_Test_1D <Mesh<Node<Point_Prototype>, Element>>;
template class Task_Test_1D <Mesh<Node<Point_2D>, Element>>;
template class Task_Test_1D <Mesh<Node<Point_3D>, Element>>;
template class Task_Test_1D <Mesh<Node<Point_2D_Polar>, Element>>;
template class Task_Test_1D <Mesh<Node_2D, Element>>;
template class Task_Test_1D <Mesh<Node_2D_Polar, Element>>;
template class Task_Test_1D <Mesh<Node_1D, Element>>;
template class Task_Test_1D <Mesh<Node_3D, Element>>;

template class Task_Test_1D <Mesh_1D_L1>;
template class Task_Test_1D <Mesh_1D_Hier>;
template class Task_Test_1D <Mesh_1D_Hermitian>;
template class Task_Test_1D <Triangular_Mesh>;
template class Task_Test_1D <Mixed_Triangular_Mesh>;
template class Task_Test_1D <Triangular_Polar_Mesh>;
template class Task_Test_1D <Triangular_Mesh_Hier>;
template class Task_Test_1D <Rectangular_Mesh>;
template class Task_Test_1D <Rectangular_Mesh_3>;
template class Task_Test_1D <Rectangular_Mesh_Spline>;
template class Task_Test_1D <Prismatic_Mesh>;
template class Task_Test_1D <Cubic_Mesh>;

template<class Type_Mesh>
void Task_Test_1D<Type_Mesh>::apply_boundary_conditions (int k_system)
{
	// set 0 function
	int iF = 0;
	int area_element = 0;
	int area = mesh_pointer->get_area (area_element);
	double coordinates[1];
	mesh_pointer->get_0_boundaries (coordinates);

	eq_matrixes[k_system].set_entry (iF, iF, 1.0);
	eq_matrixes[k_system].set_f_entry (iF, function_FCondition (k_system, coordinates, area, 0));
	eq_matrixes[k_system].clear_row (iF);

	// set N function
	area = mesh_pointer->get_area (mesh_pointer->get_n_elements () - 1);
	mesh_pointer->get_N_boundaries (coordinates);
	iF = N_functions - 1;

	eq_matrixes[k_system].set_entry (iF, iF, 1.0);
	eq_matrixes[k_system].set_f_entry (iF, function_FCondition (k_system, coordinates, area, 1));
	eq_matrixes[k_system].clear_row (iF);

	// for hermitian set functions derivative
	if (mesh_pointer->get_amount_non_zero_functions (0) == 4)
	{
		// set 1 function
		int iF = 1;
		int area = mesh_pointer->get_area (area_element);
		double coordinates[1];
		mesh_pointer->get_0_boundaries (coordinates);
		
		eq_matrixes[k_system].set_entry (iF, iF, 1.0);
		eq_matrixes[k_system].set_f_entry (iF, function_FCondition_der (k_system, coordinates, area, 0));
		eq_matrixes[k_system].clear_row (iF);

		// set N - 1, N - 2 function
		area = mesh_pointer->get_area (mesh_pointer->get_n_elements () - 1);
		mesh_pointer->get_N_boundaries (coordinates);
		iF = N_functions - 1;
		eq_matrixes[k_system].set_f_entry (iF, function_FCondition_der (k_system, coordinates, area, 1));

		iF = N_functions - 2;
		eq_matrixes[k_system].set_entry (iF, iF, 1.0);
		eq_matrixes[k_system].set_f_entry (iF, function_FCondition (k_system, coordinates, area, 1));
		eq_matrixes[k_system].clear_row (iF);
	}
}

template<class Type_Mesh>
double Task_Test_1D<Type_Mesh>::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	double x = coordinates[0];
	//return x;
	return x * x * x + x * x + x;
}

template<class Type_Mesh>
double Task_Test_1D<Type_Mesh>::function_f (int k_system, double * coordinates, int area)
{
	double x = coordinates[0];
	//return x;
	return x * x * x + x * x + x - 6.0 * x - 2.0;
}

template<class Type_Mesh>
double Task_Test_1D<Type_Mesh>::function_FCondition_der (int k_system, double * coordinates, int area, int boundary)
{
	double x = coordinates[0];
	//return 1.0;
	return 3.0 * x * x + 2.0 * x + 1.0;
}

template<class Type_Mesh>
void Task_Test_1D<Type_Mesh>::save_slice (int k_system, char * file_sour, char * file_dest, bool keep_points)
{ 
	int dim_node = mesh_pointer->get_dimentionality ();
	double * coordinates = new double[dim_node];
	double * coordinatesN = new double[dim_node];
	double value; // value of the solution
	bool exists; // flag of existence of solution's value 
	int n_amount = 0; // amount of points written in file 
	double h; // iteration step
	int iter_amount; // amount of points to iterate through
	double vstart, vend; // starting and ending values of variable 

	FILE * points = fopen (file_sour, "r");
	fscanf (points, "%lf", &h); // get step
	fclose (points);

	// get starting value
	mesh_pointer->get_0_boundaries (coordinatesN);
	vstart = coordinatesN[0];
	// get ending value
	mesh_pointer->get_N_boundaries (coordinatesN);
	vend = coordinatesN[0];
	// get amount of points to iterate through
	iter_amount = (int)((vend - vstart) / h);
	
	// open file for result input
	FILE * result_points = fopen (file_dest, "w");
	for (int i = 0; i < iter_amount + 1; i++)
	{
		coordinates[0] = vstart + i * h;
		//calculate the point
		value = 0.0;
		// get value in the point
		exists = get_solution_in_point (k_system, coordinates, &value);
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

	delete[] coordinates;
	delete[] coordinatesN;
}
