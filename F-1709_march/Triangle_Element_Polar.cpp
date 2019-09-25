#include "Triangle_Element_Polar.h"

Triangle_Polar::Triangle_Polar ()
{
}

Triangle_Polar::Triangle_Polar (const Triangle_Polar & triangle)
{
	if (n_def_nodes != triangle.n_def_nodes) // if triangle's amount of nodes doesn't equal to amount of nodes of current triangle
	{
		if (defining_nodes != NULL) // if memory was even allocated
			delete[] defining_nodes; // free it

		n_def_nodes = triangle.n_def_nodes; // set new amount of nodes
		defining_nodes = new int[n_def_nodes]; // allocate the memory

		if (base_nodes != NULL) // if memory was even allocated
			delete[] base_nodes; // free it

		n_base_nodes = triangle.n_base_nodes; // set new amount of nodes
		base_nodes = new int[n_base_nodes]; // allocate the memory
	}

	if (defining_nodes == NULL)
		defining_nodes = new int[n_def_nodes]; // allocate the memory

	for (int i = 0; i < n_def_nodes; i++) // copy nodes
		defining_nodes[i] = triangle.defining_nodes[i];

	if (base_nodes == NULL)
		base_nodes = new int[n_base_nodes]; // allocate the memory

	for (int i = 0; i < n_base_nodes; i++) // copy nodes
		base_nodes[i] = triangle.base_nodes[i];

	area = triangle.area;

	if (triangle.Alpha != NULL)
		Alpha->Full_Copy (*triangle.Alpha);
	else
		Alpha = NULL;
}

Triangle_Polar::~Triangle_Polar ()
{
}

double Triangle_Polar::get_function_value (int i, int j, int * m, double * coordinates)
{
	double val = 0.0;
	double r = coordinates[0];
	double phi = coordinates[1];

	switch (m[0])
	{
	case 0: // g(i,j)
		val += get_basis_function_derivative (i, 0, coordinates) * get_basis_function_derivative (j, 0, coordinates);
		val += (get_basis_function_derivative (i, 1, coordinates) * get_basis_function_derivative (j, 1, coordinates)) *
			   (1.0 / pow (r, 2.0));
		val *= r;
		break;
	case 1: // m(i,j)
		val += get_basis_function_value (i, coordinates) * get_basis_function_value (j, coordinates) * r;
		break;
	default:
		val = 0.0;
	}
	return val;
}
