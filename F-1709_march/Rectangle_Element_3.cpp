#include "Rectangle_Element_3.h"

Rectangle_Element_3::Rectangle_Element_3 ()
{
	n_def_nodes = 16;
	n_functions = 16;
	order = 3;
}

Rectangle_Element_3::Rectangle_Element_3 (const Rectangle_Element_3 & rect3)
{
	if (n_def_nodes != rect3.n_def_nodes) // if triangle's amount of nodes doesn't equal to amount of nodes of current triangle
	{
		if (defining_nodes != NULL) // if memory was even allocated
			delete[] defining_nodes; // free it

		n_def_nodes = rect3.n_def_nodes; // set new amount of nodes
		defining_nodes = new int[n_def_nodes]; // allocate the memory

		if (base_nodes != NULL) // if memory was even allocated
			delete[] base_nodes; // free it

		n_base_nodes = rect3.n_base_nodes; // set new amount of nodes
		base_nodes = new int[n_base_nodes]; // allocate the memory
	}

	if (defining_nodes == NULL)
		defining_nodes = new int[n_def_nodes]; // allocate the memory

	for (int i = 0; i < n_def_nodes; i++) // copy nodes
		defining_nodes[i] = rect3.defining_nodes[i];

	if (base_nodes == NULL)
		base_nodes = new int[n_base_nodes]; // allocate the memory

	for (int i = 0; i < n_base_nodes; i++) // copy nodes
		base_nodes[i] = rect3.base_nodes[i];

	// set area
	area = rect3.area;

	// set G if it exists
	if (rect3.G != NULL)
		G->Full_Copy (*rect3.G);
	else
		G = NULL;

	// set M if it exists
	if (rect3.M != NULL)
		M->Full_Copy (*rect3.M);
	else
		M = NULL;
}

Rectangle_Element_3::~Rectangle_Element_3 ()
{
}

double Rectangle_Element_3::get_basis_function_value (int k_func, double * coordinates)
{
	double x = coordinates[0];
	double y = coordinates[1];
	double hx = xN - x0;
	double hy = yN - y0;

	double xF[4];
	xF[0] = (xN - x) / hx;
	xF[1] = (x - x0) / hx;
	xF[2] = (x - x0) / hx;
	xF[3] = (x - x0) / hx;
	double yF[4];
	yF[0] = (yN - y) / hy;
	yF[1] = (y - y0) / hy;
	yF[2] = (y - y0) / hy;
	yF[3] = (y - y0) / hy;

	int m = k_func % 4;
	int n = k_func / 4;

	return xF[m] * yF[n];
}

double Rectangle_Element_3::get_basis_function_derivative (int k_func, int k_var, double * coordinates)
{
	return 0.0;
}
