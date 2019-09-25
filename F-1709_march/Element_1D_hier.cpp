#include "Element_1D_hier.h"

Element_1D_hier::Element_1D_hier ()
{
	n_def_nodes = 2;
	n_base_nodes = 2;
	n_functions = 3;
	defining_nodes = NULL;
	base_nodes = NULL;
	G = NULL;
	M = NULL;
	x0 = xN = 0.0;

	dim = 1;
	order = 2;
}

Element_1D_hier::Element_1D_hier (const Element_1D_hier & e)
{
	if (n_def_nodes != e.n_def_nodes) // if triangle's amount of nodes doesn't equal to amount of nodes of current triangle
	{
		if (defining_nodes != NULL) // if memory was even allocated
			delete[] defining_nodes; // free it

		n_def_nodes = e.n_def_nodes; // set new amount of nodes

		if (base_nodes != NULL) // if memory was even allocated
			delete[] base_nodes; // free it

		n_base_nodes = e.n_base_nodes; // set new amount of nodes
	}

	if (e.defining_nodes != NULL)
	{
		defining_nodes = new int[n_def_nodes]; // allocate the memory

		for (int i = 0; i < n_def_nodes; i++) // copy nodes
			defining_nodes[i] = e.defining_nodes[i];
	}

	if (e.base_nodes != NULL)
	{
		base_nodes = new int[n_base_nodes]; // allocate the memory

		for (int i = 0; i < n_base_nodes; i++) // copy nodes
			base_nodes[i] = e.base_nodes[i];
	}

	// set area
	area = e.area;
	n_functions = e.n_functions;
	dim = e.dim;

	// copy functions' numbers
	if (e.global_functions != NULL)
	{
		global_functions = new int[n_functions];
		for (int i = 0; i < n_functions; i++) // copy numbers
			global_functions[i] = e.global_functions[i];
	}

	// set G if it exists
	if (e.G != NULL)
		G->Full_Copy (*e.G);
	else
		G = NULL;

	// set M if it exists
	if (e.M != NULL)
		M->Full_Copy (*e.M);
	else
		M = NULL;

	x0 = e.x0;
	xN = e.xN;
}

Element_1D_hier::~Element_1D_hier ()
{
}

double Element_1D_hier::get_basis_function_value (int k_func, double * coordinates)
{
	double psi = (coordinates[0] - x0) / (xN - x0);
	double r = 0.0;

	switch (k_func)
	{
	case 0:
		r = 1.0 - psi;
		break;
	case 1: // second
		r = -4.0 * psi * (psi - 1.0);
		break;
	case 2:
		r = psi;
		break;
	}
	return r;
}

double Element_1D_hier::get_basis_function_derivative (int k_func, int k_var, double * coordinates)
{
	double h = xN - x0;
	double psi = (coordinates[0] - x0) / h;
	double r = 0.0;
	
	switch (k_func)
	{
	case 0:
		r = -1.0; 
		break;
	case 1: // second
		r = -4.0 * (psi - 1.0) - 4.0 * psi;
		break;
	case 2:
		r = 1.0;
		break;
	}
	r /= h;
	return r;
}