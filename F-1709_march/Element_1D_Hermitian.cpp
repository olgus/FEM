#include "Element_1D_Hermitian.h"

Element_1D_Hermitian::Element_1D_Hermitian ()
{
	n_def_nodes = 2;
	n_base_nodes = 2;
	n_functions = 4;
	defining_nodes = NULL;
	base_nodes = NULL;
	G = NULL;
	M = NULL;
	D = NULL;

	x0 = xN = 0.0;

	dim = 1;
	order = 3;
}

Element_1D_Hermitian::Element_1D_Hermitian (const Element_1D_Hermitian & e)
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

Element_1D_Hermitian::~Element_1D_Hermitian ()
{
}

void Element_1D_Hermitian::get_D_local_matrix (Matrix * D_matrix)
{
	D_matrix->Copy (*D);
}

void Element_1D_Hermitian::inside_prepare (const Mesh_Prototype & mesh)
{
	D = new Matrix (n_functions, n_functions);

	double ** c;
	c = new double*[n_base_nodes];
	for (int i = 0; i < n_base_nodes; i++)
	{
		c[i] = new double[1];
	}
	for (int i = 0; i < n_base_nodes; i++)
	{
		mesh.get_node_coordinates (base_nodes[i], c[i]);
	}
	x0 = c[0][0];
	xN = c[1][0];

	int m[] = { 2 };
	for (int j = 0; j < n_functions; j++)
	{
		for (int i = 0; i < n_functions; i++)
		{
			D->setElem (i, j, integrate (i, j, m));
		}
	}

	for (int i = 0; i < n_base_nodes; i++)
	{
		delete[] c[i];
	}
	delete[] c;
}

double Element_1D_Hermitian::get_basis_function_value (int k_func, double * coordinates)
{
	double h = xN - x0;
	double psi = (coordinates[0] - x0) / h;
	double r = 0.0;

	switch (k_func)
	{
	case 0:
		r = 1.0 - 3.0 * pow (psi, 2.0) + 2.0 * pow (psi, 3.0);
		break;
	case 1:
		r = h * (psi - 2.0 * pow (psi, 2.0) + pow (psi, 3.0));
		break;
	case 2:
		r = 3.0 * pow (psi, 2.0) - 2.0 * pow (psi, 3.0);
		break;
	case 3:
		r = h * (-pow (psi, 2.0) + pow (psi, 3.0));
		break;
	}
	return r;
}

double Element_1D_Hermitian::get_basis_function_derivative (int k_func, int k_var, double * coordinates)
{
	double h = xN - x0;
	double psi = (coordinates[0] - x0) / h;
	double r = 0.0;

	switch (k_func)
	{
	case 0:
		r = -6.0 * psi + 6.0 * pow (psi, 2.0);
		break;
	case 1:
		r = h * (1.0 - 4.0 * psi + 3.0 * pow (psi, 2.0));
		break;
	case 2:
		r = 6.0 * psi - 6.0 * pow (psi, 2.0);
		break;
	case 3:
		r = h * (-2.0 * psi + 3.0 * pow (psi, 2.0));
		break;
	}
	r /= h;
	return r;
}
