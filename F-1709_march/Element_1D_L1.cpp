#include "Element_1D_L1.h"

Element_1D_L1::Element_1D_L1 ()
{
	n_def_nodes = 2;
	n_base_nodes = 2;
	n_functions = 2;
	defining_nodes = NULL;
	base_nodes = NULL;
	G = NULL;
	M = NULL;
	x0 = xN = 0.0;

	dim = 1;
	order = 1;
}

Element_1D_L1::Element_1D_L1 (const Element_1D_L1 & e)
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

Element_1D_L1::~Element_1D_L1 ()
{

}

double Element_1D_L1::get_geometrical_area ()
{
	return xN - x0;
}

bool Element_1D_L1::point_inside (const Mesh_Prototype & mesh, double * coordinates)
{
	if (coordinates[0] >= x0 - ZERO_element && coordinates[0] < xN)
	{
		return true;
	}

	double cN[1];
	mesh.get_N_boundaries (cN);
	if (fabs (xN - cN[0]) < ZERO_element && fabs (coordinates[0] - cN[0]) < ZERO_element)
	{
		return true;
	}

	return false;
}

void Element_1D_L1::inside_prepare (const Mesh_Prototype & mesh)
{
	double * c;
	c = new double[n_base_nodes];
	for (int i = 0; i < n_base_nodes; i++)
	{
		mesh.get_node_coordinates (base_nodes[i], &(c[i]));
	}
	x0 = c[0];
	xN = c[1];

	delete[] c;
}

double Element_1D_L1::integrate (int i, int j, int * m)
{
	// Gauss-4 
	int n_int_nodes = 4;
	// set w(i)
	MathVector * w = new MathVector (n_int_nodes);
	w->setElem (0, (18.0 + pow (30, 0.5)) / 36.0);
	w->setElem (1, (18.0 + pow (30, 0.5)) / 36.0);
	w->setElem (2, (18.0 - pow (30, 0.5)) / 36.0);
	w->setElem (3, (18.0 - pow (30, 0.5)) / 36.0);

	// set master point's cordinates
	MathVector * xi = new MathVector (n_int_nodes);
	xi->setElem (0, -pow (3.0 / 7.0 - 2.0 * pow (6.0 / 5.0, 0.5) / 7.0, 0.5));
	xi->setElem (1, pow (3.0 / 7.0 - 2.0 * pow (6.0 / 5.0, 0.5) / 7.0, 0.5));
	xi->setElem (2, -pow (3.0 / 7.0 + 2.0 * pow (6.0 / 5.0, 0.5) / 7.0, 0.5));
	xi->setElem (3, pow (3.0 / 7.0 + 2.0 * pow (6.0 / 5.0, 0.5) / 7.0, 0.5));

	double coordinates_d[1];
	double r = 0.0;

	MathVector * coordinates = new MathVector (1);
	for (int k = 0; k < n_int_nodes; k++)
	{
		// get coordinates through master coordinates
		coordinates_d[0] = (xN - x0) * (xi->getElem (k) + 1.0) / 2.0 + x0;
		// calculate addition to integral value
		r += w->getElem (k) * get_function_value (i, j, m, coordinates_d);
	}
	// jacobian
	r *= (xN - x0) / 2.0;

	delete w;
	delete xi;
	delete coordinates;
	return r;
}

double Element_1D_L1::get_basis_function_value (int k_func, double * coordinates)
{
	double psi = (coordinates[0] - x0) / (xN - x0);
	double r = 0.0;

	switch (k_func)
	{
	case 0:
		r = 1.0 - psi;
		break;
	case 1: 
		r = psi;
		break;
	}
	return r;
}

double Element_1D_L1::get_basis_function_derivative (int k_func, int k_var, double * coordinates)
{
	double h = xN - x0;
	double r = 0.0;

	switch (k_func)
	{
	case 0:
		r = -1.0;
		break;
	case 1: 
		r = 1.0;
		break;
	}
	r /= h;
	return r;
}

int Element_1D_L1::get_isoline_points (const Mesh_Prototype & mesh, double value, double * q, double * c1, double * c2)
{
	double cn0[2];
	double cn1[2];
	mesh.get_node_coordinates (defining_nodes[0], cn0);
	mesh.get_node_coordinates (defining_nodes[1], cn1);

	c1[0] = ((xN - x0) * value - xN * q[0] + q[1] * x0) / (q[1] - q[0]);
	return 1;
}

int Element_1D_L1::amount_of_integration_points ()
{
	return 4;
}

void Element_1D_L1::integration_points (double ** points, double * weigths, double * jac)
{
	// Gauss-4 
	int n_int_nodes = 4;
	// set w(i)
	weigths[0] = (18.0 + pow (30, 0.5)) / 36.0;
	weigths[1] = (18.0 + pow (30, 0.5)) / 36.0;
	weigths[2] = (18.0 - pow (30, 0.5)) / 36.0;
	weigths[3] = (18.0 - pow (30, 0.5)) / 36.0;

	// set master point's cordinates
	MathVector * xi = new MathVector (n_int_nodes);
	xi->setElem (0, -pow (3.0 / 7.0 - 2.0 * pow (6.0 / 5.0, 0.5) / 7.0, 0.5));
	xi->setElem (1, pow (3.0 / 7.0 - 2.0 * pow (6.0 / 5.0, 0.5) / 7.0, 0.5));
	xi->setElem (2, -pow (3.0 / 7.0 + 2.0 * pow (6.0 / 5.0, 0.5) / 7.0, 0.5));
	xi->setElem (3, pow (3.0 / 7.0 + 2.0 * pow (6.0 / 5.0, 0.5) / 7.0, 0.5));

	for (int k = 0; k < n_int_nodes; k++)
	{
		// get coordinates through master coordinates
		// set points
		points[k][0] = (xN - x0) * (xi->getElem (k) + 1.0) / 2.0 + x0;
	}
	// jacobian
	*jac = (xN - x0) / 2.0;

	delete xi;
}

