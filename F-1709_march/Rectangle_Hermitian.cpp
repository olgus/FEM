#include "Rectangle_Hermitian.h"

Rectangle_Element_Hermitian::Rectangle_Element_Hermitian ()
{
	n_functions = 16;
	D = NULL;
	order = 3;
}

Rectangle_Element_Hermitian::Rectangle_Element_Hermitian (const Rectangle_Element_Hermitian & rect3)
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

	// copy n_functions
	n_functions = rect3.n_functions;

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

	// set D if it exists
	if (rect3.D != NULL)
		D->Full_Copy (*rect3.D);
	else
		D = NULL;
}

Rectangle_Element_Hermitian::~Rectangle_Element_Hermitian ()
{
	if (D != NULL)
		delete D;
}

void Rectangle_Element_Hermitian::set_global_functions ()
{
	if (global_functions != NULL)
		delete[] global_functions;
	global_functions = new int [n_functions];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			global_functions[i * 4 + j] = defining_nodes[i] * 4 + j;
		}
	}
}

bool Rectangle_Element_Hermitian::point_inside (const Mesh_Prototype & mesh, double * coordinates)
{
	// count point if in fits into [)
	// an angle point will only count for an element where its left down node
	// if there is no such element, the point will be on the boudary, so look down
	// TEST -ZERO_Rectangle_Element_hermitian
	if (coordinates[0] >= x0 - ZERO_Rectangle_Element_hermitian && coordinates[0] < xN)
	{
		if (coordinates[1] >= y0 - ZERO_Rectangle_Element_hermitian && coordinates[1] < yN)
		{
			return true;
		}
	}

	double cN[2];
	mesh.get_N_boundaries (cN);
	// or ], if it is on the right boundary
	// check if the element is also boundary one
	// fits into x, but the element is on the boundary and value too
	if (coordinates[0] >= x0 - ZERO_Rectangle_Element_hermitian && coordinates[0] < xN)
	{
		if (fabs(coordinates[1] - cN[1]) < ZERO_Rectangle_Element_hermitian && fabs (yN - cN[1]) < ZERO_Rectangle_Element_hermitian)
		{
			return true;
		}
	}
	// fits into y, but the element is on the boundary and value too
	if (coordinates[1] >= y0 - ZERO_Rectangle_Element_hermitian && coordinates[1] < yN)
	{
		if (fabs (coordinates[0] - cN[0]) < ZERO_Rectangle_Element_hermitian && fabs (xN - cN[0]) < ZERO_Rectangle_Element_hermitian)
		{
			return true;
		}
	}
	// right upper point
	if (fabs (xN - cN[0]) < ZERO_Rectangle_Element_hermitian && fabs (yN - cN[1]) < ZERO_Rectangle_Element_hermitian && fabs (coordinates[0] - cN[0]) < ZERO_Rectangle_Element_hermitian && fabs (coordinates[1] - cN[1]) < ZERO_Rectangle_Element_hermitian)
	{
		return true;
	}

	return false;
}

void Rectangle_Element_Hermitian::get_D_local_matrix (Matrix * D_matrix)
{
	D_matrix->Copy (*D);
}

double Rectangle_Element_Hermitian::get_basis_function_value (int k_func, double * coordinates)
{
	double x = coordinates[0];
	double y = coordinates[1];
	double hx = xN - x0;
	double hy = yN - y0;

	double xi, etta;
	xi = (x - x0) / hx;
	etta = (y - y0) / hy;

	double xF[4];
	xF[0] = 1.0 - 3.0 * pow (xi, 2.0) + 2.0 * pow (xi, 3.0);
	xF[1] = hx * (xi - 2.0 * pow (xi, 2.0) + pow (xi, 3.0));
	xF[2] = 3.0 * pow (xi, 2.0) - 2.0 * pow (xi, 3.0);
	xF[3] = hx * (- pow (xi, 2.0) + pow (xi, 3.0));
	double yF[4];
	yF[0] = 1.0 - 3.0 * pow (etta, 2.0) + 2.0 * pow (etta, 3.0);
	yF[1] = hy * (etta - 2.0 * pow (etta, 2.0) + pow (etta, 3.0));
	yF[2] = 3.0 * pow (etta, 2.0) - 2.0 * pow (etta, 3.0);
	yF[3] = hy * (-pow (etta, 2.0) + pow (etta, 3.0));

	int m, n;
	m = 2 * ((k_func / 4) % 2) + k_func % 2;
	n = 2 * (k_func / 8) + ((k_func / 2) % 2);

	return xF[m] * yF[n];
}

double Rectangle_Element_Hermitian::get_basis_function_derivative (int k_func, int k_var, double * coordinates)
{
	double r = 0.0;
	double x = coordinates[0];
	double y = coordinates[1];
	double hx = xN - x0;
	double hy = yN - y0;
	double xi, etta;
	xi = (x - x0) / hx;
	etta = (y - y0) / hy;

	double xF[4];
	xF[0] = 1.0 - 3.0 * pow (xi, 2.0) + 2.0 * pow (xi, 3.0);
	xF[1] = hx * (xi - 2.0 * pow (xi, 2.0) + pow (xi, 3.0));
	xF[2] = 3.0 * pow (xi, 2.0) - 2.0 * pow (xi, 3.0);
	xF[3] = hx * (-pow (xi, 2.0) + pow (xi, 3.0));
	double yF[4];
	yF[0] = 1.0 - 3.0 * pow (etta, 2.0) + 2.0 * pow (etta, 3.0);
	yF[1] = hy * (etta - 2.0 * pow (etta, 2.0) + pow (etta, 3.0));
	yF[2] = 3.0 * pow (etta, 2.0) - 2.0 * pow (etta, 3.0);
	yF[3] = hy * (-pow (etta, 2.0) + pow (etta, 3.0));
	double xderF[4];
	xderF[0] = - 6.0 * xi + 6.0 * pow (xi, 2.0);
	xderF[1] = hx * (1.0 - 4.0 * xi + 3.0 * pow (xi, 2.0));
	xderF[2] = 6.0 * xi - 6.0 * pow (xi, 2.0);
	xderF[3] = hx * (-2.0 * xi + 3.0 * pow (xi, 2.0));
	double yderF[4];
	yderF[0] = - 6.0 * etta + 6.0 * pow (etta, 2.0);
	yderF[1] = hy * (1.0 - 4.0 * etta + 3.0 * pow (etta, 2.0));
	yderF[2] = 6.0 * etta - 6.0 * pow (etta, 2.0);
	yderF[3] = hy * (-2.0 * etta + 3.0 * pow (etta, 2.0));

	int m, n;
	m = 2 * ((k_func / 4) % 2) + k_func % 2;
	n = 2 * (k_func / 8) + ((k_func / 2) % 2);

	for (int i = 0; i < 4; i++)
	{
		xderF[i] /= hx;
		yderF[i] /= hy;
	}

	r = 1.0;
	switch (k_var)
	{
	case 0:
		r *= xderF[m];
		r *= yF[n];
		break;
	case 1:
		r *= yderF[n];
		r *= xF[m];
		break;
	case 10: // partial derivative (?)
		r *= xderF[m];
		r *= yderF[n];
		break;
	default:
		break;
	}
	return r; 
}

double Rectangle_Element_Hermitian::get_basis_function_derivative_second (int k_func, int k_var, double * coordinates)
{
	double r;
	double x = coordinates[0];
	double y = coordinates[1];
	double hx = xN - x0;
	double hy = yN - y0;
	double xi, etta;
	xi = (x - x0) / hx;
	etta = (y - y0) / hy;

	double xF[4];
	xF[0] = 1.0 - 3.0 * pow (xi, 2.0) + 2.0 * pow (xi, 3.0);
	xF[1] = hx * (xi - 2.0 * pow (xi, 2.0) + pow (xi, 3.0));
	xF[2] = 3.0 * pow (xi, 2.0) - 2.0 * pow (xi, 3.0);
	xF[3] = hx * (-pow (xi, 2.0) + pow (xi, 3.0));
	double yF[4];
	yF[0] = 1.0 - 3.0 * pow (etta, 2.0) + 2.0 * pow (etta, 3.0);
	yF[1] = hy * (etta - 2.0 * pow (etta, 2.0) + pow (etta, 3.0));
	yF[2] = 3.0 * pow (etta, 2.0) - 2.0 * pow (etta, 3.0);
	yF[3] = hy * (-pow (etta, 2.0) + pow (etta, 3.0));
	double xderF[4];
	xderF[0] = -6.0 + 12.0 * xi;
	xderF[1] = hx * (- 4.0 + 6.0 * xi);
	xderF[2] = 6.0 - 12.0 * xi;
	xderF[3] = hx * (-2.0 + 6.0 * xi);
	double yderF[4];
	yderF[0] = -6.0 + 12.0 * etta;
	yderF[1] = hy * (- 4.0 + 6.0 * etta);
	yderF[2] = 6.0 - 12.0 * etta;
	yderF[3] = hy * (-2.0 + 6.0 * etta);

	int m, n;
	m = 2 * ((k_func / 4) % 2) + k_func % 2;
	n = 2 * (k_func / 8) + ((k_func / 2) % 2);

	for (int i = 0; i < 4; i++)
	{
		xderF[i] /= hx;
		xderF[i] /= hx;
		yderF[i] /= hy;
		yderF[i] /= hy;
	}

	r = 1.0;
	switch (k_var)
	{
	case 0:
		r *= xderF[m];
		r *= yF[n];
		break;
	case 1:
		r *= yderF[n];
		r *= xF[m];
		break;
	default:
		break;
	}
	return r;
}

double Rectangle_Element_Hermitian::get_function_value (int i, int j, int * m, double * coordinates)
{
	double r = 0.0;

	switch (m[0])
	{
	case 0: // g(i,j)
		r += get_basis_function_derivative (i, 0, coordinates) * get_basis_function_derivative (j, 0, coordinates);
		r += get_basis_function_derivative (i, 1, coordinates) * get_basis_function_derivative (j, 1, coordinates);
		break;
	case 1: // m(i,j)
		r += get_basis_function_value (i, coordinates) * get_basis_function_value (j, coordinates);
		break;
	case 2: // d(i,j)
		r += get_basis_function_derivative_second (i, 0, coordinates) * get_basis_function_derivative_second (j, 0, coordinates);
		r += get_basis_function_derivative_second (i, 1, coordinates) * get_basis_function_derivative_second (j, 1, coordinates);
		break;
	default:
		break;
	}
	return r;
}

int Rectangle_Element_Hermitian::get_amount_non_zero_functions ()
{
	return n_functions; 
}

void Rectangle_Element_Hermitian::inside_prepare (const Mesh_Prototype & mesh)
{
	D = new Matrix (n_functions, n_functions);

	double ** c;
	c = new double*[n_base_nodes];
	for (int i = 0; i < n_base_nodes; i++)
	{
		c[i] = new double[2];
	}
	for (int i = 0; i < n_base_nodes; i++)
	{
		mesh.get_node_coordinates (base_nodes[i], c[i]);
	}
	x0 = c[0][0];
	xN = c[1][0];
	y0 = c[0][1];
	yN = c[2][1];

	int m[] = {2};
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
