#include "Rectangle_Element.h"

Rectangle_Element::Rectangle_Element ()
{
	n_def_nodes = 4;
	n_base_nodes = 4; 
	n_functions = 4;
	defining_nodes = NULL;
	base_nodes = NULL;
	G = NULL;
	M = NULL;
	x0 = xN = y0 = yN = 0.0;

	dim = 2;
}

Rectangle_Element::Rectangle_Element (const Rectangle_Element & Rectangle_Element)
{
	if (n_def_nodes != Rectangle_Element.n_def_nodes) // if triangle's amount of nodes doesn't equal to amount of nodes of current triangle
	{
		if (defining_nodes != NULL) // if memory was even allocated
			delete[] defining_nodes; // free it

		n_def_nodes = Rectangle_Element.n_def_nodes; // set new amount of nodes

		if (base_nodes != NULL) // if memory was even allocated
			delete[] base_nodes; // free it

		n_base_nodes = Rectangle_Element.n_base_nodes; // set new amount of nodes
	}

	if (Rectangle_Element.defining_nodes != NULL)
	{
		defining_nodes = new int[n_def_nodes]; // allocate the memory

		for (int i = 0; i < n_def_nodes; i++) // copy nodes
			defining_nodes[i] = Rectangle_Element.defining_nodes[i];
	}

	if (Rectangle_Element.base_nodes != NULL)
	{
		base_nodes = new int[n_base_nodes]; // allocate the memory

		for (int i = 0; i < n_base_nodes; i++) // copy nodes
			base_nodes[i] = Rectangle_Element.base_nodes[i];
	}

	// set area
	area = Rectangle_Element.area;
	dim = Rectangle_Element.dim;

	// set G if it exists
	if (Rectangle_Element.G != NULL)
		G->Full_Copy (*Rectangle_Element.G);
	else
		G = NULL;

	// set M if it exists
	if (Rectangle_Element.M != NULL)
		M->Full_Copy (*Rectangle_Element.M);
	else
		M = NULL;

	x0 = Rectangle_Element.x0;
	xN = Rectangle_Element.xN;
	y0 = Rectangle_Element.y0;
	yN = Rectangle_Element.yN;
}

Rectangle_Element::~Rectangle_Element ()
{
}

double Rectangle_Element::get_geometrical_area ()
{
	return (xN - x0) * (yN - y0);
}


bool Rectangle_Element::point_inside (const Mesh_Prototype & mesh, double * coordinates)
{
	if (coordinates[0] >= x0 && coordinates[0] <= xN)
	{
		if (coordinates[1] >= y0 && coordinates[1] <= yN)
		{
			return true;
		}
	}
	return false;
}

double Rectangle_Element::integrate (int i, int j, int * m)
{
	// 7 
	int n_int_nodes = 12;
	double a, b, c;
	double wa, wb, wc;
	a = pow ((114.0 - 3.0 * pow (583, 0.5)) / 287.0, 0.5);
	b = pow ((114.0 + 3.0 * pow (583, 0.5)) / 287.0, 0.5);
	c = pow (6.0 / 7.0, 0.5);
	wa = 307.0 / 810.0 + 923.0 / (270.0 * pow (583.0, 0.5));
	wb = 307.0 / 810.0 - 923.0 / (270.0 * pow (583.0, 0.5));
	wc = 98.0 / 405.0;

	// set master point's cordinates
	MathVector * xi = new MathVector (n_int_nodes);
	xi->setElem (0, -c);
	xi->setElem (1, c);
	xi->setElem (2, 0.0);
	xi->setElem (3, 0.0);
	xi->setElem (4, -a);
	xi->setElem (5, a);
	xi->setElem (6, -a);
	xi->setElem (7, a);
	xi->setElem (8, -b);
	xi->setElem (9, b);
	xi->setElem (10, -b);
	xi->setElem (11, b);

	MathVector * etta = new MathVector (n_int_nodes);
	etta->setElem (0, 0.0);
	etta->setElem (1, 0.0);
	etta->setElem (2, -c);
	etta->setElem (3, c);
	etta->setElem (4, -a);
	etta->setElem (5, -a);
	etta->setElem (6, a);
	etta->setElem (7, a);
	etta->setElem (8, -b);
	etta->setElem (9, -b);
	etta->setElem (10, b);
	etta->setElem (11, b);

	MathVector * w = new MathVector (n_int_nodes);
	w->setElem (0, wc);
	w->setElem (1, wc);
	w->setElem (2, wc);
	w->setElem (3, wc);
	w->setElem (4, wa);
	w->setElem (5, wa);
	w->setElem (6, wa);
	w->setElem (7, wa);
	w->setElem (8, wb);
	w->setElem (9, wb);
	w->setElem (10, wb);
	w->setElem (11, wb);

	double coordinates_d[2];
	double r = 0.0;
	for (int k = 0; k < n_int_nodes; k++)
	{
		coordinates_d[0] = (xN - x0) * (xi->getElem(k) + 1.0) / 2.0 + x0;
		coordinates_d[1] = (yN - y0) * (etta->getElem (k) + 1.0) / 2.0 + y0;
		r += get_function_value (i, j, m, coordinates_d) * w->getElem (k);
	}
	double jac = (xN - x0) * (yN - y0) / 4.0;
	r *= jac;

	delete xi;
	delete etta;
	delete w;
	return r;
}

void Rectangle_Element::inside_prepare (const Mesh_Prototype & mesh)
{
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

	for (int i = 0; i < n_base_nodes; i++)
	{
		delete[] c[i];
	}
	delete[] c;
}

double Rectangle_Element::get_basis_function_value (int k_func, double * coordinates)
{
	double x = coordinates[0];
	double y = coordinates[1];
	double hx = xN - x0;
	double hy = yN - y0;

	double xF[2];
	xF[0] = (xN - x) / hx;
	xF[1] = (x - x0) / hx;
	double yF[2];
	yF[0] = (yN - y) / hy;
	yF[1] = (y - y0) / hy;

	int m = k_func % 2; 
	int n = k_func / 2;
	
	return xF[m] * yF[n];
}

double Rectangle_Element::get_basis_function_derivative (int k_func, int k_var, double * coordinates)
{
	double x = coordinates[0];
	double y = coordinates[1];
	double hx = xN - x0;
	double hy = yN - y0;
	double xF[2];
	xF[0] = (xN - x) / hx;
	xF[1] = (x - x0) / hx;
	double yF[2];
	yF[0] = (yN - y) / hy;
	yF[1] = (y - y0) / hy;

	int m = k_func % 2;
	int n = k_func / 2;

	double r = 1.0;
	switch (k_var)
	{
	case 0:
		r /= hx;
		r *= pow (-1.0, m + 1);
		r *= yF[n];
		break;
	case 1:
		r /= hy;
		r *= pow (-1.0, n + 1);
		r *= xF[m];
		break;
	}
	return r;
}

double Rectangle_Element::get_function_value (int i, int j, int * m, double * coordinates)
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
	default:
		break;
	}
	return r;
}

int Rectangle_Element::get_amount_second_condition ()
{
	return 2;
}

void Rectangle_Element::get_edge_functions (const Mesh_Prototype & mesh, int n1, int n2, int * functions)
{
	int n1_local;
	int n2_local;
	int n;

	for (int i = 0; i < 4; i++)
	{
		if (n1 == base_nodes[i])
			n1_local = i;
		if (n2 == base_nodes[i])
			n2_local = i;
	}
	if (n1_local > n2_local)
	{
		n = n1_local;
		n1_local = n2_local;
		n2_local = n;
	}

	functions[0] = n1_local;
	functions[1] = n2_local;
}

bool Rectangle_Element::edge_exists (int n1, int n2)
{
	int n1_local = -1;
	int n2_local = -1;
	int n;

	for (int i = 0; i < 4; i++)
	{
		if (n1 == base_nodes[i])
			n1_local = i;
		if (n2 == base_nodes[i])
			n2_local = i;
	}
	if (n1_local == -1 || n2_local == -1)
		return false;
	if (n1_local > n2_local)
	{
		n = n1_local;
		n1_local = n2_local;
		n2_local = n;
	}
	if ((n2_local == 0 && n2_local == 3) || (n2_local == 1 && n2_local == 2))
		return false;
	return true;
}

int Rectangle_Element::get_amount_edges ()
{
	return 4;
}

int Rectangle_Element::amount_of_integration_points ()
{
	return 12;
}

void Rectangle_Element::integration_points (double ** points, double * weigths, double * jac)
{
	// 7 
	int n_int_nodes = 12;
	double a, b, c;
	double wa, wb, wc;
	a = pow ((114.0 - 3.0 * pow (583, 0.5)) / 287.0, 0.5);
	b = pow ((114.0 + 3.0 * pow (583, 0.5)) / 287.0, 0.5);
	c = pow (6.0 / 7.0, 0.5);
	wa = 307.0 / 810.0 + 923.0 / (270.0 * pow (583.0, 0.5));
	wb = 307.0 / 810.0 - 923.0 / (270.0 * pow (583.0, 0.5));
	wc = 98.0 / 405.0;

	// set master point's cordinates
	MathVector * xi = new MathVector (n_int_nodes);
	xi->setElem (0, -c);
	xi->setElem (1, c);
	xi->setElem (2, 0.0);
	xi->setElem (3, 0.0);
	xi->setElem (4, -a);
	xi->setElem (5, a);
	xi->setElem (6, -a);
	xi->setElem (7, a);
	xi->setElem (8, -b);
	xi->setElem (9, b);
	xi->setElem (10, -b);
	xi->setElem (11, b);

	MathVector * etta = new MathVector (n_int_nodes);
	etta->setElem (0, 0.0);
	etta->setElem (1, 0.0);
	etta->setElem (2, -c);
	etta->setElem (3, c);
	etta->setElem (4, -a);
	etta->setElem (5, -a);
	etta->setElem (6, a);
	etta->setElem (7, a);
	etta->setElem (8, -b);
	etta->setElem (9, -b);
	etta->setElem (10, b);
	etta->setElem (11, b);

	weigths[0] = wc;
	weigths[1] = wc;
	weigths[2] = wc;
	weigths[3] = wc;
	weigths[4] = wa;
	weigths[5] = wa;
	weigths[6] = wa;
	weigths[7] = wa;
	weigths[8] = wb;
	weigths[9] = wb;
	weigths[10] = wb;
	weigths[11] = wb;

	for (int k = 0; k < n_int_nodes; k++)
	{
		points[k][0] = (xN - x0) * (xi->getElem (k) + 1.0) / 2.0 + x0;
		points[k][1] = (yN - y0) * (etta->getElem (k) + 1.0) / 2.0 + y0;
	}
	*jac = (xN - x0) * (yN - y0) / 4.0;

	delete xi;
	delete etta;
}

int Rectangle_Element::edge_amount_of_integration_points ()
{
	return 4;
}

void Rectangle_Element::edge_integration_points (const Mesh_Prototype & mesh, int n1, int n2, double ** points, double * weigths, double * jac)
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

	double coordinates_1[2];
	double coordinates_2[2];
	double r = 0.0;
	mesh.get_node_coordinates (n1, coordinates_1);
	mesh.get_node_coordinates (n2, coordinates_2);

	double L = pow (pow (coordinates_2[0] - coordinates_1[0], 2.0) + pow (coordinates_2[1] - coordinates_1[1], 2.0), 0.5);
	double l;
	for (int k = 0; k < n_int_nodes; k++)
	{
		// get coordinates through master coordinates
		l = xi->getElem (k) + 1.0;
		// save points
		points[k][0] = (coordinates_2[0] - coordinates_1[0]) * l / 2.0 + coordinates_1[0];
		points[k][1] = (coordinates_2[1] - coordinates_1[1]) * l / 2.0 + coordinates_1[1];
	}
	// jacobian (half the edge length?)
	*jac = L / 2.0;

	delete xi;
}