#include "Triangle_Element.h"

Triangle::Triangle ()
{
	n_def_nodes = 3;
	n_base_nodes = 3;
	defining_nodes = NULL;
	base_nodes = NULL;
	Alpha = NULL;
	Determinant_M = NULL;
	n_functions = 3;
	determinant = 0.0;

	D = NULL;
	DDx = NULL;
	DDy = NULL;

	dim = 2;
}

Triangle::Triangle (const Triangle & triangle)
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

	n_functions = triangle.n_functions;
	determinant = triangle.determinant;
	area = triangle.area;
	dim = triangle.dim;

	// set Alpha if it exists
	if (triangle.Alpha != NULL)
	{
		if (Alpha != NULL)
			delete Alpha;
		Alpha = new Matrix (triangle.Alpha->Size0 (), triangle.Alpha->Size1 ());
		Alpha->Full_Copy (*triangle.Alpha);
	}
	else
		Alpha = NULL;

	// set G if it exists
	if (triangle.G != NULL)
	{
		if (G != NULL)
			delete G;
		G = new Matrix (triangle.G->Size0 (), triangle.G->Size1 ());
		G->Full_Copy (*triangle.G);
	}
	else
		G = NULL;

	// set M if it exists
	if (triangle.M != NULL)
	{
		if (M != NULL)
			delete M;
		M = new Matrix (triangle.M->Size0 (), triangle.M->Size1 ());
		M->Full_Copy (*triangle.M);
	}
	else
		M = NULL;

	// set Determinant_M if it exists
	if (triangle.Determinant_M != NULL)
	{
		if (Determinant_M != NULL)
			delete Determinant_M;
		Determinant_M = new Matrix (triangle.Determinant_M->Size0 (), triangle.Determinant_M->Size1 ());
		Determinant_M->Full_Copy (*triangle.Determinant_M);
	}
	else
		Determinant_M = NULL;

	// set D if it exists
	if (triangle.D != NULL)
	{
		if (D != NULL)
			delete D;
		D = new Matrix (triangle.D->Size0 (), triangle.D->Size1 ());
		D->Full_Copy (*triangle.D);
	}
	else
		D = NULL;

	if (triangle.DDx != NULL)
	{
		if (DDx != NULL)
			delete DDx;
		DDx = new MathVector (triangle.DDx->getSize ());
		DDx->Copy (*triangle.DDx);
	}
	else
		DDx = NULL;

	if (triangle.DDy != NULL)
	{
		if (DDy != NULL)
			delete DDy;
		DDy = new MathVector (triangle.DDy->getSize ());
		DDy->Copy (*triangle.DDy);
	}
	else
		DDy = NULL;
}

Triangle::~Triangle ()
{
	if (Alpha != NULL)
		delete Alpha;
	Alpha = NULL;
	if (Determinant_M != NULL)
		delete Determinant_M;
	Determinant_M = NULL;

	if (D != NULL)
		delete D;
	D = NULL;

	if (DDx != NULL)
		delete DDx;
	DDx = NULL;

	if (DDy != NULL)
		delete DDy;
	DDy = NULL;
}

double Triangle::get_geometrical_area ()
{
	return fabs (determinant) / 2.0;
}

bool Triangle::point_inside (const Mesh_Prototype & mesh, double * coordinates)
{
	double tr_coord[3][2];
	// get element's coordinates
	for (int i = 0; i < 3; i++)
	{
		mesh.get_node_coordinates (base_nodes[i], tr_coord[i]);
	}

	// get its geometrical area
	double tr_S = get_geometrical_area ();

	// calculate geometrical areas of triangles made with point (coordinates)
	double square = 0.0;
	for (int i = 0; i < 3; i++)
	{
		square += fabs ((tr_coord[(i + 1) % 3][0] - tr_coord[i][0]) * (coordinates[1] - tr_coord[i][1]) -
			(coordinates[0] - tr_coord[i][0]) * (tr_coord[(i + 1) % 3][1] - tr_coord[i][1]));
	}
	square /= 2.0;

	// if half of their sum is equal to element's geometrical area
	if (abs (square - tr_S) < ZERO_element)
	{
		// point is in the element
		return true;
	}
	return false;
}

void Triangle::inside_prepare (const Mesh_Prototype & mesh)
{
	Alpha = new Matrix (3, 3);
	Determinant_M = new Matrix (3, 3);
	double coordinates[2];

	for (int i = 0; i < 3; i++)
	{
		mesh.get_node_coordinates (base_nodes[i], coordinates);
		Determinant_M->setElem (0, i, 1.0);
		Determinant_M->setElem (1, i, coordinates[0]);
		Determinant_M->setElem (2, i, coordinates[1]);
	}

	int n1, n2;
	double coordinates_2[2];
	for (int i = 0; i < 3; i++)
	{
		n1 = i / 2;
		n2 = (i + 3) / 2;
		mesh.get_node_coordinates (base_nodes[n1], coordinates);
		mesh.get_node_coordinates (base_nodes[n2], coordinates_2);
		length[i] = sqrt (pow (coordinates[0] - coordinates_2[0], 2) + pow (coordinates[1] - coordinates_2[1], 2));
	}
	determinant = Determinant_M->get_determinant_3 ();

	Determinant_M->inverse_matrix_size_3 (Alpha);

	// set D
	{
		if (D != NULL)
			delete D;
		D = new Matrix (n_functions, n_functions);
		double r;
		int m[] = { 10 * (0 + 1) };
		for (int i = 0, i_end = D->Size0 (); i < i_end; i++)
		{
			for (int j = 0, j_end = D->Size1 (); j < j_end; j++)
			{
				r = integrate (i, j, m);
				D->setElem (i, j, r);
			}
		}
	}
	// set DDx
	{
		if (DDx != NULL)
			delete D;
		DDx = new MathVector (n_functions * n_functions * n_functions);
		int m[] = { 9, 0, 1, 0 };
		for (int i = 0; i < n_functions; i++)
		{
			for (int j = 0; j < n_functions; j++)
			{
				for (int k = 0; k < n_functions; k++)
				{
					m[1] = k;
					DDx->setElem (i * n_functions * n_functions + j * n_functions + k, integrate (i, j, m));
				}
			}
		}
	}
	// set DDy
	{
		if (DDy != NULL)
			delete D;
		DDy = new MathVector (n_functions * n_functions * n_functions);
		int m[] = { 9, 0, 0, 1 };
		for (int i = 0; i < n_functions; i++)
		{
			for (int j = 0; j < n_functions; j++)
			{
				for (int k = 0; k < n_functions; k++)
				{
					m[1] = k;
					DDy->setElem (i * n_functions * n_functions + j * n_functions + k, integrate (i, j, m));
				}
			}
		}
	}
}

double Triangle::integrate (int i, int j, int * m)
{
	// Gauss-4 
	int n_int_nodes = 4;
	// set w(i)
	MathVector * w = new MathVector(n_int_nodes);
	w->setElem (0, -9.0 / 32.0);
	w->setElem (1, 25.0 / 96.0);
	w->setElem (2, 25.0 / 96.0);
	w->setElem (3, 25.0 / 96.0);

	// set master point's cordinates
	MathVector * xi = new MathVector (n_int_nodes);
	xi->setElem (0, 1.0 / 3.0);
	xi->setElem (1, 3.0 / 5.0);
	xi->setElem (2, 1.0 / 5.0);
	xi->setElem (3, 1.0 / 5.0);

	MathVector * etta = new MathVector (n_int_nodes);
	etta->setElem (0, 1.0 / 3.0);
	etta->setElem (1, 1.0 / 5.0);
	etta->setElem (2, 3.0 / 5.0);
	etta->setElem (3, 1.0 / 5.0);

	double coordinates_d[2];
	double r = 0.0;
	
	MathVector * coordinates = new MathVector (3);
	MathVector * master_coordinates = new MathVector (3);
	for (int k = 0; k < n_int_nodes; k++)
	{
		master_coordinates->setElem (0, 1.0 - xi->getElem (k) - etta->getElem (k));
		master_coordinates->setElem (1, xi->getElem (k));
		master_coordinates->setElem (2, etta->getElem (k));
		// get coordinates through master coordinates
		Determinant_M->MultiplyMatrixByVector (*master_coordinates, coordinates);
		coordinates_d[0] = coordinates->getElem (1);
		coordinates_d[1] = coordinates->getElem (2);
		// calculate addition to integral value
		r += w->getElem(k) * get_function_value (i, j, m, coordinates_d);
	}
	// jacobian. Why absolute value though?
	r *= fabs (determinant);

	delete w;
	delete xi;
	delete etta;
	delete coordinates;
	delete master_coordinates;
	return r;
}

double Triangle::get_basis_function_value (int k_func, double * coordinates)
{
	double r = 0.0;
	r += Alpha->Elem (k_func, 0);
	r += Alpha->Elem (k_func, 1) * coordinates[0];
	r += Alpha->Elem (k_func, 2) * coordinates[1];
	return r;
}

double Triangle::get_basis_function_derivative (int k_func, int k_var, double * coordinates)
{
	double r = 0.0;
	r += Alpha->Elem (k_func, k_var + 1);
	return r;
}

bool Triangle::edge_exists (int n1, int n2)
{
	// for triangle if n1 and n2 are both in base nodes, then the edge exists
	int amount = 0;
	for (int i = 0; i < n_base_nodes; i++)
	{
		if (n1 == base_nodes[i] || n2 == base_nodes[i])
		{
			amount++;
		}
	}
	if (amount == 2)
		return true;
	else
		return false;
}

int Triangle::get_amount_second_condition ()
{
	return 2;
}

void Triangle::get_edge_functions (const Mesh_Prototype & mesh, int n1, int n2, int * functions)
{
	// two functions for simple triangle element
	int node;
	int counter = 0;
	for (int i = 0; i < n_def_nodes; i++)
	{
		node = defining_nodes[i];
		if (node == n1 || node == n2)
		{
			functions[counter] = i;
			counter++;
		}
	}
}

void Triangle::get_edges (int ** edges)
{
	edges[0][0] = base_nodes[0];
	edges[0][1] = base_nodes[1];
	edges[1][0] = base_nodes[0];
	edges[1][1] = base_nodes[2];
	edges[2][0] = base_nodes[1];
	edges[2][1] = base_nodes[2];
}

int Triangle::get_amount_edges ()
{
	return 3;
}

void Triangle::get_D (int var_der, Matrix * D_copy)
{
	D_copy->Copy (*D);
}

void Triangle::get_DDF (int var_der1, int var_der2, MathVector * DDvd)
{
	if ((var_der1 == 1) && (var_der2 == 0))
		DDvd->Copy (*DDx);

	if ((var_der1 == 0) && (var_der2 == 1))
		DDvd->Copy (*DDy);
}

void Triangle::get_DDF (int k, int var_der1, int var_der2, Matrix * D_copy)
{
	// "normal" DDF

	double v;
	int m[] = { 9, k, var_der1, var_der2 };

	for (int i = 0; i < n_functions; i++)
	{
		for (int j = 0; j < n_functions; j++)
		{
			v = integrate (i, j, m);
			D_copy->setElem (i, j, v);
		}
	}
}

int Triangle::amount_of_integration_points ()
{
	return 4;
}

void Triangle::integration_points (double ** points, double * weigths, double * jac)
{
	// Gauss-4 
	int n_int_nodes = 4;
	// set w(i)
	weigths[0] = -9.0 / 32.0;
	weigths[1] = 25.0 / 96.0;
	weigths[2] = 25.0 / 96.0;
	weigths[3] = 25.0 / 96.0;

	// set master point's cordinates
	MathVector * xi = new MathVector (n_int_nodes);
	xi->setElem (0, 1.0 / 3.0);
	xi->setElem (1, 3.0 / 5.0);
	xi->setElem (2, 1.0 / 5.0);
	xi->setElem (3, 1.0 / 5.0);

	MathVector * etta = new MathVector (n_int_nodes);
	etta->setElem (0, 1.0 / 3.0);
	etta->setElem (1, 1.0 / 5.0);
	etta->setElem (2, 3.0 / 5.0);
	etta->setElem (3, 1.0 / 5.0);

	double r = 0.0;
	MathVector * coordinates = new MathVector (3);
	MathVector * master_coordinates = new MathVector (3);
	for (int k = 0; k < n_int_nodes; k++)
	{
		master_coordinates->setElem (0, 1.0 - xi->getElem (k) - etta->getElem (k));
		master_coordinates->setElem (1, xi->getElem (k));
		master_coordinates->setElem (2, etta->getElem (k));
		// get coordinates through master coordinates
		Determinant_M->MultiplyMatrixByVector (*master_coordinates, coordinates);
		points[k][0] = coordinates->getElem (1);
		points[k][1] = coordinates->getElem (2);
	}
	// jacobian. Why absolute value though?
	*jac = fabs (determinant);

	delete xi;
	delete etta;
	delete coordinates;
	delete master_coordinates;
}

void Triangle::edge_integration_points (const Mesh_Prototype & mesh, int n1, int n2, double ** points, double * weigths, double * jac)
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

void Triangle::segment_integration_points (const Mesh_Prototype & mesh, double * coordinates_1, double * coordinates_2, double ** points, double * weigths, double * jac)
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

	double r = 0.0;

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

int Triangle::edge_amount_of_integration_points ()
{
	return 4;
}

int Triangle::get_isoline_points (const Mesh_Prototype & mesh, double value, double * q, double * c1, double * c2)
{
	int counter = 0;
	double cn0[2];
	double cn1[2];

	// check three sides
	for (int i = 0; i < 3 && counter < 2; i++)
	{
		int k0 = i;
		int k1 = (i + 1) % 3;
		// get nodes' points
		mesh.get_node_coordinates (defining_nodes[k0], cn0);
		mesh.get_node_coordinates (defining_nodes[k1], cn1);
		// get function's equation
		double A = q[k0] * Alpha->Elem (k0, 1) + q[k1] * Alpha->Elem (k1, 1);
		double B = q[k0] * Alpha->Elem (k0, 2) + q[k1] * Alpha->Elem (k1, 2);
		double C = q[k0] * Alpha->Elem (k0, 0) + q[k1] * Alpha->Elem (k1, 0);
		C -= value;
		
		double x = cn0[0];
		double y = cn0[1];
		{
			// if x = const
			if (fabs (cn0[0] - cn1[0]) < 1e-9)
			{
				y = -(A * x + C) / B;
			}
			else
			{
				// if y = const
				if (fabs (cn0[1] - cn1[1]) < 1e-9)
				{
					x = -(B * y + C) / A;
				}
				else // if y = kx + b
				{
					// get side's equation
					double k = (cn1[1] - cn0[1]) / (cn1[0] - cn0[0]);
					double b = (cn0[1] * cn1[0] - cn0[0] * cn1[1]) / (cn1[0] - cn0[0]);

					// get point of intersection
					x = -(B * b + C) / (A + B * k);
					y = k * x + b;
				}
			}
		}

		// check that the point is on the side of the triangle	
		bool intersection = false;
		double sp = (x - cn1[0]) * (x - cn0[0]) + (y - cn1[1]) * (y - cn0[1]);
		if (sp < 1e-7)
			intersection = true;
		// if it is, increase counter and save the point
		if (intersection)
		{
			if (counter == 0)
			{
				c1[0] = x;
				c1[1] = y;
				counter++;
			}
			if (counter == 1)
			{
				// if it is not the same point
				if (fabs (x - c1[0]) > 1e-7 || fabs (y - c1[1]) > 1e-7)
				{
					c2[0] = x;
					c2[1] = y;
					counter++;
				}
			}
		}
	}
	return counter;
}

int Triangle::get_function_by_node (int n)
{
	for (int i = 0; i < n_def_nodes; i++)
	{
		if (defining_nodes[i] == n)
		{
			return global_functions[i];
		}
	}
	return -1;
}

int Triangle::get_function_by_edge (int n1, int n2)
{
	return -1;
}

