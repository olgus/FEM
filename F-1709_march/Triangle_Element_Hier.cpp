#include "Triangle_Element_Hier.h"

Triangle_Hier::Triangle_Hier ()
{
	n_functions = 6;
	order = 2;
}

Triangle_Hier::Triangle_Hier (const Triangle_Hier & triangle)
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
	area = triangle.area;

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
}

Triangle_Hier::~Triangle_Hier ()
{
}

double Triangle_Hier::get_basis_function_value (int k_func, double * coordinates)
{
	// first 3 basis functions are linear node ones
	// L[i](x, y) = al0[i] + al1[i] * x + al2[i] * y

	MathVector * L = new MathVector (n_functions);
	for (int i = 0; i < 3; i++)
	{
		L->setElem (i, Alpha->Elem (i, 0) + Alpha->Elem (i, 1) * coordinates[0] + Alpha->Elem (i, 2) * coordinates[1]);
	}

	// next 3 basis functions are square edge ones
	L->setElem (3, L->getElem (0) * L->getElem (1));
	L->setElem (4, L->getElem (0) * L->getElem (2));
	L->setElem (5, L->getElem (1) * L->getElem (2));

	double val = L->getElem (k_func);
	delete L;
	
	return val;
}

double Triangle_Hier::get_basis_function_derivative (int k_func, int k_var, double * coordinates)
{
	MathVector * L = new MathVector (n_functions);
	for (int i = 0; i < 3; i++)
	{
		L->setElem (i, Alpha->Elem (i, 0) + Alpha->Elem (i, 1) * coordinates[0] + Alpha->Elem (i, 2) * coordinates[1]);
	}

	switch (k_var)
	{
	case 0:
		L->setElem (3, Alpha->Elem (0, 1) * L->getElem (1) + L->getElem (0) * Alpha->Elem (1, 1));
		L->setElem (4, Alpha->Elem (0, 1) * L->getElem (2) + L->getElem (0) * Alpha->Elem (2, 1));
		L->setElem (5, Alpha->Elem (1, 1) * L->getElem (2) + L->getElem (1) * Alpha->Elem (2, 1));
		break;
	case 1:
		L->setElem (3, Alpha->Elem (0, 2) * L->getElem (1) + L->getElem (0) * Alpha->Elem (1, 2));
		L->setElem (4, Alpha->Elem (0, 2) * L->getElem (2) + L->getElem (0) * Alpha->Elem (2, 2));
		L->setElem (5, Alpha->Elem (1, 2) * L->getElem (2) + L->getElem (1) * Alpha->Elem (2, 2));
		break;
	}

	for (int i = 0; i < 3; i++)
	{
		switch (k_var)
		{
		case 0:
			L->setElem (i, Alpha->Elem (i, 1));
			break;
		case 1:
			L->setElem (i, Alpha->Elem (i, 2));
			break;
		}
	}

	double val = L->getElem (k_func);
	delete L;

	return val;
}

int Triangle_Hier::get_function_nodes (int k_func, int * nodes)
{
	int n = 0;
	// 1 for first 3 functions
	if (k_func < 3)
	{
		n = 1;
		nodes[0] = k_func;
	}
	// 2 for the second three functions
	if (k_func >= 3)
	{
		n = 2;
		nodes[0] = (k_func + 1) / 6;
		nodes[1] = (k_func + 2) / 3;
	}
	return n;
}

int Triangle_Hier::get_amount_second_condition ()
{
	return 3;
}

void Triangle_Hier::get_edge_functions (const Mesh_Prototype & mesh, int n1, int n2, int * functions)
{
	// three functions for hierarchical triangle element
	int node;
	int counter = 0;
	int n1_local, n2_local;
	n1_local = n2_local = -1;
	for (int i = 0; i < n_def_nodes; i++)
	{
		node = defining_nodes[i];
		if (node == n1 || node == n2)
		{
			functions[counter] = i;
			counter++;
		}
		if (node == n1 )
		{
			n1_local = i;
		}
		if (node == n2)
		{
			n2_local = i;
		}
	}

	if (n1_local == 0)
	{
		if (n2_local == 1)
			functions[counter] = 3;
		else
			functions[counter] = 4;
	}
	else
	{
		functions[counter] = 5;
	}
}

int Triangle_Hier::function_boundary (const Mesh_Prototype & mesh, int k_func)
{
	// by default check if def_nodes[k_func] is on the boundary

	int dim = mesh.get_dimentionality ();
	double * coord_node = new double[dim];
	double * coord_0 = new double[dim];
	double * coord_N = new double[dim];
	int boundary = -1;

	if (k_func < 3)
	{
		int node = defining_nodes[k_func];
		mesh.get_node_coordinates (node, coord_node);
		mesh.get_0_boundaries (coord_0);
		mesh.get_N_boundaries (coord_N);
		for (int d = 0; d < dim; d++)
		{
			if (abs (coord_node[d] - coord_0[d]) < ZERO_boundary)
			{
				boundary = d * dim;
			}

			if (abs (coord_node[d] - coord_N[d]) < ZERO_boundary)
			{
				boundary = d * dim + 1;
			}
		}
	}
	else
	{
		int node1, node2;
		if (k_func == 5)
		{
			node1 = defining_nodes[1];
		}
		else
		{
			node1 = defining_nodes[0];
		}
		if (k_func == 3)
		{
			node2 = defining_nodes[1];
		}
		else
		{
			node2 = defining_nodes[2];
		}
		double * coord_node2 = new double[dim];
		double * coord_node_edge = new double[dim];

		mesh.get_node_coordinates (node1, coord_node);
		mesh.get_node_coordinates (node2, coord_node2);
		mesh.get_0_boundaries (coord_0);
		mesh.get_N_boundaries (coord_N);

		for (int i = 0; i < dim; i++)
		{
			coord_node_edge[i] = coord_node2[i] - 0.5 * (coord_node2[i] - coord_node[i]);
		}

		for (int d = 0; d < dim; d++)
		{
			if (abs (coord_node_edge[d] - coord_0[d]) < ZERO_boundary)
			{
				boundary = d * dim;
			}

			if (abs (coord_node_edge[d] - coord_N[d]) < ZERO_boundary)
			{
				boundary = d * dim + 1;
			}
		}
		delete[] coord_node2;
		delete[] coord_node_edge;
	}

	delete[] coord_node;
	delete[] coord_0;
	delete[] coord_N;
	return boundary;
}

int Triangle_Hier::get_isoline_points (const Mesh_Prototype & mesh, double value, double * q, double * c1, double * c2)
{
	int counter = 0;
	// nodes' coordinates
	double cn0[2];
	double cn1[2];
	// intersections
	double ic[2][2];

	// check three sides
	for (int i = 0; i < 3 && counter < 2; i++)
	{
		int k0 = i / 2;
		int k1 = (i + 3) / 2;
		int k2 = i + 3;
		// get nodes' points
		mesh.get_node_coordinates (defining_nodes[k0], cn0);
		mesh.get_node_coordinates (defining_nodes[k1], cn1);
		// get function's equation
		double a[] = { (Alpha->Elem (k0, 1) * Alpha->Elem (k1, 1)) * q[k2],
			(Alpha->Elem (k0, 2) * Alpha->Elem (k1, 2)) * q[k2],
			(Alpha->Elem (k0, 1) * Alpha->Elem (k1, 2) + Alpha->Elem (k0, 2) * Alpha->Elem (k1, 1)) * q[k2],
			Alpha->Elem (k0, 1) * q[k0] + Alpha->Elem (k1, 1) * q[k1] + (Alpha->Elem (k0, 0) * Alpha->Elem (k1, 1) + Alpha->Elem (k1, 0) * Alpha->Elem (k0, 1)) * q[k2],
			Alpha->Elem (k0, 2) * q[k0] + Alpha->Elem (k1, 2) * q[k1] + (Alpha->Elem (k0, 0) * Alpha->Elem (k1, 2) + Alpha->Elem (k1, 0) * Alpha->Elem (k0, 2)) * q[k2],
			Alpha->Elem (k0, 0) * q[k0] + Alpha->Elem (k1, 0) * q[k1] + Alpha->Elem (k0, 0) * Alpha->Elem (k1, 0) * q[k2] - value};

		double A, B, C;
		double x = cn0[0];
		double y = cn0[1];
		{
			// if x = const
			if (fabs (cn0[0] - cn1[0]) < 1e-9)
			{
				A = a[1];
				B = a[2] * x + a[4];
				C = a[0] * x * x + a[3] * x + a[5];

				ic[0][0] = ic[1][0] = x;
				if (fabs (A) > 1e-7)
				{
					ic[0][1] = (-B + sqrt (B * B - 4.0 * A * C)) / (2.0 * A);
					ic[1][1] = (-B - sqrt (B * B - 4.0 * A * C)) / (2.0 * A);
				}
				else
				{
					ic[0][1] = ic[1][1] = -C / B;
				}
			}
			else
			{
				// if y = const
				if (fabs (cn0[1] - cn1[1]) < 1e-9)
				{
					A = a[0];
					B = a[2] * y + a[3];
					C = a[1] * y * y + a[4] * y + a[5];
					ic[0][1] = ic[1][1] = y;
					if (fabs (A) > 1e-7)
					{
						ic[0][0] = (-B + sqrt (B * B - 4.0 * A * C)) / (2.0 * A);
						ic[1][0] = (-B - sqrt (B * B - 4.0 * A * C)) / (2.0 * A);
					}
					else
					{
						ic[0][0] = ic[1][0] = -C / B;
					}
				}
				else // if y = kx + b
				{
					// get side's equation
					double k = (cn1[1] - cn0[1]) / (cn1[0] - cn0[0]);
					double b = (cn0[1] * cn1[0] - cn0[0] * cn1[1]) / (cn1[0] - cn0[0]);

					A = a[0] + a[1] * k * k + a[2] * k;
					B = 2.0 * k * b * a[1] + a[4] * k + a[2] * b + a[3];
					C = a[5] + a[1] * b * b + a[4] * b;
					// get points of intersection
					if (fabs (A) > 1e-7)
					{
						ic[0][0] = (-B + sqrt (B * B - 4.0 * A * C)) / (2.0 * A);
						ic[1][0] = (-B - sqrt (B * B - 4.0 * A * C)) / (2.0 * A);
					}
					else
					{
						ic[0][0] = ic[1][0] = -C / B;
					}
					ic[0][1] = k * ic[0][0] + b;
					ic[1][1] = k * ic[1][0] + b;
				}
			}
		}
		// we got two points, find only one

		// check that the point is on the side of the triangle	
		int k_intersection = -1;
		if (((ic[0][0] - cn1[0]) * (ic[0][0] - cn0[0]) + (ic[0][1] - cn1[1]) * (ic[0][1] - cn0[1])) < 1e-7)
		{
			k_intersection = 0;
		}
		if (((ic[1][0] - cn1[0]) * (ic[1][0] - cn0[0]) + (ic[1][1] - cn1[1]) * (ic[1][1] - cn0[1])) < 1e-7)
		{
			k_intersection = 1;
		}

		// if it is, increase counter and save the point
		if (k_intersection != -1)
		{
			if (counter == 0)
			{
				c1[0] = ic[k_intersection][0];
				c1[1] = ic[k_intersection][1];
				counter++;
			}
			if (counter == 1)
			{
				// if it is not the same point
				if (fabs (ic[k_intersection][0] - c1[0]) > 1e-7 || fabs (ic[k_intersection][1] - c1[1]) > 1e-7)
				{
					c2[0] = ic[k_intersection][0];
					c2[1] = ic[k_intersection][1];
					counter++;
				}
			}
		}
	}
	return counter;
}

int Triangle_Hier::get_function_by_node (int n)
{
	// three functions for hierarchical triangle element
	for (int i = 0; i < n_def_nodes; i++)
	{
		if (defining_nodes[i] == n)
		{
			return global_functions[i];
		}
	}
	return -1;
}

int Triangle_Hier::get_function_by_edge (int n1, int n2)
{
	int node;
	int n1_local, n2_local;
	n1_local = n2_local = -1;
	for (int i = 0; i < n_def_nodes; i++)
	{
		node = defining_nodes[i];
		if (node == n1)
		{
			n1_local = i;
		}
		if (node == n2)
		{
			n2_local = i;
		}
	}

	if (n1_local == 0)
	{
		if (n2_local == 1)
			return global_functions[3];
		else			
			return global_functions[4];
	}
	else
	{
		return global_functions[5];
	}
}
