#include "Triangle_Vector.h"

double Triangle_Vector::get_scalar_product (int k1, int k2, double * coordinates)
{
	double f1[2], f2[2];
	get_basis_function_value (k1, coordinates, f1);
	get_basis_function_value (k2, coordinates, f2);

	double r = 0.0;
	for (int i = 0; i < 2; i++)
	{
		r += f1[i] * f2[i];
	}
	return r;
}

Triangle_Vector::Triangle_Vector ()
{
	n_functions = 6;
	order = 1;
}

Triangle_Vector::Triangle_Vector (const Triangle_Vector & triangle)
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

Triangle_Vector::~Triangle_Vector ()
{
}

double Triangle_Vector::get_basis_function_value (int k_func, double * coordinates)
{
	printf ("ERROR: nothing-function in Triangle_Vector\n");
	return 0.0;
}

void Triangle_Vector::get_basis_function_value (int k_func, double * coordinates, double * f_value)
{
	// L[i](x, y) = al0[i] + al1[i] * x + al2[i] * y
	MathVector * L = new MathVector (3);
	for (int i = 0; i < 3; i++)
	{
		L->setElem (i, Alpha->Elem (i, 0) + Alpha->Elem (i, 1) * coordinates[0] + Alpha->Elem (i, 2) * coordinates[1]);
	}

	int c[2];
	switch (k_func)
	{
	case 0:
		c[0] = 1;
		c[1] = 0;
		break;
	case 1:
		c[0] = 0;
		c[1] = 1;
		break;
	case 2:
		c[0] = 2;
		c[1] = 0;
		break;
	case 3:
		c[0] = 0;
		c[1] = 2;
		break;
	case 4:
		c[0] = 2;
		c[1] = 1;
		break;
	case 5:
		c[0] = 1;
		c[1] = 2;
		break;
	}
	int e = k_func % 3;

	for (int i = 0; i < 2; i++)
	{
		f_value[i] = L->getElem (c[0]) * Alpha->Elem (c[1], i + 1);
	}
}

double Triangle_Vector::get_basis_function_derivative (int k_func, int k_var, double * coordinates)
{
	// return rot
	// it's only z component

	int c[2];
	switch (k_func)
	{
	case 0:
		c[0] = 1;
		c[1] = 0;
		break;
	case 1:
		c[0] = 0;
		c[1] = 1;
		break;
	case 2:
		c[0] = 2;
		c[1] = 0;
		break;
	case 3:
		c[0] = 0;
		c[1] = 2;
		break;
	case 4:
		c[0] = 2;
		c[1] = 1;
		break;
	case 5:
		c[0] = 1;
		c[1] = 2;
		break;
	}
	double r = Alpha->Elem (c[0], 1) * Alpha->Elem (c[1], 2);
	r -= Alpha->Elem (c[0], 2) * Alpha->Elem (c[1], 1);
		
	return r;
}

double Triangle_Vector::get_function_value (int i, int j, int * m, double * coordinates)
{
	double r = 0.0;

	switch (m[0])
	{
	case 0: // g(i,j)
		r += get_basis_function_derivative (i, 0, coordinates) * get_basis_function_derivative (j, 0, coordinates);
		break;
	case 1: // m(i,j)
		r += get_scalar_product (i, j, coordinates);
		break;
	default:
		r = 0.0;
	}
	return r;
}

void Triangle_Vector::get_edge_functions (const Mesh_Prototype & mesh, int n1, int n2, int * functions)
{
	// two functions for vector triangle element
	functions[0] = functions[1] = -1;

	int edge = -1;
	for (int i = 0; i < 3; i++)
	{
		if (n1 == defining_nodes[i] || n2 == defining_nodes[i])
			edge += i;
	}

	if (edge != -1)
	{
		functions[0] = edge * 2;
		functions[1] = edge * 2 + 1;
	}
}
