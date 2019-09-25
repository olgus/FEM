#include "Prototype_Element.h"

Element::Element ()
{
	n_def_nodes = 0;
	n_base_nodes = 0;
	defining_nodes = NULL;
	base_nodes = NULL;
	n_functions = 0;
	global_functions = NULL;
	dim = 0;
	order = 1;

	G = NULL;
	M = NULL;
}

Element::Element (const Element & element)
{
	if (n_def_nodes != element.n_def_nodes) // if triangle's amount of nodes doesn't equal to amount of nodes of current triangle
	{
		if (defining_nodes != NULL && n_def_nodes > 0) // if memory was even allocated
			delete[] defining_nodes; // free it
		
		n_def_nodes = element.n_def_nodes; // set new amount of nodes
		defining_nodes = new int[n_def_nodes]; // allocate the memory
		
		if (base_nodes != NULL && n_base_nodes > 0) // if memory was even allocated
			delete[] base_nodes;

		n_base_nodes = element.n_base_nodes; // set new amount of nodes
		base_nodes = new int[n_base_nodes]; // allocate the memory
	}

	if (global_functions != NULL && n_functions > 0) // if memory was even allocated
		delete[] global_functions; // free it

	if (n_functions > 0 && element.global_functions != NULL)
	{
		global_functions = new int[n_functions]; // allocate the memory
		for (int i = 0; i < n_functions; i++) // copy functions
			global_functions[i] = element.global_functions[i];
	}

	for (int i = 0; i < n_def_nodes; i++) // copy nodes
		defining_nodes[i] = element.defining_nodes[i];

	for (int i = 0; i < n_base_nodes; i++) // copy nodes
		base_nodes[i] = element.base_nodes[i];

	area = element.area; // copy area

	// set G if it exists
	if (element.G != NULL)
		G->Full_Copy (*element.G);
	else
		G = NULL;

	// set M if it exists
	if (element.M != NULL)
		M->Full_Copy (*element.M);
	else
		M = NULL;
}

Element::~Element ()
{
	if (defining_nodes != NULL && n_def_nodes > 0)
		delete[] defining_nodes;

	if (base_nodes != NULL && n_base_nodes > 0)
		delete[] base_nodes;

	if (global_functions != NULL && n_functions > 0)
		delete[] global_functions;

	if (G != NULL)
		delete G;
	if (M != NULL)
		delete M;
}

void Element::set_element (int n, int * nodes)
{
	if (defining_nodes != NULL) // if element has been set before
		delete[] defining_nodes; // free the memory

	n_def_nodes = n; // set amount of nodes
	defining_nodes = new int[n_def_nodes]; // allocate memory

	for (int i = 0; i < n_def_nodes; i++) // fill nodes
		defining_nodes[i] = nodes[i];
}

void Element::set_area (int a)
{
	area = a;
}

void Element::set_global_functions (int * funcs)
{
	if (global_functions != NULL)
		delete[] global_functions; // free the memory

	global_functions = new int[n_functions];

	for (int i = 0; i < n_functions; i++)
	{
		global_functions[i] = funcs[i];
	}
}

void Element::set_global_functions ()
{
	// global functions by default are numbers of nodes
	if (global_functions != NULL) // if memory hasn't been allocated 
		delete[] global_functions;
	global_functions = new int[n_functions]; // do it
	for (int i = 0; i < n_functions; i++) // fill nodes
		global_functions[i] = defining_nodes[i];
}

void Element::set_global_function (int k_func_local, int k_func_global)
{
	if (k_func_local < n_functions)
	{
		global_functions[k_func_local] = k_func_global;
	}
}

void Element::get_mass_center (const Mesh_Prototype & mesh, double * coordinates)
{
	for (int i = 0; i < dim; i++)
	{
		coordinates[i] = 0.0;
	}

	double * c_base_node;
	c_base_node = new double[dim];
	for (int k = 0; k < n_base_nodes; k++)
	{
		mesh.get_node_coordinates (base_nodes[k], c_base_node);
		for (int i = 0; i < dim; i++)
		{
			coordinates[i] += c_base_node[i];
		}
	}

	for (int i = 0; i < dim; i++)
	{
		coordinates[i] /= n_base_nodes;
	}
	delete[] c_base_node;
}

int Element::get_amount_of_def_nodes () const
{
	return n_def_nodes;
}

int Element::get_amount_of_base_nodes () const
{
	return n_base_nodes;
}

void Element::get_def_nodes (int * n)
{
	for (int i = 0; i < n_def_nodes; i++)
		n[i] = defining_nodes[i];
}

int Element::get_def_node (int k_node)
{
	return defining_nodes[k_node];
}

void Element::get_base_nodes (int * n)
{
	for (int i = 0; i < n_base_nodes; i++)
		n[i] = base_nodes[i];
}

int Element::get_area () const
{
	return area;
}

int Element::get_node (int k_node) const
{
	return defining_nodes[k_node];
}

bool Element::node_belongs (int k_node)
{
	bool result = false;
	for (int i = 0; i < n_def_nodes && !result; i++)
	{
		if (defining_nodes[i] == k_node)
			result = true;
	}
	return result;
}

void Element::sort_def_nodes ()
{
	std::sort (defining_nodes, defining_nodes + n_def_nodes);
	std::sort (base_nodes, base_nodes + n_base_nodes);
}

int Element::get_function_global_number (int k_func)
{
	if (k_func < n_functions)
		return global_functions[k_func];
	else
		return -1;
}

int Element::get_order ()
{
	return order;
}

int Element::get_amount_second_condition ()
{
	return 0;
}

double Element::get_geometrical_area ()
{
	return 0.0;
}

void Element::get_G_local_matrix (Matrix * G_matrix)
{
	G_matrix->Copy (*G);
}

void Element::get_M_local_matrix (Matrix * M_matrix)
{
	M_matrix->Copy (*M);
}

void Element::get_edge_functions (const Mesh_Prototype & mesh, int n1, int n2, int * functions)
{
}

bool Element::point_inside (const Mesh_Prototype & mesh, double * coordinates)
{
	return false;
}

void Element::get_local_function_values (double * coordinates, MathVector * v)
{
	for (int i = 0; i < n_functions; i++)
	{
		v->setElem (i, get_basis_function_value (i, coordinates));
	}
}

void Element::get_local_function_values (double * coordinates, Matrix * v)
{
	double * value = new double[dim];
	for (int i = 0; i < n_functions; i++)
	{
		get_basis_function_value (i, coordinates, value);
		for (int j = 0; j < dim; j++)
		{
			v->setElem (i, j, value[j]);
		}
	}
	delete[] value;
}

void Element::get_local_function_first_derivative_values (double * coordinates, int k_var, MathVector * v)
{
	for (int i = 0; i < n_functions; i++)
	{
		v->setElem (i, get_basis_function_derivative (i, k_var, coordinates));
	}
}

void Element::prepare (const Mesh_Prototype & mesh)
{
	inside_prepare (mesh);
	prepare_GM ();
}

void Element::inside_prepare (const Mesh_Prototype & mesh)
{
}

void Element::prepare_GM ()
{
	if (G == NULL)
		G = new Matrix (n_functions, n_functions);
	if (M == NULL)
		M = new Matrix (n_functions, n_functions);

	int v[1];
	for (int j = 0; j < n_functions; j++)
	{
		for (int i = 0; i < n_functions; i++)
		{
			v[0] = 0;
			G->setElem (i, j, integrate (i, j, v));
			v[0] = 1;
			M->setElem (i, j, integrate (i, j, v));
		}
	}
}

double Element::integrate (int i, int j, int * m)
{
	return 0.0;
}

double Element::get_function_value (int i, int j, int * m, double * coordinates)
{
	double r = 0.0;

	switch (m[0])
	{
	case 0: // g(i,j)
		for (int k = 0; k < dim; k++)
		{
			r += get_basis_function_derivative (i, k, coordinates) * get_basis_function_derivative (j, k, coordinates);
		}
		break;
	case 1: // m(i,j)
		r += get_basis_function_value (i, coordinates) * get_basis_function_value (j, coordinates);
		break;
	case 3: // bf i * bf j * bf k
		r += get_basis_function_value (i, coordinates) * get_basis_function_value (j, coordinates) * get_basis_function_value (m[1], coordinates);
		break;
	case 9:
	{
		int k, v1, v2;

		k = m[1];
		v1 = m[2];
		v2 = m[3];

		r += get_basis_function_derivative (k, v1, coordinates) * get_basis_function_derivative (j, v2, coordinates) * get_basis_function_value (i, coordinates);
		break;
	}
	case 10: // D(i,j) / dm
		r += get_basis_function_value (i, coordinates) * get_basis_function_derivative (j, 0, coordinates);
		break;
	case 20: // D(i,j) / dm
		r += get_basis_function_value (i, coordinates) * get_basis_function_derivative (j, 1, coordinates);
		break;
	default:
		r = 0.0;
	}
	return r;
}

double Element::get_basis_function_value (int k_func, double * coordinates)
{
	return 0.0;
}

void Element::get_basis_function_value (int k_func, double * coordinates, double * f_value)
{
}

double Element::get_basis_function_derivative (int k_func, int k_var, double * coordinates)
{
	return 0.0;
}

int Element::function_boundary (const Mesh_Prototype & mesh, int k_func)
{
	// by default check if def_nodes[k_func] is on the boundary
	int node = defining_nodes[k_func];
	int dim = mesh.get_dimentionality ();
	double * coord_node = new double[dim];
	double * coord_0 = new double[dim];
	double * coord_N = new double [dim];
	mesh.get_node_coordinates (node, coord_node);
	mesh.get_0_boundaries (coord_0);
	mesh.get_N_boundaries (coord_N);
	int boundary = -1;
	for (int d = 0; d < dim; d++)
	{
		if (abs (coord_node[d] - coord_0[d]) < ZERO_boundary)
		{
			boundary = d * 2;
		} 

		if (abs (coord_node[d] - coord_N[d]) < ZERO_boundary)
		{
			boundary = d * 2 + 1;
		}
	}

	delete[] coord_node;
	delete[] coord_0;
	delete[] coord_N;
	return boundary;
}

int Element::get_amount_non_zero_functions ()
{
	return n_functions;
}

bool Element::edge_exists (int n1, int n2)
{
	return false;
}

int Element::get_function_nodes (int k_func, int * nodes)
{
	// by default return base_nodes[k_func]
	nodes[0] = k_func;
	return 1;
}

void Element::get_edges (int ** edges)
{
}

int Element::get_amount_edges ()
{
	return 0;
}

void Element::get_D (int var_der, Matrix * D) // returns integrate (bfi * d bfj/ d var_der)
{
	double r;
	int m[] = { 10 * (var_der + 1) };
	for (int i = 0, i_end = D->Size0 (); i < i_end; i++)
	{
		for (int j = 0, j_end = D->Size1 (); j < j_end; j++)
		{
			r = integrate (i, j, m);
			D->setElem (i, j, r);
		}
	}
}

void Element::get_DDD (MathVector * D)
{
	int m[] = {3, 0};
	for (int i = 0; i < n_functions; i++)
	{
		for (int j = 0; j < n_functions; j++)
		{
			for (int k = 0; k < n_functions; k++)
			{
				m[1] = k;
				D->setElem (i * n_functions * n_functions + j * n_functions + k, integrate (i, j, m));
			}
		}
	}
}

void Element::get_DDF (int var_der1, int var_der2, MathVector * D)
{
	int m[] = { 9, 0, var_der1, var_der2 };
	for (int i = 0; i < n_functions; i++)
	{
		for (int j = 0; j < n_functions; j++)
		{
			for (int k = 0; k < n_functions; k++)
			{
				m[1] = k;
				D->setElem (i * n_functions * n_functions + j * n_functions + k, integrate (i, j, m));
			}
		}
	}
}

void Element::get_DF (int var_der, Matrix * D)
{
	get_D (var_der, D); // get D
	D->Transpose (); // transpose it
}

void Element::get_DDF (int k, int var_der1, int var_der2, Matrix * D)
{
	double v;
	int m[] = { 9, k, var_der1, var_der2 };

	for (int i = 0; i < n_functions; i++)
	{
		for (int j = 0; j < n_functions; j++)
		{
			v = integrate (i, j, m);
			D->setElem (i, j, v);
		}
	} 
}

int Element::amount_of_integration_points ()
{
	return 0;
}

void Element::integration_points (double ** points, double * weigths, double * jac)
{
} 

int Element::edge_amount_of_integration_points ()
{
	return 0;
}

void Element::edge_integration_points (const Mesh_Prototype & mesh, int n1, int n2, double ** points, double * weigths, double * jac)
{
	// nothing for 1 and 3 - dim tasks
}

void Element::segment_integration_points (const Mesh_Prototype & mesh, double * coordinates_1, double * coordinates_2, double ** points, double * weigths, double * jac)
{
}

int Element::get_isoline_points (const Mesh_Prototype & mesh, double value, double * q, double * c1, double * c2)
{
	return 0;
}

int Element::get_function_by_node (int n)
{
	return 0;
}

int Element::get_function_by_edge (int n1, int n2)
{
	return 0;
}

void Element::recalc_matrices ()
{
	prepare_GM (); // TEST if it works correctly for elements with overriden prepare_GM
}

void Element::set_amount_def_nodes (int n)
{
	n_def_nodes = n;
}

void Element::set_element (int * nodes) // set nodes for the elements
{
	if (n_def_nodes > 0)
	{
		if (defining_nodes == NULL) // if memory hasn't been allocated 
			defining_nodes = new int[n_def_nodes]; // do it
		for (int i = 0; i < n_def_nodes; i++) // fill nodes
			defining_nodes[i] = nodes[i];
	}
}

void Element::set_base_nodes (int * nodes)
{
	if (n_base_nodes > 0)
	{
		if (base_nodes == NULL) // if memory hasn't been allocated 
			base_nodes = new int[n_base_nodes]; // do it
		for (int i = 0; i < n_base_nodes; i++) // fill nodes
			base_nodes[i] = nodes[i];
	}
}

void Element::set_nodes (int * nodes)
{
	if (n_def_nodes > 0)
	{
		if (defining_nodes == NULL) // if memory hasn't been allocated 
			defining_nodes = new int[n_def_nodes]; // do it
		for (int i = 0; i < n_def_nodes; i++) // fill nodes
			defining_nodes[i] = nodes[i];
	}
}
