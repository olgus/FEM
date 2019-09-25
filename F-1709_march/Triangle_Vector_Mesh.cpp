#include "Triangle_Vector_Mesh.h"

Triangle_Vector_Mesh::Triangle_Vector_Mesh ()
{
	dim = 2;

	coord0 = new double[dim];
	coordN = new double[dim];
	n_axis = new int[dim];
}

Triangle_Vector_Mesh::Triangle_Vector_Mesh (const Triangle_Vector_Mesh & mesh)
{
	// if memory has been allocated, free it
	if (coord0 != NULL)
		delete[] coord0;
	if (coordN != NULL)
		delete[] coordN;
	if (n_axis != NULL)
		delete[] n_axis;

	// reset dimentionality
	dim = mesh.dim;

	// allocate the memory
	// and copy triangular_Mesh's data
	coord0 = new double[dim];
	for (int i = 0; i < dim; i++)
		coord0[i] = mesh.coord0[i];

	coordN = new double[dim];
	for (int i = 0; i < dim; i++)
		coordN[i] = mesh.coordN[i];

	n_axis = new int[dim];
	for (int i = 0; i < dim; i++)
		n_axis[i] = mesh.n_axis[i];

	// clear nodes and copy triangular_Mesh's nodes
	nodes.clear ();
	n_nodes = mesh.n_nodes;
	nodes.insert (nodes.begin (), mesh.nodes.begin (), mesh.nodes.end ());

	// clear elements and copy triangular_Mesh's elements
	elements.clear ();
	n_elements = mesh.n_elements;
	for (const auto& e : mesh.elements)
		elements.push_back (std::make_unique<Triangle_Vector> (*e));
}

Triangle_Vector_Mesh::~Triangle_Vector_Mesh ()
{
}

void Triangle_Vector_Mesh::get_basis_function_value (int k_element, int k_func, double * coordinates, double * f_value)
{
	elements[k_element]->get_basis_function_value (k_func, coordinates, f_value);
}

int Triangle_Vector_Mesh::numerate_functions ()
{
	int functions[6];

	int n_edges = elements[0]->get_amount_edges ();
	int * local_edges = new int[n_edges];
	
	for (int i = 0; i < n_elements; i++)
	{
		// get numbers of edges of the element
		get_element_edges (i, local_edges);
		// number of the function is number of the edge multiplied by two and another + 1
		// since there are two functions associated with each edge
		// edges are numerated by minimum node's number
		// nodes within triangles are also sorted by number in increasing order
		// so it should be fine, maybe
		for (int j = 0; j < n_edges; j++)
		{
			functions[j * 2] = local_edges[j] * 2;
			functions[j * 2 + 1] = local_edges[j] * 2 + 1;
		}

		elements[i]->set_global_functions (functions);
	}
	return 2 * edges->get_n_entries ();

	delete[] local_edges;
}
