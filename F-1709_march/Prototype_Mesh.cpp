#include "Prototype_Mesh.h"

template class Mesh <Node<Point_Prototype>, Element>;
template class Mesh <Node<Point_2D>, Element>;
template class Mesh <Node<Point_2D_Polar>, Element>;
template class Mesh <Node<Point_3D>, Element>;
template class Mesh <Node_2D, Element>;
template class Mesh <Node_2D_Polar, Element>;
template class Mesh <Node_1D, Element>;
template class Mesh <Node_3D, Element>;

template class Mesh <Node<Point_Prototype>, Triangle>;
template class Mesh <Node<Point_2D>, Triangle>;
template class Mesh <Node<Point_2D_Polar>, Triangle>;
template class Mesh <Node<Point_3D>, Triangle>;
template class Mesh <Node_2D, Triangle>;
template class Mesh <Node_2D_Polar, Triangle>;
template class Mesh <Node_1D, Triangle>;
template class Mesh <Node_3D, Triangle>;

template class Mesh <Node<Point_Prototype>, Triangle_Polar>;
template class Mesh <Node<Point_2D>, Triangle_Polar>;
template class Mesh <Node<Point_2D_Polar>, Triangle_Polar>;
template class Mesh <Node<Point_3D>, Triangle_Polar>;
template class Mesh <Node_2D, Triangle_Polar>;
template class Mesh <Node_2D_Polar, Triangle_Polar>;
template class Mesh <Node_1D, Triangle_Polar>;
template class Mesh <Node_3D, Triangle_Polar>;

template class Mesh <Node<Point_Prototype>, Triangle_Hier>;
template class Mesh <Node<Point_2D>, Triangle_Hier>;
template class Mesh <Node<Point_2D_Polar>, Triangle_Hier>;
template class Mesh <Node<Point_3D>, Triangle_Hier>;
template class Mesh <Node_2D, Triangle_Hier>;
template class Mesh <Node_2D_Polar, Triangle_Hier>;
template class Mesh <Node_1D, Triangle_Hier>;
template class Mesh <Node_3D, Triangle_Hier>;

template class Mesh <Node<Point_Prototype>, Triangle_Vector>;
template class Mesh <Node<Point_2D>, Triangle_Vector>;
template class Mesh <Node<Point_2D_Polar>, Triangle_Vector>;
template class Mesh <Node<Point_3D>, Triangle_Vector>;
template class Mesh <Node_2D, Triangle_Vector>;
template class Mesh <Node_2D_Polar, Triangle_Vector>;
template class Mesh <Node_1D, Triangle_Vector>;
template class Mesh <Node_3D, Triangle_Vector>;

template class Mesh <Node<Point_Prototype>, Rectangle_Element>;
template class Mesh <Node<Point_2D>, Rectangle_Element>;
template class Mesh <Node<Point_2D_Polar>, Rectangle_Element>;
template class Mesh <Node<Point_3D>, Rectangle_Element>;
template class Mesh <Node_2D, Rectangle_Element>;
template class Mesh <Node_2D_Polar, Rectangle_Element>;
template class Mesh <Node_1D, Rectangle_Element>;
template class Mesh <Node_3D, Rectangle_Element>;

template class Mesh <Node<Point_Prototype>, Rectangle_Element_3>;
template class Mesh <Node<Point_2D>, Rectangle_Element_3>;
template class Mesh <Node<Point_2D_Polar>, Rectangle_Element_3>;
template class Mesh <Node<Point_3D>, Rectangle_Element_3>;
template class Mesh <Node_2D, Rectangle_Element_3>;
template class Mesh <Node_2D_Polar, Rectangle_Element_3>;
template class Mesh <Node_1D, Rectangle_Element_3>;
template class Mesh <Node_3D, Rectangle_Element_3>;

template class Mesh <Node<Point_Prototype>, Rectangle_Element_Hermitian>;
template class Mesh <Node<Point_2D>, Rectangle_Element_Hermitian>;
template class Mesh <Node<Point_2D_Polar>, Rectangle_Element_Hermitian>;
template class Mesh <Node<Point_3D>, Rectangle_Element_Hermitian>;
template class Mesh <Node_2D, Rectangle_Element_Hermitian>;
template class Mesh <Node_2D_Polar, Rectangle_Element_Hermitian>;
template class Mesh <Node_1D, Rectangle_Element_Hermitian>;
template class Mesh <Node_3D, Rectangle_Element_Hermitian>;

template class Mesh <Node<Point_Prototype>, Element_1D_hier>;
template class Mesh <Node<Point_2D>, Element_1D_hier>;
template class Mesh <Node<Point_2D_Polar>, Element_1D_hier>;
template class Mesh <Node<Point_3D>, Element_1D_hier>;
template class Mesh <Node_2D, Element_1D_hier>;
template class Mesh <Node_2D_Polar, Element_1D_hier>;
template class Mesh <Node_1D, Element_1D_hier>;
template class Mesh <Node_3D, Element_1D_hier>;

template class Mesh <Node<Point_Prototype>, Element_1D_L1>;
template class Mesh <Node<Point_2D>, Element_1D_L1>;
template class Mesh <Node<Point_2D_Polar>, Element_1D_L1>;
template class Mesh <Node<Point_3D>, Element_1D_L1>;
template class Mesh <Node_2D, Element_1D_L1>;
template class Mesh <Node_2D_Polar, Element_1D_L1>;
template class Mesh <Node_1D, Element_1D_L1>;
template class Mesh <Node_3D, Element_1D_L1>;

template class Mesh <Node<Point_Prototype>, Element_1D_Hermitian>;
template class Mesh <Node<Point_2D>, Element_1D_Hermitian>;
template class Mesh <Node<Point_2D_Polar>, Element_1D_Hermitian>;
template class Mesh <Node<Point_3D>, Element_1D_Hermitian>;
template class Mesh <Node_2D, Element_1D_Hermitian>;
template class Mesh <Node_2D_Polar, Element_1D_Hermitian>;
template class Mesh <Node_1D, Element_1D_Hermitian>;
template class Mesh <Node_3D, Element_1D_Hermitian>;

template class Mesh <Node<Point_Prototype>, Prism>;
template class Mesh <Node<Point_2D>, Prism>;
template class Mesh <Node<Point_2D_Polar>, Prism>;
template class Mesh <Node<Point_3D>, Prism>;
template class Mesh <Node_2D, Prism>;
template class Mesh <Node_2D_Polar, Prism>;
template class Mesh <Node_1D, Prism>;
template class Mesh <Node_3D, Prism>;

template class Mesh <Node<Point_Prototype>, Cube>;
template class Mesh <Node<Point_2D>, Cube>;
template class Mesh <Node<Point_2D_Polar>, Cube>;
template class Mesh <Node<Point_3D>, Cube>;
template class Mesh <Node_2D, Cube>;
template class Mesh <Node_2D_Polar, Cube>;
template class Mesh <Node_1D, Cube>;
template class Mesh <Node_3D, Cube>;

template <class Node_Type, class Element_Type>
Mesh<Node_Type, Element_Type>::Mesh ()
{
	dim = 0;

	n_nodes = n_elements = 0;
	n_material_areas = 1;

	coord0 = NULL;
	coordN = NULL;

	n_axis = NULL;
	edges = NULL;
}

template <class Node_Type, class Element_Type>
Mesh<Node_Type, Element_Type>::Mesh (const Mesh & mesh)
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
	coord0 = new double[dim];
	for (int i = 0; i < dim; i++)
		coord0[i] = mesh.coord0[i]; // and copy triangular_Mesh's data

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
	//elements.insert (elements.begin (), mesh.elements.begin (), mesh.elements.end ());

	if (mesh.edges != NULL)
	{
		delete edges;
		edges = new Edge_Structure;
		edges->Copy (*(mesh.edges));
	}
}

template <class Node_Type, class Element_Type>
Mesh<Node_Type, Element_Type>::~Mesh ()
{
	if (coord0 != NULL)
		delete[] coord0;
	if (coordN != NULL)
		delete[] coordN;
	if (n_axis != NULL)
		delete[] n_axis;
	if (edges != NULL)
		delete edges;
	// delete nodes
	nodes.clear ();
	// delete elements
	elements.clear ();
}

template <class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::get_0_boundaries (double * boundaries) const
{
	for (int i = 0; i < dim; i++)
		boundaries[i] = coord0[i];
}

template <class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::get_N_boundaries (double * boundaries) const
{
	for (int i = 0; i < dim; i++)
		boundaries[i] = coordN[i];
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::get_area_0_boundaries (int area, double * c0) const
{
	int * e_nodes = NULL;
	int n_nodes;
	double * coordinates = new double [dim];

	for (int i = 0; i < dim; i++)
		c0[i] = 1e+20;

	for (int k_element = 0, k_element_end = get_n_elements (); k_element < k_element_end; k_element++) // go through elements
	{
		if (elements[k_element]->get_area () == area)
		{
			n_nodes = elements[k_element]->get_amount_of_base_nodes ();
			e_nodes = new int[n_nodes];
			elements[k_element]->get_base_nodes (e_nodes);
			for (int k = 0; k < n_nodes; k++)
			{
				nodes[e_nodes[k]].get_coordinates (coordinates);
				for (int i = 0; i < dim; i++)
				{
					if (coordinates[i] < c0[i])
						c0[i] = coordinates[i];
				}
			}
			delete[] e_nodes;
			e_nodes = NULL;
		}
	}

	if (e_nodes != NULL)
		delete[] e_nodes;
	delete[] coordinates;
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::get_area_N_boundaries (int area, double * cn) const
{
	int * e_nodes = NULL;
	int n_nodes;
	double * coordinates = new double[dim];

	for (int i = 0; i < dim; i++)
		cn[i] = -1e+20;

	for (int k_element = 0, k_element_end = get_n_elements (); k_element < k_element_end; k_element++) // go through elements
	{
		if (elements[k_element]->get_area () == area)
		{
			n_nodes = elements[k_element]->get_amount_of_base_nodes ();
			e_nodes = new int[n_nodes];
			elements[k_element]->get_base_nodes (e_nodes);
			for (int k = 0; k < n_nodes; k++)
			{
				nodes[e_nodes[k]].get_coordinates (coordinates);
				for (int i = 0; i < dim; i++)
				{
					if (coordinates[i] > cn[i])
						cn[i] = coordinates[i];
				}
			}
			delete[] e_nodes;
			e_nodes = NULL;
		}
	}

	if (e_nodes != NULL)
		delete[] e_nodes;
	delete[] coordinates;
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::get_amount_of_areas () const
{
	return n_material_areas;
}

template <class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::get_dimentionality () const
{
	return dim;
}

template<class Node_Type, class Element_Type>
bool Mesh<Node_Type, Element_Type>::edge_exists (int k_element, int n1, int n2)
{
	return elements[k_element]->edge_exists (n1, n2);
}

template<class Node_Type, class Element_Type>
bool Mesh<Node_Type, Element_Type>::get_isoline_section (int k_element, double * q, double value, double * c1, double * c2)
{
	return false;
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::move_node (int k_node, double * coordinates)
{
	nodes[k_node].set_coordinates (coordinates);
}

template<class Node_Type, class Element_Type>
double Mesh<Node_Type, Element_Type>::nodes_distance (int k1_node, int k2_node)
{
	return nodes[k1_node].distance (nodes[k2_node]);
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::refresh_mesh ()
{
	for (int i = 0, i_end = (int)elements.size (); i < i_end; i++)
		elements[i]->recalc_matrices ();
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::reset_n_nodes ()
{
	n_nodes = (int)nodes.size ();
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::reset_n_elements ()
{
	n_elements = (int)elements.size ();
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::copy_nodes (const Mesh_Prototype & mesh)
{
	double * c = new double[dim];
	for (int i = 0; i < mesh.get_n_nodes (); i++)
	{
		mesh.get_node_coordinates (i, c);
		Node_Type node;
		node.set_coordinates (c);
		nodes.push_back (node);
	}
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::reset_elements ()
{
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::get_amount_of_def_nodes (int k_element) const
{
	return elements[k_element]->get_amount_of_def_nodes ();
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::get_node_coordinates (int k_node, double * c) const
{
	nodes[k_node].get_coordinates (c);
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::prepare_element (int k_element)
{
	 elements[k_element]->prepare (*this);
}

template<class Node_Type, class Element_Type>
bool Mesh<Node_Type, Element_Type>::point_inside (int k_element, double * coordinates)
{
	return elements[k_element]->point_inside (*this, coordinates);
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::get_local_function_values (int k_element, double * coordinates, MathVector * v)
{
	 elements[k_element]->get_local_function_values (coordinates, v);
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::get_local_function_values (int k_element, double * coordinates, Matrix * v)
{
	elements[k_element]->get_local_function_values (coordinates, v);
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::get_local_function_first_derivative_values (int k_element, int k_var, double * coordinates, MathVector * v)
{
	 elements[k_element]->get_local_function_first_derivative_values (coordinates, k_var, v);
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::get_G_local_matrix (int k_element, Matrix * G_matrix)
{
	 elements[k_element]->get_G_local_matrix (G_matrix);
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::get_M_local_matrix (int k_element, Matrix * G_matrix)
{
	 elements[k_element]->get_M_local_matrix (G_matrix);
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::get_DDD_local_matrix (int k_element, MathVector * DDD)
{
	elements[k_element]->get_DDD (DDD);
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::get_edge_functions (int k_element, int n1, int n2, int * functions)
{
	 elements[k_element]->get_edge_functions (*this, n1, n2, functions);
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::get_def_nodes (int k_element, int * d_nodes) const
{
	 elements[k_element]->get_def_nodes (d_nodes);
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::get_base_nodes (int k_element, int * b_nodes) const
{
	 elements[k_element]->get_base_nodes (b_nodes);
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::get_amount_of_base_nodes (int k_element)
{
	return elements[k_element]->get_amount_of_base_nodes ();
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::get_amount_non_zero_functions (int k_element)
{
	return elements[k_element]->get_amount_non_zero_functions ();
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::belonging_element (int k_node)
{
	int * nodes = NULL;
	int n_nodes;

	for (int k_element = 0, k_element_end = get_n_elements (); k_element < k_element_end; k_element++) // go through elements
	{
		n_nodes = elements[k_element]->get_amount_of_def_nodes ();
		nodes = new int[n_nodes];
		elements[k_element]->get_def_nodes (nodes);
		for (int k = 0; k < n_nodes; k++)
		{
			if (nodes[k] == k_node)
			{
				delete[] nodes;
				return k_element;
			}
		}
		if (nodes != NULL)
			delete[] nodes;
		nodes = NULL;
	}
	if (nodes != NULL)
		delete[] nodes;
	return 0;
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::function_boundary (int k_element, int k_func)
{
	return elements[k_element]->function_boundary (*this, k_func);
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::point_inside (double * coordinates)
{
	int p_element = -1;
	bool found = false;
	for (int i = 0; i < n_elements && !found; i++)
	{
		if (elements[i]->point_inside (*this, coordinates))
		{
			p_element = i;
			found = true;
		}
	}
	return p_element;
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::get_amount_of_folds ()
{
	return 0;
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::get_fold_coordinates (int k_fold, double * c0, double * cN)
{
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::refine (std::vector<int> elements_to_refine)
{
	return 0;
}

template<class Node_Type, class Element_Type>
bool Mesh<Node_Type, Element_Type>::element_inside_area (int k_element, double * c0, double * cN)
{
	int n_nodes = elements[k_element]->get_amount_of_base_nodes ();
	int * base_nodes = new int[n_nodes];
	double ZERO = 1e-12;
	bool found = true;
	// get base nodes of the element
	 elements[k_element]->get_base_nodes (base_nodes);

	double * coordinates = new double[dim];
	for (int i = 0; i < n_nodes && !found; i++)
	{
		nodes[base_nodes[i]].get_coordinates (coordinates);

		// if at least one of them is not inside the area
		for (int k = 0; k < dim && !found; k++)
		{ 
			if (!(c0[k] - ZERO < coordinates[k] && coordinates[k] < cN[k] + ZERO))
				found = false;
		}
	}
	delete[] base_nodes;
	delete[] coordinates;
	return found;
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::function_boundary_edge (int n1, int n2)
{
	int boundary = -1;
	// get coordinates of both nodes
	double * coordinates1 = new double[dim];
	double * coordinates2 =  new double [dim];
	get_node_coordinates (n1, coordinates1);
	get_node_coordinates (n2, coordinates2);

	double * coordinates_middle = new double[dim];

	// get the middle of the vector
	for (int i = 0; i < dim; i++)
	{
		coordinates_middle[i] = coordinates2[i] - 0.5 * (coordinates2[i] - coordinates1[i]);
	}
	// if it is on the boundary, then return the boundary number
	for (int d = 0; d < dim; d++)
	{
		if (abs (coordinates_middle[d] - coord0[d]) < ZERO_boundary_edge)
		{
			boundary = d * dim;
		}

		if (abs (coordinates_middle[d] - coordN[d]) < ZERO_boundary_edge)
		{
			boundary = d * dim + 1;
		}
	}
	delete[] coordinates1;
	delete[] coordinates2;
	delete[] coordinates_middle;
	return boundary;
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::belonging_element (int n1, int n2)
{
	int i_element = -1;
	bool found = false;
	for (int i = 0; i < n_elements && !found; i++)
	{
		if (elements[i]->edge_exists (n1, n2))
		{
			i_element = i;
			found = true;
		}
	}
	return i_element;
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::numerate_edges ()
{
	// for dim = 2
	if (dim == 2)
	{
		std::vector<int> listbeg;
		// get amount of nodes
		for (int i = 0; i < n_nodes; i++)
		{
			listbeg.push_back (-1);
		}
		int listsize = -1;
		std::vector<int> list1;
		std::vector<int> list2;

		int Link1;
		int n_base_nodes;
		int link1, link2; // global numbers of linked nodes
		int iaddr;
		// cycle through elements
		for (int i_element = 0; i_element < n_elements; i_element++)
		{
			// get amount of base nodes of the element
			n_base_nodes = get_amount_of_base_nodes (i_element);
			// cycle through them
			for (int i_node = 0; i_node < n_base_nodes; i_node++)
			{
				// get global number of i_node
				Link1 = get_node_number (i_element, i_node);
				// for following nodes on the element
				for (int j_node = i_node + 1; j_node < n_base_nodes; j_node++)
				{
					link1 = Link1;
					// get global number of j_node
					link2 = get_node_number (i_element, j_node);
					// if node exists
					if (edge_exists (i_element, link1, link2))
					{
						// sort link1 and link2, so link2 < link1 
						// TEST which of > <
						if (link1 < link2)
						{
							link1 = link2;
							link2 = Link1;
						}
						// get place where list for link2 starts
						iaddr = listbeg[link2];
						// if the list hasn't existed
						if (iaddr == -1)
						{
							listsize++;
							listbeg[link2] = listsize;
							list1.push_back (link1);
							list2.push_back (-1);
						}
						else
						{
							// list exists 
							while (list1[iaddr] < link1 && list2[iaddr] > -1) // look for link1 in it
							{
								iaddr = list2[iaddr];
							}
							// next link is bigger than link1
							if (list1[iaddr] > link1)
							{
								listsize++;
								list1.push_back (list1[iaddr]);
								list2.push_back (list2[iaddr]);
								list1[iaddr] = link1;
								list2[iaddr] = listsize;
							}
							else  // end of list
							{
								if (list1[iaddr] < link1) // to skip already existing link
								{
									listsize++;
									list2[iaddr] = listsize;
									list1.push_back (link1);
									list2.push_back (-1);
								}
							}
						}
					}
				}
			}
		}
		listsize++;
		// make portrait
		// set ig, jg
		int * ig = new int[n_nodes + 1];
		int * jg = new int[listsize];

		int i = 0;
		ig[0] = 0;
		// go by global functions
		for (i = 0; i < n_nodes; i++)
		{
			// code for test purposes
			//printf ("\n%i:\t", i);

			ig[i + 1] = ig[i];
			iaddr = listbeg[i];
			// till the end of links for current function
			while (iaddr != -1)
			{
				// code for test purposes
				//printf ("%i ", list1[iaddr]);

				// save the link
				jg[ig[i + 1]] = list1[iaddr];

				// increase ig
				ig[i + 1]++;
				// go to the next link
				iaddr = list2[iaddr];
			}
		}

		// set edges
		if (edges != NULL)
			delete edges;
		edges = new Edge_Structure;
		edges->set_size (n_nodes, listsize);
		edges->set_ig_jg (ig, jg);

		// clear memory
		delete[] ig;
		delete[] jg;
	}
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::get_elements (int k_edge, int * k_element)
{
	int n1, n2;
	int counter = 0;
	
	if (edges != NULL)
	{
		edges->get_edge_nodes (k_edge, &n1, &n2);

		for (int i_element = 0; i_element < n_elements; i_element++)
		{
			if (elements[i_element]->node_belongs (n1) && elements[i_element]->node_belongs (n2))
			{
				k_element[counter] = i_element;
				counter++;
			}
		}
	}
	return counter;
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::get_amount_edges (int k_element)
{
	return elements[k_element]->get_amount_edges ();
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::get_element_edges (int k_element, int * element_edges)
{
	int n_edges = elements[k_element]->get_amount_edges ();

	int ** edges_nodes = new int *[n_edges];
	for (int i = 0; i < n_edges; i++)
	{
		edges_nodes[i] = new int[dim];
	}
	 elements[k_element]->get_edges (edges_nodes);
	for (int i = 0; i < n_edges; i++)
	{
		element_edges[i] = edges->get_edge_number (edges_nodes[i][0], edges_nodes[i][1]);
	}

	for (int i = 0; i < n_edges; i++)
	{
		delete[] edges_nodes[i];
	}
	delete[] edges_nodes;
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::get_area (int k_element) const
{
	return elements[k_element]->get_area();
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::get_area (double * coord) 
{
	int k_element = point_inside (coord);
	return get_area (k_element);
}

template<class Node_Type, class Element_Type>
double Mesh<Node_Type, Element_Type>::get_geometrical_area (int k_element) const
{
	return elements[k_element]->get_geometrical_area ();
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::get_amount_second_condition (int k_element)
{
	return elements[k_element]->get_amount_second_condition ();
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::amount_of_integration_points (int k_element)
{
	return elements[k_element]->amount_of_integration_points ();
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::integration_points (int k_element, double ** points, double * weigths, double * jac)
{
	elements[k_element]->integration_points (points, weigths, jac);
}

template<class Node_Type, class Element_Type>
double Mesh<Node_Type, Element_Type>::get_basis_function_value (int k_element, int k_func, double * coordinates)
{
	return elements[k_element]->get_basis_function_value (k_func, coordinates);
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::get_mass_center (int k_element, double * c)
{
	elements[k_element]->get_mass_center (*this, c);
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::get_order (int k_element)
{
	return elements[k_element]->get_order ();
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::get_function_global_number (int k_element, int k_func)
{
	return elements[k_element]->get_function_global_number (k_func);
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::get_function_by_node (int k_element, int n)
{
	return elements[k_element]->get_function_by_node (n);
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::get_function_by_edge (int k_element, int n1, int n2)
{
	return elements[k_element]->get_function_by_edge (n1, n2);
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::get_n_axis (int * axis) const
{
	for (int i = 0; i < dim; i++)
		axis[i] = n_axis[i];
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::get_n_nodes () const
{
	return n_nodes;
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::get_n_elements () const
{ 
	return n_elements;
}

template<class Node_Type, class Element_Type>
bool Mesh<Node_Type, Element_Type>::need_to_numerate_edges ()
{
	return true;
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::get_elements_by_node (int node, int * elements_by_node)
{
	int counter = 0;
	for (int i = 0; i < n_elements; i++)
	{
		if (elements[i]->node_belongs (node))
		{
			elements_by_node[counter] = i;
			counter++;
		}
	}
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::get_area_dividers (std::vector <std::pair<int, int>> * dividers_nodes)
{
	if (edges != NULL)
	{
		int area;
		int edge_elements[16];
		int n_elements;
		bool add_edge;
		int n1, n2;

		// go by edges
		for (int i = 0, i_end = edges->get_n_entries (); i < i_end; i++)
		{
			n_elements = get_elements (i, edge_elements);
			// if those elements have different areas, add the edge into divider
			area = elements[edge_elements[0]]->get_area ();
			add_edge = false;
			for (int j = 1; j < n_elements; j++)
			{
				if (elements[edge_elements[j]]->get_area () != area)
					add_edge = true;
			}
			if (add_edge)
			{
				edges->get_edge_nodes (i, &n1, &n2);
				dividers_nodes->push_back (std::pair<int, int> (n1, n2));
			}
		}
	}
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::wrap_material (int k_material, double * c0, double * cN)
{
	int n_base_nodes = elements[0]->get_amount_of_base_nodes ();
	int * base_nodes = new int[n_base_nodes];
	int material;

	for (int i = 0; i < dim; i++)
	{
		c0[i] = 1e+15;
		cN[i] = -1e+15;
	}

	double * coordinates = new double[dim];
	for (int i = 0; i < n_elements; i++)
	{
		material = elements[i]->get_area ();

		if (material == k_material)
		{
			elements[i]->get_base_nodes (base_nodes);
			for (int j = 0; j < n_base_nodes; j++)
			{
				nodes[base_nodes[j]].get_coordinates (coordinates);
				for (int k = 0; k < dim; k++)
				{
					if (coordinates[k] < c0[k])
					{
						c0[k] = coordinates[k];
					}
					if (coordinates[k] > cN[k])
					{
						cN[k] = coordinates[k];
					}
				}
			}
		}
	}
	delete[] coordinates;
	delete[] base_nodes;
}

template<class Node_Type, class Element_Type> 
int Mesh<Node_Type, Element_Type>::get_node_number (Node_Type node)
{
	int number = -1;
	bool found = false;
	for (unsigned int i = 0; i < nodes.size () && !found; i++)
	{
		if (node == nodes[i])
		{
			found = true;
			number = i;
		}
	}
	return number;
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::get_node_number (int k_element, int i)
{
	return elements[k_element]->get_node (i);
}

template<class Node_Type, class Element_Type>
Mesh<Node_Type, Element_Type> & Mesh<Node_Type, Element_Type>::operator=(const Mesh<Node_Type, Element_Type> & mesh)
{
	if (this != &mesh)
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
		coord0 = new double[dim];
		for (int i = 0; i < dim; i++)
			coord0[i] = mesh.coord0[i]; // and copy triangular_Mesh's data

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
			elements.push_back (std::make_unique<Element_Type> (*e));

		// copy edges
		if (mesh.edges != NULL)
		{
			delete edges;
			edges = new Edge_Structure;
			edges->Copy (*(mesh.edges));
		}
	}
	return *this;
}

template<class Node_Type, class Element_Type>
bool Mesh<Node_Type, Element_Type>::build_Mesh (char * file_name, bool save, int * N_functions)
{
	printf ("\tscanning mesh data\n");
	input_mesh_data (file_name);
	printf ("\tbuilding a mesh\n");
	if (make_init_Mesh ())
	{
		printf ("\trenumerating mesh\n");
		renumerate ();

		if (need_to_numerate_edges ())
			numerate_edges ();

		*N_functions = numerate_functions ();
		fold_mesh ();

		if (save)
		{
			printf ("\tsaving mesh data for easier reference for a human.\n\tthere is really no need for it, you just don't trust me. *sigh*\n");
			output ();
		}
		return true;
	}
	else
	{
		printf ("ERROR: mesh data is inadeaquate.\ni knew you would fail to give proper instructions\n");
		return false;
	}
}

template<class Node_Type, class Element_Type>
bool Mesh<Node_Type, Element_Type>::build_Mesh (char * file_name, char * file_name_nodes, char * file_name_elements, int * N_functions)
{
	printf ("\tscanning mesh data\n");
	input_mesh_data (file_name);
	printf ("\tbuilding a mesh\n");
	if (make_init_Mesh ())
	{
		printf ("\trenumerating mesh\n");
		//renumerate ();

		if (need_to_numerate_edges ())
			numerate_edges ();

		*N_functions = numerate_functions ();
		fold_mesh ();

		printf ("\tsaving mesh data for easier reference for a human.\n\tthere is really no need for it, you just don't trust me. *sigh*\n");
		output (file_name_nodes, file_name_elements);

		return true;
	}
	else
	{
		printf ("ERROR: mesh data is inadeaquate.\ni knew you would fail to give proper instructions\n");
		return false;
	}
}

template<class Node_Type, class Element_Type>
bool Mesh<Node_Type, Element_Type>::build_Mesh (char * file_name_nodes, char * file_name_elements, int * N_functions)
{
	for (int i = 0; i < dim; i++)
	{
		coord0[i] = 1e+15;
		coordN[i] = -1e+15;
	}

	// get nodes
	FILE * file = fopen (file_name_nodes, "r");
	double * coordinates = new double[dim];
	Node_Type node;
	while (!feof (file))
	{
		for (int i = 0; i < dim; i++)
		{
			fscanf (file, "%lf", &coordinates[i]);
			if (coordinates[i] > coordN[i])
			{
				coordN[i] = coordinates[i];
			}
			if (coordinates[i] < coord0[i])
			{
				coord0[i] = coordinates[i];
			}
		}
		node.set_coordinates (coordinates);
		nodes.push_back (node);
	}
	nodes.pop_back ();
	n_nodes = (int)nodes.size ();
	fclose (file);
	printf ("\tgot %i nodes\n", n_nodes);

	std::set <int> material_areas;
	
	// get elements
	file = fopen (file_name_elements, "r");

	while (!feof (file))
	{
		std::unique_ptr <Element_Type> element = std::make_unique<Element_Type> ();
		int n_def_nodes = element->get_amount_of_def_nodes ();
		int area;
		int * element_nodes = new int[n_def_nodes];

		for (int i = 0; i < n_def_nodes; i++)
		{
			fscanf (file, "%i", &element_nodes[i]);
		}
		fscanf (file, "%i", &area);
		material_areas.insert (area);

		element->set_base_nodes (element_nodes);
		element->set_element (element_nodes);
		element->set_area (area);
		elements.push_back (std::move(element));
		delete[] element_nodes;

		// scanfs order which is not needed typically, except mixed_tr_mesh
		fscanf (file, "%i", &area);
	}
	elements.pop_back ();
	n_elements = (int)elements.size ();
	n_material_areas = (int)material_areas.size ();
	fclose (file);
	printf ("\tgot %i elements\n", n_elements);
	printf ("\tgot %i areas\n", n_material_areas);
		
	if (need_to_numerate_edges ())
		numerate_edges ();
	*N_functions = numerate_functions ();
	fold_mesh ();

	delete[] coordinates;
	return true;
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::input_mesh_data (char * file_name)
{
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::input_mesh_data (char * file_name, char * file_density)
{
}

template<class Node_Type, class Element_Type>
bool Mesh<Node_Type, Element_Type>::make_init_Mesh ()
{
	return false;
}

template<class Node_Type, class Element_Type>
bool Mesh<Node_Type, Element_Type>::make_init_Mesh (char * file_density)
{
	return false;
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::renumerate ()
{
}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::numerate_functions ()
{
	for (size_t i = 0, i_end = elements.size (); i < i_end; i++)
	{
		elements[i]->set_global_functions ();
	}

	return n_nodes;
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::fold_mesh ()
{

}

template<class Node_Type, class Element_Type>
int Mesh<Node_Type, Element_Type>::find_maximal_node ()
{
	double l, r;
	l = -1e+15;
	int res = 0;

	for (unsigned int i = 0; i < nodes.size (); i++)
	{
		r = nodes[i].distance_origin ();
		if (r > l)
		{
			res = (int)i;
			l = r;
		}
	}
	return res;
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::output ()
{
}

template<class Node_Type, class Element_Type>
void Mesh<Node_Type, Element_Type>::output (char * file_nodes, char * file_triangles)
{
}
