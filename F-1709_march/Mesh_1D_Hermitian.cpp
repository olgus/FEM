#include "Mesh_1D_Hermitian.h"

int Mesh_1D_Hermitian::numerate_functions ()
{
	int functions[4];
	for (size_t i = 0, i_end = elements.size (); i < i_end; i++)
	{
		functions[0] = (int)(i * 2);
		functions[1] = (int)(i * 2 + 1);
		functions[2] = (int)((i + 1) * 2);
		functions[3] = (int)((i + 1) * 2 + 1);

		elements[i]->set_global_functions (functions);
	}
	return n_nodes * 2;
}

void Mesh_1D_Hermitian::input_mesh_data (char * file_name)
{
	FILE * file = fopen (file_name, "r");

	fscanf (file, "%i", &n_axis[0]);
	fscanf (file, "%i", &lvl);

	// get sections' boundaries values
	boundaries = new double[n_axis[0] + 1];
	for (int i = 0; i < n_axis[0] + 1; i++)
	{
		fscanf (file, "%lf", &boundaries[i]);
	}

	// get amount of subsections
	n_sections = new int[n_axis[0]];
	for (int i = 0; i < n_axis[0]; i++)
	{
		fscanf (file, "%i", &n_sections[i]);
		n_sections[i] *= (int)round (pow (2.0, lvl));
	}

	// get coefficients of subsections
	coef = new double[n_axis[0]];
	for (int i = 0; i < n_axis[0]; i++)
	{
		fscanf (file, "%lf", &coef[i]);
	}

	// get directions of subsections
	int direc;
	for (int i = 0; i < n_axis[0]; i++)
	{
		fscanf (file, "%i", &direc);
		if (direc == -1)
			coef[i] = 1.0 / coef[i];
	}

	// get materials 
	materials = new int[n_axis[0]];
	for (int i = 0; i < n_axis[0]; i++)
	{
		fscanf (file, "%i", &materials[i]);
	}
	fclose (file);
}

bool Mesh_1D_Hermitian::make_init_Mesh ()
{
	Node_1D node;
	double coordinates[1];

	double l0; // first lenght
	double lcur; // length of current segment
	double GPS; // geometric sum
	double L; // length of a section
	int el_nodes[2];

	n_elements = 0;
	n_nodes = 0;

	// push first node
	coordinates[0] = boundaries[0];
	coord0[0] = boundaries[0];
	// push new node
	node.set_coordinates (coordinates);
	nodes.push_back (node);

	// go through sections
	for (int i = 0; i < n_axis[0]; i++)
	{
		// get GPS
		if (fabs (coef[i] - 1.0) < ZERO_mesh_1D_L1)
		{
			GPS = n_sections[i];
		}
		else
		{
			GPS = (1.0 - pow (coef[i], n_sections[i])) / (1.0 - coef[i]); // TEST
		}
		// get L
		L = (boundaries[i + 1] - boundaries[i]);
		// get l0
		l0 = L / GPS;
		// go through subsections
		// set starting node as boundary[i]
		coordinates[0] = boundaries[i];
		for (int j = 0; j < n_sections[i]; j++)
		{
			// current length
			lcur = l0 * pow (coef[i], j);
			// add that lenght to previous node
			coordinates[0] += lcur;
			// push new node
			node.set_coordinates (coordinates);
			nodes.push_back (node);
			n_nodes++;

			// add new element
			std::unique_ptr <Element_1D_Hermitian> element = std::make_unique<Element_1D_Hermitian> ();
			el_nodes[0] = n_nodes - 1;
			el_nodes[1] = n_nodes;
			element->set_area (materials[i]);
			element->set_base_nodes (el_nodes);
			element->set_element (el_nodes);
			elements.push_back (std::move (element));
		}

		n_elements += n_sections[i];
	}
	n_nodes++;

	coordN[0] = boundaries[n_axis[0]];

	return true;
}

bool Mesh_1D_Hermitian::need_to_numerate_edges ()
{
	return false;
}

bool Mesh_1D_Hermitian::build_Mesh (std::vector<double> coordinates, int * n_functions)
{
	// coordinates are nodes
	Node_1D node;
	double c[1];
	int el_nodes[2];

	// push first node
	c[0] = coordinates[0];
	coord0[0] = c[0];
	// push new node
	node.set_coordinates (c);
	nodes.push_back (node);

	n_elements = n_nodes = 0;

	for (size_t i = 1, i_end = coordinates.size (); i < i_end; i++)
	{
		// next node is next data point
		c[0] = coordinates[i];
		// make a node
		node.set_coordinates (c);
		// push it
		nodes.push_back (node);
		n_nodes++;

		// add new element
		std::unique_ptr <Element_1D_Hermitian> element = std::make_unique<Element_1D_Hermitian> ();

		el_nodes[0] = n_nodes - 1;
		el_nodes[1] = n_nodes;

		element->set_area (1);
		element->set_base_nodes (el_nodes);
		element->set_element (el_nodes);
		elements.push_back (std::move (element));
		n_elements++;
	}

	coordN[0] = c[0];

	*n_functions = numerate_functions ();
	return true;
}

bool Mesh_1D_Hermitian::build_Mesh (double c0, double cN, int n_x, int * n_functions)
{
	n_axis[0] = n_x;
	coord0[0] = c0;
	coordN[0] = cN;

	int material = 1;
	double h = (cN - c0) / n_x;
	double ce0[1], ceN[1];
	ce0[0] = coord0[0];
	Node_1D node;
	node.set_coordinates (ce0);
	nodes.push_back (node);
	int cur_node = 0;
	int el_nodes[2];

	for (int i = 0; i < n_axis[0]; i++)
	{
		ce0[0] = coord0[0] + i * h;
		ceN[0] = coord0[0] + (i + 1) * h;

		// push new node
		node.set_coordinates (ceN);
		nodes.push_back (node);
		cur_node++;

		// add new element
		std::unique_ptr <Element_1D_Hermitian> element = std::make_unique<Element_1D_Hermitian> ();
		el_nodes[0] = cur_node - 1;
		el_nodes[1] = cur_node;
		element->set_area (material);
		element->set_base_nodes (el_nodes);
		element->set_element (el_nodes);
		elements.push_back (std::move (element));
	}

	n_nodes = cur_node + 1;
	n_elements = (int)elements.size ();

	*n_functions = numerate_functions ();
	return true;
}

Mesh_1D_Hermitian::Mesh_1D_Hermitian ()
{
	dim = 1;

	coord0 = new double[dim];
	coordN = new double[dim];
	n_axis = new int[dim];

	materials = NULL;
	n_sections = NULL;
	coef = NULL;
	boundaries = NULL;
}

Mesh_1D_Hermitian::Mesh_1D_Hermitian (const Mesh_1D_Hermitian & mesh)
{
	// if memory has been allocated, free it
	if (coord0 != NULL)
		delete[] coord0;
	if (coordN != NULL)
		delete[] coordN;
	if (n_axis != NULL)
		delete[] n_axis;
	if (n_sections != NULL)
		delete[] coord0;
	if (coef != NULL)
		delete[] coordN;
	if (boundaries != NULL)
		delete[] n_axis;

	// reset dimentionality
	dim = mesh.dim;

	// allocate the memory
	// and copy mesh's data
	coord0 = new double[dim];
	if (mesh.coord0 != NULL)
		for (int i = 0; i < dim; i++)
			coord0[i] = mesh.coord0[i];

	coordN = new double[dim];
	if (mesh.coord0 != NULL)
		for (int i = 0; i < dim; i++)
			coordN[i] = mesh.coordN[i];

	n_axis = new int[dim];
	if (mesh.n_axis != NULL)
		for (int i = 0; i < dim; i++)
			n_axis[i] = mesh.n_axis[i];

	n_sections = new int[n_axis[0]];
	if (mesh.n_sections != NULL)
		for (int i = 0; i < n_axis[0]; i++)
			n_sections[i] = mesh.n_sections[i];

	coef = new double[n_axis[0]];
	if (mesh.coef != NULL)
		for (int i = 0; i < n_axis[0]; i++)
			coef[i] = mesh.coef[i];

	boundaries = new double[n_axis[0] + 1];
	if (mesh.boundaries != NULL)
		for (int i = 0; i < n_axis[0] + 1; i++)
			boundaries[i] = mesh.boundaries[i];

	materials = new int[n_axis[0]];
	if (mesh.materials != NULL)
		for (int i = 0; i < n_axis[0]; i++)
			materials[i] = mesh.materials[i];

	// clear nodes and copy triangular_Mesh's nodes
	nodes.clear ();
	n_nodes = mesh.n_nodes;
	nodes.insert (nodes.begin (), mesh.nodes.begin (), mesh.nodes.end ());

	// clear elements and copy triangular_Mesh's elements
	elements.clear ();
	n_elements = mesh.n_elements;
	for (const auto& e : mesh.elements)
		elements.push_back (std::make_unique<Element_1D_Hermitian> (*e));
}

Mesh_1D_Hermitian::~Mesh_1D_Hermitian ()
{
	delete[] materials;
	delete[] n_sections;
	delete[] coef;
	delete[] boundaries;
}

void Mesh_1D_Hermitian::output ()
{
	printf ("Nodes:\t%i\tElements:\t%i\n", n_nodes, n_elements);
	FILE * log = fopen ("MeshData\\Nodes_1D.txt", "w");

	double coordinates[1];
	for (int i = 0; i < n_nodes; i++)
	{
		nodes[i].get_coordinates (coordinates);
		fprintf (log, "%.16lf\n", coordinates[0]);
	}

	fclose (log);

	log = fopen ("MeshData\\Elements_1D.txt", "w");

	int * nodes;
	int j_end;
	for (int i = 0; i < n_elements; i++)
	{
		j_end = elements[i]->get_amount_of_def_nodes ();
		nodes = new int[j_end];
		elements[i]->get_def_nodes (nodes);
		for (int j = 0; j < j_end; j++)
		{
			fprintf (log, "%i ", nodes[j]);
		}
		fprintf (log, "%i\n", elements[i]->get_area ());
	}

	fclose (log);
}
