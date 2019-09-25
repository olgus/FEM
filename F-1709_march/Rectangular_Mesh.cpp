#include "Rectangular_Mesh.h"

void Rectangular_Mesh::input_mesh_data (char * file_name)
{
	FILE * file = fopen (file_name, "r");

	fscanf (file, "%i %i", &n_axis[0], &n_axis[1]);
	fscanf (file, "%i", &lvl);

	n_axis[0] *= (int)round (pow (2.0, lvl));
	n_axis[1] *= (int)round (pow (2.0, lvl));

	tetra_nodes = new Point_2D *[2];
	for (int i = 0; i < 2; i++)
		tetra_nodes[i] = new Point_2D[2];

	// get nodes
	double c0, c1;
	for (int j = 0; j < 2; j++)
	{
		for (int i = 0; i < 2; i++)
		{
			fscanf (file, "%lf %lf", &c0, &c1);
			tetra_nodes[i][j].set_dim (2);
			tetra_nodes[i][j].set_point (c0, c1);
		}
	}
	// reset boundaries' values
	coord0[0] = tetra_nodes[0][0].X ();
	coord0[1] = tetra_nodes[0][0].Y ();
	coordN[0] = tetra_nodes[1][1].X ();
	coordN[1] = tetra_nodes[1][1].Y ();

	fscanf (file, "%i", &material);
	fscanf (file, "%lf %lf", &coef[0], &coef[1]);
	fscanf (file, "%i %i", &direc[0], &direc[1]);

	fclose (file);
}

bool Rectangular_Mesh::make_init_Mesh ()
{
	// go from the bottom
	double c[2]; // x, y coordinates of the new node
	double unit[2]; // unit vector on the side
	double GPS[2]; // geometry progression sum
	double L[2]; // full length of the side
	double l[2]; // current length
	double q[2]; // current coefficient
	double l0[2]; // starting length
	Node_2D node; // new node
	int rect_base_nodes[4]; // nodes that define the rectange
	int cur_node;

	// X
	L[0] = sqrt (pow (tetra_nodes[0][0].X () - tetra_nodes[1][1].X (), 2.0));
	// set q
	q[0] = coef[0];
	if (direc[0] == -1)
		q[0] = 1.0 / q[0];
	// get geometric sum
	if (fabs (q[0] - 1.0) < ZERO_rectangular_mesh) // or if coef == 1, just the amount of sections
		GPS[0] = n_axis[0];
	else
		GPS[0] = (1.0 - pow (q[0], n_axis[0])) / (1.0 - q[0]);
	l0[0] = L[0] / GPS[0];
	unit[0] = (tetra_nodes[1][1].X () - tetra_nodes[0][0].X ()) / L[0];

	// Y
	L[1] = sqrt (pow (tetra_nodes[0][0].Y () - tetra_nodes[1][1].Y (), 2.0));
	// set q
	q[1] = coef[1];
	if (direc[1] == -1)
		q[1] = 1.0 / q[1];
	// get geometric sum
	if (fabs (q[1] - 1.0) < ZERO_rectangular_mesh) // or if coef == 1, just the amount of sections
		GPS[1] = n_axis[1];
	else
		GPS[1] = (1.0 - pow (q[1], n_axis[1])) / (1.0 - q[1]);
	l0[1] = L[1] / GPS[1];
	unit[1] = (tetra_nodes[1][1].Y () - tetra_nodes[0][0].Y ()) / L[1];
	
	l[0] = l0[0];
	l[1] = l0[1];

	c[1] = tetra_nodes[0][0].Y ();
	for (int j = 0; j < n_axis[1]; j++) // go through sections
	{		
		// start with x0 and lenght0
		c[0] = tetra_nodes[0][0].X ();
		l[0] = l0[0]; 
		for (int i = 0; i < n_axis[0]; i++)
		{
			// calculate down left node
			node.set_coordinates (c);
			// add it
			if (std::find (nodes.begin (), nodes.end (), node) == nodes.end ())
			{
				nodes.push_back (node);
				rect_base_nodes[0] = n_nodes;
				n_nodes++;
			}
			else
			{
				cur_node = get_node_number (node);
				rect_base_nodes[0] = cur_node;
			}
			// calculate down right node
			c[0] += l[0];
			node.set_coordinates (c);
			// add it
			if (std::find (nodes.begin (), nodes.end (), node) == nodes.end ())
			{
				nodes.push_back (node);
				rect_base_nodes[1] = n_nodes;
				n_nodes++;
			}
			else
			{
				cur_node = get_node_number (node);
				rect_base_nodes[1] = cur_node;
			}
			// calculate up left node
			c[0] -= l[0];
			c[1] += l[1];
			node.set_coordinates (c);
			// add it
			if (std::find (nodes.begin (), nodes.end (), node) == nodes.end ())
			{
				nodes.push_back (node);
				rect_base_nodes[2] = n_nodes;
				n_nodes++;
			}
			else
			{
				cur_node = get_node_number (node);
				rect_base_nodes[2] = cur_node;
			}
			// calculate up right node
			c[0] += l[0];
			node.set_coordinates (c);
			// add it
			if (std::find (nodes.begin (), nodes.end (), node) == nodes.end ())
			{
				nodes.push_back (node);
				rect_base_nodes[3] = n_nodes;
				n_nodes++;
			}
			else
			{
				cur_node = get_node_number (node);
				rect_base_nodes[3] = cur_node;
			}
			c[1] -= l[1];
			
			// set base nodes
			std::unique_ptr <Rectangle_Element> rectangle_Element = std::make_unique<Rectangle_Element> ();
			//auto rectangle_Element = std::make_unique <Rectangle_Element> ();
			rectangle_Element->set_base_nodes (rect_base_nodes);
			rectangle_Element->set_element (rect_base_nodes);
			rectangle_Element->set_area (material);
			// add element
			elements.push_back (std::move (rectangle_Element));

			n_elements++;

			if (fabs (q[0] - 1.0) > ZERO_rectangular_mesh) // calculate new lenght by x
				l[0] = l[0] * q[0];

		}
		c[1] += l[1];
		// calculate new lenght by y
		l[1] = l[1] * q[1];
	}
	return true;
}

void Rectangular_Mesh::output ()
{
	FILE * log = fopen ("Result Files Extra//log_mesh.txt", "w");

	fprintf (log, "%i %i\n", n_nodes, n_elements);
	printf ("Nodes:\t%i\tElements:\t%i\n", n_nodes, n_elements);

	fclose (log);

	log = fopen ("MeshData\\Nodes_rect.txt", "w");

	double coordinates[2];
	for (int i = 0; i < n_nodes; i++)
	{
		nodes[i].get_coordinates (coordinates);
		fprintf (log, "%.16lf %.16lf\n", coordinates[0], coordinates[1]);
	}

	fclose (log);

	log = fopen ("MeshData\\Rectangle_Elements.txt", "w");

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

Rectangular_Mesh::Rectangular_Mesh ()
{
	dim = 2;

	coord0 = new double[dim];
	coordN = new double[dim];
	n_axis = new int[dim];
}

Rectangular_Mesh::Rectangular_Mesh (const Rectangular_Mesh & rectangular_mesh)
{
	// if memory has been allocated, free it
	if (coord0 != NULL)
		delete[] coord0;
	if (coordN != NULL)
		delete[] coordN;
	if (n_axis != NULL)
		delete[] n_axis;

	// reset dimentionality
	dim = rectangular_mesh.dim;

	// allocate the memory
	// and copy triangular_Mesh's data
	coord0 = new double[dim];
	for (int i = 0; i < dim; i++)
		coord0[i] = rectangular_mesh.coord0[i];

	coordN = new double[dim];
	for (int i = 0; i < dim; i++)
		coordN[i] = rectangular_mesh.coordN[i];

	n_axis = new int[dim];
	for (int i = 0; i < dim; i++)
		n_axis[i] = rectangular_mesh.n_axis[i];

	// clear nodes and copy triangular_Mesh's nodes
	nodes.clear ();
	n_nodes = rectangular_mesh.n_nodes;
	nodes.insert (nodes.begin (), rectangular_mesh.nodes.begin (), rectangular_mesh.nodes.end ());

	// clear elements and copy triangular_Mesh's elements
	elements.clear ();
	n_elements = rectangular_mesh.n_elements;
	for (const auto& e : rectangular_mesh.elements)
		elements.push_back (std::make_unique<Rectangle_Element> (*e));
}

Rectangular_Mesh::~Rectangular_Mesh ()
{
}

void Rectangular_Mesh::make_uniform_mesh (int n1, int n2, double * c0, double * cN, int Lvl)
{
	for (int i = 0; i < 2; i++)
	{
		direc[i] = 1;
		coef[i] = 1.0;
		coord0[i] = c0[i];
		coordN[i] = cN[i];
	}
	material = 1;
	lvl = Lvl;

	n_axis[0] = n1;
	n_axis[1] = n2;

	n_axis[0] *= (int)round (pow (2.0, lvl));
	n_axis[1] *= (int)round (pow (2.0, lvl));

	tetra_nodes = new Point_2D *[2];
	for (int i = 0; i < 2; i++)
		tetra_nodes[i] = new Point_2D[2];

	tetra_nodes[0][0].set_dim (2);
	tetra_nodes[0][0].set_point (c0[0], c0[1]);
	tetra_nodes[0][1].set_dim (2);
	tetra_nodes[0][1].set_point (cN[0], c0[1]);
	tetra_nodes[1][0].set_dim (2);
	tetra_nodes[1][0].set_point (c0[0], cN[1]);
	tetra_nodes[1][1].set_dim (2);
	tetra_nodes[1][1].set_point (cN[0], cN[1]);
}
