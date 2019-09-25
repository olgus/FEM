#include "Cubic_Mesh.h"

void Cubic_Mesh::input_mesh_data (char * file_name)
{
	FILE * file = fopen (file_name, "r");
	for (int i = 0; i < 3; i++)
	{
		fscanf (file, "%lf", &coord0[i]);
		fscanf (file, "%lf", &coordN[i]);
		fscanf (file, "%i", &n_axis[i]);
	}
	fclose (file);
}

bool Cubic_Mesh::make_init_Mesh ()
{
	double h[3];
	for (int i = 0; i < 3; i++)
	{
		h[i] = (coordN[i] - coord0[i]) / n_axis[i];
	}

	// nodes
	Node_3D node;
	double coordinates[3];
	int counter = 0;
	for (int k = 0; k < n_axis[2] + 1; k++)
	{
		coordinates[2] = coord0[2] + h[2] * k;
		for (int j = 0; j < n_axis[1] + 1; j++)
		{
			coordinates[1] = coord0[1] + h[1] * j;
			for (int i = 0; i < n_axis[0] + 1; i++)
			{
				coordinates[0] = coord0[0] + h[0] * i;
				node.set_coordinates (coordinates);
				nodes.push_back (node);
				counter++;
			}
		}
	}
	n_nodes = counter;

	// elements
	counter = 0;
	int nodes_element[8];
	for (int k = 0; k < n_axis[2]; k++)
	{
		for (int j = 0; j < n_axis[1]; j++)
		{
			for (int i = 0; i < n_axis[0]; i++)
			{
				nodes_element[0] = i + j * (n_axis[0] + 1) + k * (n_axis[0] + 1) * (n_axis[1] + 1);
				nodes_element[1] = i + 1 + j * (n_axis[0] + 1) + k * (n_axis[0] + 1) * (n_axis[1] + 1);
				nodes_element[2] = i + (j + 1) * (n_axis[0] + 1) + k * (n_axis[0] + 1) * (n_axis[1] + 1);
				nodes_element[3] = i + 1 + (j + 1) * (n_axis[0] + 1) + k * (n_axis[0] + 1) * (n_axis[1] + 1);
				nodes_element[4] = i + j * (n_axis[0] + 1) + (k + 1) * (n_axis[0] + 1) * (n_axis[1] + 1);
				nodes_element[5] = i + 1 + j * (n_axis[0] + 1) + (k + 1) * (n_axis[0] + 1) * (n_axis[1] + 1);
				nodes_element[6] = i + (j + 1) * (n_axis[0] + 1) + (k + 1) * (n_axis[0] + 1) * (n_axis[1] + 1);
				nodes_element[7] = i + 1 + (j + 1) * (n_axis[0] + 1) + (k + 1) * (n_axis[0] + 1) * (n_axis[1] + 1);
				
				std::unique_ptr <Cube> cube = std::make_unique<Cube> ();
				cube->set_nodes (nodes_element);
				cube->set_base_nodes (nodes_element);
				cube->set_area (1);
				elements.push_back (std::move (cube));
				counter++;
			}
		}
	}
	n_elements = counter;
	// prepare
	for (int i = 0; i < n_elements; i++)
	{
		elements[i]->prepare (*this);
	}
	return true;
}

void Cubic_Mesh::output ()
{
	FILE * log = fopen ("MeshData\\Nodes_cube.txt", "w");

	double coordinates[3];
	for (int i = 0; i < n_nodes; i++)
	{
		nodes[i].get_coordinates (coordinates);
		fprintf (log, "%.16lf %.16lf %.16lf\n", coordinates[0], coordinates[1], coordinates[2]);
	}
	fclose (log);

	log = fopen ("MeshData\\Cubes.txt", "w");

	int nodes[8];
	for (int i = 0; i < n_elements; i++)
	{
		elements[i]->get_def_nodes (nodes);
		fprintf (log, "%i %i %i %i %i %i %i %i\n", nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7]);
	}

	fclose (log);
}

bool Cubic_Mesh::check_neighbouring_state (int e1, int e2)
{
	// 4 out of 8 nodes are the same?
	if (e1 != e2)
	{
		int overlap = 0;
		int e1_nodes[8];
		get_base_nodes (e1, e1_nodes);
		int e2_nodes[8];
		get_base_nodes (e2, e2_nodes);
		bool found;
		for (int i = 0; i < 8; i++)
		{
			found = false;
			for (int j = 0; j < 8 && !found; j++)
			{
				if (e1_nodes[i] == e2_nodes[j])
				{
					overlap++;
					found = true;
				}
			}
		}
		if (overlap == 4)
			return true;
	}
	return false;
}

bool Cubic_Mesh::need_to_numerate_edges ()
{
	return false;
}

bool Cubic_Mesh::build_Mesh (double * C0, double * CN, int * N_axis)
{
	for (int i = 0; i < 3; i++)
	{
		coord0[i] = C0[i];
		coordN[i] = CN[i];
		n_axis[i] = N_axis[i];
	}

	printf ("\tbuilding a mesh\n");
	make_init_Mesh ();
	output ();

	return true;
}

void Cubic_Mesh::reset_areas (char * file_areas)
{
	int N_areas = 0;
	FILE * file = fopen (file_areas, "r");
	fscanf (file, "%i", &N_areas);

	int area;
	double ac0[3], acN[3];
	double mass_center[3];

	// set everything to 0 for covinience
	for (int j = 0; j < n_elements; j++)
	{
		elements[j]->set_area (0);
	}

	for (int i = 0; i < N_areas; i++)
	{
		fscanf (file, "%i", &area);
		fscanf (file, "%lf %lf %lf", &ac0[0], &ac0[1], &ac0[2]);
		fscanf (file, "%lf %lf %lf", &acN[0], &acN[1], &acN[2]);
		// if mass center of the element is inside area, change it's area value
		for (int j = 0; j < n_elements; j++)
		{
			elements[j]->get_mass_center (*this, mass_center);
			int counter = 0;
			for (int k = 0; k < 3; k++)
			{
				if (ac0[k] - 1e-10 < mass_center[k] && mass_center[k] < acN[k] + 1e-10)
				{
					counter++;
				}
			}
			if (counter == 3)
			{
				elements[j]->set_area (area);
			}
		}

	}
	fclose (file);
}

Cubic_Mesh::Cubic_Mesh ()
{
	dim = 3;

	coord0 = new double[dim];
	coordN = new double[dim];
	n_axis = new int[dim];
}

Cubic_Mesh::Cubic_Mesh (const Cubic_Mesh & cubic_Mesh)
{
	// if memory has been allocated, free it
	if (coord0 != NULL)
		delete[] coord0;
	if (coordN != NULL)
		delete[] coordN;
	if (n_axis != NULL)
		delete[] n_axis;

	// reset dimentionality
	dim = cubic_Mesh.dim;

	// allocate the memory
	// and copy triangular_Mesh's data
	coord0 = new double[dim];
	for (int i = 0; i < dim; i++)
		coord0[i] = cubic_Mesh.coord0[i];

	coordN = new double[dim];
	for (int i = 0; i < dim; i++)
		coordN[i] = cubic_Mesh.coordN[i];

	n_axis = new int[dim];
	for (int i = 0; i < dim; i++)
		n_axis[i] = cubic_Mesh.n_axis[i];

	// clear nodes and copy triangular_Mesh's nodes
	nodes.clear ();
	n_nodes = cubic_Mesh.n_nodes;
	nodes.insert (nodes.begin (), cubic_Mesh.nodes.begin (), cubic_Mesh.nodes.end ());

	// clear elements and copy triangular_Mesh's elements
	elements.clear ();
	n_elements = cubic_Mesh.n_elements;
	for (const auto& e : cubic_Mesh.elements)
		elements.push_back (std::make_unique<Cube> (*e));
}

Cubic_Mesh::~Cubic_Mesh ()
{
}
