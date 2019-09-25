#include "Prismatic_Mesh.h"

void Prismatic_Mesh::input_mesh_data (char * file_name)
{
	char file_name_mesh[512];
	strcpy (file_name_mesh, file_name);
	strcat (file_name_mesh, "xy.txt");
	mesh_triangular = new Triangular_Mesh ();

	int n_functions;
	mesh_triangular->build_Mesh (file_name_mesh, false, &n_functions);

	strcpy (file_name_mesh, file_name);
	strcat (file_name_mesh, "z.txt");
	mesh_1d = new Mesh_1D_L1 ();
	mesh_1d->build_Mesh (file_name_mesh, false, &n_functions);

	// save n_axis
	mesh_1d->get_n_axis (n_axis);
	n_axis[2] = n_axis[0];
	mesh_triangular->get_n_axis (n_axis);

	// save coord0
	mesh_1d->get_0_boundaries (coord0);
	coord0[2] = coord0[0];
	mesh_triangular->get_0_boundaries (coord0);

	// save coordN
	mesh_1d->get_N_boundaries (coordN);
	coordN[2] = coordN[0];
	mesh_triangular->get_N_boundaries (coordN);
}

bool Prismatic_Mesh::make_init_Mesh ()
{
	int n_crosssection = mesh_triangular->get_n_nodes (); // amount of nodes in triangular mesh
	int n_levels = mesh_1d->get_n_nodes (); // amount of nodes in 1d_mesh

	double coordinates_1d[1];
	double coordinates_3d[3];
	Node_3D node;

	// form nodes
	// iterating through Z using 1d_mesh
	for (int j = 0; j < n_levels; j++)
	{
		mesh_1d->get_node_coordinates (j, coordinates_1d);
		// go through triangles
		for (int i = 0; i < n_crosssection; i++)
		{
			mesh_triangular->get_node_coordinates (i, coordinates_3d);
			coordinates_3d[2] = coordinates_1d[0];
			node.set_coordinates (coordinates_3d);
			nodes.push_back (node);
		}
	}
	n_nodes = n_crosssection * n_levels; 

	// form elements
	int prism_nodes[6];
	int n_triangles = mesh_triangular->get_n_elements (); // amount of elements in triangular mesh
	int n_z = mesh_1d->get_n_elements (); // amount of elements in triangular mesh
	for (int i = 0; i < n_triangles; i++)
	{
		mesh_triangular->get_def_nodes (i, prism_nodes);
		// go through Z
		for (int j = 0; j < n_z; j++)
		{
			// add top level
			for (int k = 0; k < 3; k++)
			{
				prism_nodes[k + 3] = prism_nodes[k] + n_crosssection;
			}

			// add prism into elements
			std::unique_ptr <Prism> prism = std::make_unique<Prism> ();
			prism->set_area (mesh_triangular->get_area (i));
			prism->set_base_nodes (prism_nodes);
			prism->set_element (prism_nodes);
			elements.push_back (std::move (prism));

			// save top level for the next prism
			for (int k = 0; k < 3; k++)
			{
				prism_nodes[k] = prism_nodes[k + 3];
			}
		}
	}
	n_elements = n_triangles * n_z;
	return true;
}

void Prismatic_Mesh::output ()
{
	FILE * log = fopen ("Result Files Extra//log_mesh.txt", "w");

	// print mesh data
	fprintf (log, "%i %i\n", n_nodes, n_elements);
	printf ("Nodes:\t%i\tElements:\t%i\n", n_nodes, n_elements);

	fclose (log);

	log = fopen ("MeshData\\Nodes_prism.txt", "w");

	double coordinates[3];
	for (int i = 0; i < n_nodes; i++)
	{
		nodes[i].get_coordinates (coordinates);
		fprintf (log, "%.16lf %.16lf %.16lf\n", coordinates[0], coordinates[1], coordinates[2]);
	}

	fclose (log);

	log = fopen ("MeshData\\Prisms.txt", "w");

	int nodes[6];
	for (int i = 0; i < n_elements; i++)
	{
		elements[i]->get_def_nodes (nodes);
		fprintf (log, "%i %i %i %i %i %i %i\n", nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], elements[i]->get_area ());
	}

	fclose (log);
}

bool Prismatic_Mesh::need_to_numerate_edges ()
{
	return false;
}

Prismatic_Mesh::Prismatic_Mesh ()
{
	dim = 3;

	coord0 = new double[dim];
	coordN = new double[dim];
	n_axis = new int[dim];
}

Prismatic_Mesh::Prismatic_Mesh (const Prismatic_Mesh & prismatic_Mesh)
{
	// if memory has been allocated, free it
	if (coord0 != NULL)
		delete[] coord0;
	if (coordN != NULL)
		delete[] coordN;
	if (n_axis != NULL)
		delete[] n_axis;

	// reset dimentionality
	dim = prismatic_Mesh.dim;

	// allocate the memory
	// and copy triangular_Mesh's data
	coord0 = new double[dim];
	for (int i = 0; i < dim; i++)
		coord0[i] = prismatic_Mesh.coord0[i];

	coordN = new double[dim];
	for (int i = 0; i < dim; i++)
		coordN[i] = prismatic_Mesh.coordN[i];

	n_axis = new int[dim];
	for (int i = 0; i < dim; i++)
		n_axis[i] = prismatic_Mesh.n_axis[i];

	// clear nodes and copy triangular_Mesh's nodes
	nodes.clear ();
	n_nodes = prismatic_Mesh.n_nodes;
	nodes.insert (nodes.begin (), prismatic_Mesh.nodes.begin (), prismatic_Mesh.nodes.end ());

	// clear elements and copy triangular_Mesh's elements
	elements.clear ();
	n_elements = prismatic_Mesh.n_elements;
	for (const auto& e : prismatic_Mesh.elements)
		elements.push_back (std::make_unique<Prism> (*e));
}

Prismatic_Mesh::~Prismatic_Mesh ()
{
}
