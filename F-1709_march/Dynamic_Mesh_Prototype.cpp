#include "Dynamic_Mesh_Prototype.h"

void Dynamic_Mesh_Builder::build_task_matrix ()
{
}

void Dynamic_Mesh_Builder::fix_boundary_nodes ()
{
}

double Dynamic_Mesh_Builder::monitor_function (double * coordinates)
{
	return 0.0;
}

void Dynamic_Mesh_Builder::adapt_phase_front (Task_pointer * task, Mesh_Prototype * mesh)
{
	double * v = new double[mesh->get_dimentionality ()];
	double * c = new double[mesh->get_dimentionality ()];
	// go through phase front nodes
	for (int i = 0; i < (int) phase_front_nodes.size (); i++)
	{
		// for each node get velocity
		task->get_node_velocity (phase_front_nodes[i], v);
		// calculate new position
		mesh->get_node_coordinates (phase_front_nodes[i], c);
		for (int k = 0, k_end = mesh->get_dimentionality (); k < k_end; k++)
			c[i] += parameters[0] * v[i];
		// move the node
		mesh->move_node (phase_front_nodes[i], c);
	}
	delete[] v;
	delete[] c;
}

Dynamic_Mesh_Builder::Dynamic_Mesh_Builder ()
{
}

Dynamic_Mesh_Builder::Dynamic_Mesh_Builder (const Dynamic_Mesh_Builder & dm)
{
}

Dynamic_Mesh_Builder::~Dynamic_Mesh_Builder ()
{
}

void Dynamic_Mesh_Builder::Adapt_Mesh (Task_pointer * task, Mesh_Prototype * new_mesh, Mesh_Prototype * cur_mesh, double time_layer)
{
	// copy current mesh
	new_mesh->copy (*cur_mesh);

	// change phase front
	adapt_phase_front (task, new_mesh);

	// move nodes
	
	// refresh mesh
	new_mesh->refresh_mesh ();
	// refold
}

void Dynamic_Mesh_Builder::set_phase_front (std::vector<int> nodes)
{
	phase_front_nodes.insert (phase_front_nodes.begin (), nodes.begin (), nodes.end ());
}

std::vector<int> Dynamic_Mesh_Builder::get_phase_front ()
{
	return phase_front_nodes;
}

void Dynamic_Mesh_Builder::set_monitor_function_parameter (int k, double value)
{
	if ((int)parameters.size() <= k)
		parameters.resize (k + 1);
	parameters[k] = value;
}
