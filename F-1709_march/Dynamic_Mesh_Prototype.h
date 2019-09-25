#pragma once
#include <stdio.h>
#include <vector>
#include <memory>

#include "Mesh_Pointer.h"
#include "Triangular_Mesh.h"
#include "Task_Pointer.h"

// dynamic mesh
// inherits from no one
// has:
// triangular mesh as computational space 
// monitor function
// works like:
// gets a new time layer value
// builds matrix for the discretized task
// solves task for new tlv
// takes a mesh_pointer to an existing mesh
// writes new mesh into that mesh_pointer
  
class Dynamic_Mesh_Builder
{
private:
	// matrix for the task
	void build_task_matrix (); // builds matrix for the discretized task
	void fix_boundary_nodes (); // fix boundary nodes
protected:
	double monitor_function (double * coordinates);

	// function to adjust the phase front
	void adapt_phase_front (Task_pointer * task, Mesh_Prototype * new_mesh); // moves nodes that are the closest to the current phase front
	std::vector<int> phase_front_nodes; // keep numbers of nodes that are phase front
	std::vector<double> parameters; // parameters for different functions
	// parameters[0] = del_t
public: 
	Dynamic_Mesh_Builder ();
	Dynamic_Mesh_Builder (const Dynamic_Mesh_Builder & dm);
	~Dynamic_Mesh_Builder ();

	void Adapt_Mesh (Task_pointer * task, Mesh_Prototype * new_mesh, Mesh_Prototype * cur_mesh, double time_layer);
	void set_phase_front (std::vector<int> nodes);
	std::vector<int> get_phase_front ();

	void set_monitor_function_parameter (int k, double value);
};