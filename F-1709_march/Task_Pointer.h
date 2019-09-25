#pragma once
#include "Mesh_Pointer.h"
#include "Sparse_Matrixes.h"

class Task_pointer
{
private:
protected:
public:
	Task_pointer () {};
	Task_pointer (const Task_pointer & task_pointer) {};
	~Task_pointer () {};

	compressed_matrix * eq_matrixes; // array of system matrixes
	Mesh_Prototype * mesh_pointer; // pointer to mesh *sigh*

	virtual double get_Q (int k_system, int time_layer, int k_func) { return 0.0; };

	virtual int get_function_global_number (int k_element, int k_func) { return 0; }; // returns global number of k_func local function on k_element
	virtual bool get_solution_in_point (int k_system, double * coordinates, double * value) { return false; }; // returns k_system solution on current time layer in point defined by coordinates
	virtual bool get_solution_in_point (int k_system, int k_time_layer, double * coordinates, double * value) { return false; }; // returns k_system solution on current time layer in point defined by coordinates
	virtual void get_mesh_boundaries (double * c0, double * cN) {}; // returns mesh boundaries for any kind of mesh
	virtual int get_mesh_dim () { return 0; }; // returns mesh dimentionality
	virtual Mesh_Prototype * get_mesh_pointer () { return NULL; };
	virtual double get_derivative_first (int k_system, int k_var, double * coordinates) { return 0.0; }; // returns value of first derivative by k_var
	virtual bool get_derivative (int k_system, int k_var, double * coordinates, double *value) { return 0.0; }; // returns value of first derivative by k_var
	virtual void get_full_solution (int k_system, MathVector * v) {};
	virtual void calc_gradient (int k_system) {};
	virtual bool get_isoline_section (int k_system, int var, int k_element, double value, double * c1, double * c2) { return false; };
	virtual int is_symmetrical () { return 0; };
	virtual void get_min_max (int k_system, double * min, double * max) {};

	// dynamic mesh section
	void get_node_velocity (int k_node, double * v) {};
};
