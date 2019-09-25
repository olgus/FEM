#pragma once
#include "Prototype_Task.h"
#include "Point_source_2D.h"
#include "Task_Additional_Field.h"

class Task_2D_source_w_field_selection: public Task <Triangular_Mesh>
{
private:
	int n_sources; // amount of point sources
	point_source_2D * sources; // sources data
	int n_insertions; // amount of insertions
	Task_2D_Additional_field * tasks_additional_field;
	virtual double function_f (int k_system, double * coordinates, int area) override; // returns value of function f in coordinates

	virtual void apply_boundary_conditions (int k_system) override;
	double get_main_field (double * coordinates); // u0
public:
	Task_2D_source_w_field_selection ();
	~Task_2D_source_w_field_selection ();

	void prepare (int N_insertions, char * file_name_sources, char * file_name_lambdas, char * file_name_insertion_nodes[], char * file_name_insertion_triangles[]); // prepare: read sources, read insertions, read lambdas, get meshes for all insertions
	void solve_task (); // solve
	
	virtual bool get_solution_in_point (int k_system, double * coordinates, double * value) override; // get solution
	virtual int get_mesh_dim () override; // returns mesh dimentionality
};