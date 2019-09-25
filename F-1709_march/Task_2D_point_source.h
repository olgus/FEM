#pragma once
#include "Prototype_Task.h"
#include "Point_source_2D.h"

class Task_2D_point_source : public Task <Triangular_Mesh>
{
private:
	int n_sources; // amount of point sources
	point_source_2D * sources; // sources data

	double function_FCondition (int k_system, double * coordinates, int area, int boundary) override; // returns value for first condition in k-node
	virtual double function_f (int k_system, double * coordinates, int area) override; // returns value of function f in coordinates
	virtual double function_Theta (int k_system, double * coordinates, int area, int boundary) override; // returns value for first condition in coordinates

	//virtual void get_local_B (int k_system, int k_element, MathVector * B);
	virtual void apply_boundary_conditions (int k_system) override;
public:
	Task_2D_point_source ();
	~Task_2D_point_source ();

	virtual void prepare (char * mesh_file_name, char * bound_file_name, bool save) override;
	virtual void prepare (char * file_name_nodes, char * file_name_elements, char * bound_file_name) override;
	virtual void prepare (char * file_name_nodes, char * file_name_elements, char * bound_file_name, char * file_sources) ;
};