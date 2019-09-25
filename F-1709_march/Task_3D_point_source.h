#pragma once
#include "Task_3D_Prism.h"

struct point_source
{
	double x, y, z;
	double power;
};

class Task_3D_point_source : public Task_3D_Prism
{
private:
	int n_sources; // amount of point sources
	point_source * sources; // sources data

	double function_FCondition (int k_system, double * coordinates, int area, int boundary) override; // returns value for first condition in k-node

	virtual void get_local_B (int k_system, int k_element, MathVector * B) override;
	virtual void apply_boundary_conditions (int k_system) override;
	virtual double Gamma (int k_system, int k_function, int area) override;
public:
	Task_3D_point_source ();
	~Task_3D_point_source ();

	virtual void prepare (char * mesh_file_name, char * bound_file_name, bool save) override;
	virtual void prepare (char * file_name_nodes, char * file_name_elements, char * bound_file_name) override;
};