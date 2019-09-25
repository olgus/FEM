#pragma once
#include "Prototype_Task.h"

class Task_2D_Additional_field : public Task <Triangular_Mesh>
{
private:
	int k_insertion;
	Task_pointer * task_main_field;
	virtual void get_local_B (int k_system, int k_element, MathVector * B) override;
	virtual double function_FCondition (int k_system, double * coordinates, int area, int boundary) override; // returns value for first condition in coordinates
public:
	Task_2D_Additional_field ();
	~Task_2D_Additional_field ();

	void prepare (int K_insertion, Task_pointer * Task_main_field, char * file_name_nodes, char * file_name_triangles, std::map <int, double> Lambda_val); // set uniform mesh and task
};