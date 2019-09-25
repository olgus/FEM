#pragma once
#include "Prototype_Task.h"

class Task_Triangular_Hierarchical : public Task <Triangular_Mesh_Hier>
{
private:
	double function_f (int k_system, double * coordinates, int area) override;
	double function_FCondition (int k_system, double * coordinates, int area, int boundary) override;
	double function_Theta (int k_system, double * coordinates, int area, int boundary) override; // returns value for first condition in coordinates

public:
	Task_Triangular_Hierarchical () {};
	Task_Triangular_Hierarchical (const Task_Triangular_Hierarchical & ctt) {};
	~Task_Triangular_Hierarchical () {};
}; 