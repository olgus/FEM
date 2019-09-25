#pragma once
#include "Prototype_Task.h"
#include "Mixed_Triangular_Mesh.h"

class Task_Mixed_Tr_Mesh : public Task <Mixed_Triangular_Mesh>
{
private:
	double function_f (int k_system, double * coordinates, int area) override; // returns value of function f in k-node
	double function_FCondition (int k_system, double * coordinates, int area, int boundary) override; // returns value for first condition in k-node
	double function_Theta (int k_system, double * coordinates, int area, int boundary) override; // returns value for first condition in coordinates
public:
	Task_Mixed_Tr_Mesh () {};
	Task_Mixed_Tr_Mesh (const Task_Mixed_Tr_Mesh & tmtm) {};
	~Task_Mixed_Tr_Mesh () {};
};