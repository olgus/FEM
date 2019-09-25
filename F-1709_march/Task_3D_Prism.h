#pragma once

#include "Task_3D_simple.h"

class Task_3D_Prism : public Task_3D <Prismatic_Mesh>
{
private:
	double function_FCondition (int k_system, double * coordinates, int area, int boundary) override; // returns value for first condition in k-node
	double function_f (int k_system, double * coordinates, int area) override; // returns value of function f in k-node
public:
	Task_3D_Prism () {};
	Task_3D_Prism (const Task_3D_Prism & ctt) {};
	~Task_3D_Prism () {};
};
