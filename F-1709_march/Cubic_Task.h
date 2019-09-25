#pragma once

#include "Task_3D_simple.h"

class Cubic_Task : public Task_3D <Cubic_Mesh>
{
private:
	double function_f (int k_system, double * coordinates, int area) override; // returns value of function f in k-node
	double function_FCondition (int k_system, double * coordinates, int area, int boundary) override; // returns value for first condition in coordinates
public:
	Cubic_Task () {};
	Cubic_Task (const Cubic_Task & ct) {};
	~Cubic_Task () {};
};