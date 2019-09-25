#pragma once
#include "Prototype_Task.h"

class Polar_Task : public Task <Triangular_Polar_Mesh>
{
private:
	double function_f (int k_system, double * coordinates, int area) override; // returns value of function f in k-node
	double function_FCondition (int k_system, double * coordinates, int area, int boundary) override; // returns value for first condition in k-node
	double function_Theta (int k_system, double * coordinates, int area, int boundary) override; // returns value for first condition in coordinates

public:
	Polar_Task () {};
	Polar_Task (const Polar_Task & ctt) {};
	~Polar_Task () {};
};
