#pragma once
#include "Prototype_Task.h"

class Cartesian_Test_Task : public Task <Triangular_Mesh>
{
private:
	double function_f (int k_system, double * coordinates, int area) override; 
	double function_FCondition (int k_system, double * coordinates, int area, int boundary) override; 
	double function_Theta (int k_system, double * coordinates, int area, int boundary) override; // returns value for first condition in coordinates
																								 // systems' parameters
	double Lambda (int k_system, int k_function, int area) override; // returns lambda value, 1.0 by default
	double Gamma (int k_system, int k_function, int area) override; // returns gamma value, 1.0 by default

public:
	Cartesian_Test_Task () {};
	Cartesian_Test_Task (const Cartesian_Test_Task & ctt) {};
	~Cartesian_Test_Task () {};
};