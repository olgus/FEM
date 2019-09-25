#pragma once
#include "Prototype_Task.h"

class Task_Vector_2D : public Task <Triangle_Vector_Mesh>
{
private:
protected:
	// build local B 
	void get_local_B (int k_system, int k_element, MathVector * B) override;
	void get_FCondition_edge (int k_system, int k_element, int n1, int n2, int * functions, int boundary, MathVector * FCondition) override; // get vector q for first condition
	// second
	void get_Theta_edge (int k_system, int k_element, int n1, int n2, int * functions, int boundary, MathVector * Theta) override; // get vector of theta values

	// task functions
	void function_f (int k_system, double * coordinates, int area, double * f_value);
	void function_FCondition (int k_system, double * coordinates, int area, int boundary, double * value); 
	// systems' parameters
	double Lambda (int k_system, int k_function, int area) override; // returns lambda value, 1.0 by default
	double Gamma (int k_system, int k_function, int area) override; // returns gamma value, 1.0 by default
	double function_Theta (int k_system, double * coordinates, int area, int boundary) override; // returns value for first condition in coordinates


public:
	Task_Vector_2D () {};
	Task_Vector_2D (const Task_Vector_2D & task) {};
	~Task_Vector_2D () {};

	bool get_solution_in_point (int k_system, double * coordinates, double * value) override; 
	void save_evenly_2D (int k_system, double * h, char * file_dest, bool keep_points) override; // save values in points that are calculated from c0 to cN with h step

};