#pragma once
#include "Prototype_Task.h"

template <class Type_Mesh>
class Task_Test_1D : public Task <Type_Mesh>
{
private:
	virtual void apply_boundary_conditions (int k_system) override; // applies boundary conditions
protected:
	double function_FCondition (int k_system, double * coordinates, int area, int boundary) override; // returns value for first condition in k-node
	double function_f (int k_system, double * coordinates, int area) override;
	double function_FCondition_der (int k_system, double * coordinates, int area, int boundary);
public:
	Task_Test_1D () {};
	Task_Test_1D (const Task_Test_1D & ctt) {};
	virtual ~Task_Test_1D () {};

	void save_slice (int k_system, char * file_sour, char * file_dest, bool keep_points) override; // saves slice of the solution 
};