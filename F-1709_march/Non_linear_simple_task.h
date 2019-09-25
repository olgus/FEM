#pragma once
#include "Prototype_Task.h"
#define NLST_NL_PRECIS 1e-13

class Non_linear_simple_task : public Task <Triangular_Mesh>
{
private:
	double get_function_value (int k_system, int k_function, double * coordinates, int area) override;

	double function_f (int k_system, double * coordinates, int area) override;
	double function_sigma (int k_system, double * coordinates, int area);
	double function_FCondition (int k_system, double * coordinates, int area, int boundary) override;

	void build_system (int k_system) override; // set of functions for the equation

	void print_solutions (int knliter);
public:
	Non_linear_simple_task () {};
	Non_linear_simple_task (const Non_linear_simple_task & ctt) {};
	~Non_linear_simple_task () {};

	// its own solve function
	virtual bool solve_task (int solver_param[][5]) override;
	// and relax
	virtual double relax (int k_system, MathVector * prev_solution) override;
};