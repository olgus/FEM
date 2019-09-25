#pragma once
#include "Time_Dependent_Task.h"
#define NLNSTT_NL_PRECIS 1e-13

class Nonlinear_nonstationary_task : public Time_Dependent_Task <Triangular_Mesh>
{
private:
	// base matrix for the parts that don't change
	compressed_matrix * BM; 
	void build_BM (int k_system);

	double get_function_value (int k_system, int k_function, double * coordinates, int area) override;

	double function_f (int k_system, double * coordinates, int area) override;
	double function_starting_condition (int k_system, double * coordinates, int area) override;
	double function_sigma (int k_system, double * coordinates, int area);
	double function_FCondition (int k_system, double * coordinates, int area, int boundary) override;
	using Task::get_function_coefficients;
	void get_function_coefficients (int k_system, int function, MathVector * coef);
	void build_system (int k_system) override; // set of functions for the equation

	void print_solutions (int knliter);
	void print_solutions ();
public:
	Nonlinear_nonstationary_task () {};
	Nonlinear_nonstationary_task (const Nonlinear_nonstationary_task & ctt) {};
	~Nonlinear_nonstationary_task () {};

	// its own solve function
	virtual bool solve_task (int solver_param[][5]) override;
	// and relax
	virtual double relax (int k_system, double w, const MathVector & prev_non_linear_layer_solution);
	virtual double relax (int k_system, double w_left, double w_right, const MathVector & prev_non_linear_layer_solution);
	virtual double relax (int k_system, const MathVector & prev_non_linear_layer_solution);

	void sigma_test ();
};