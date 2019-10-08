#pragma once
#include "Time_Dependent_Task.h"
#include "Spline_1D_Smooth.h"

#define TCONVEC_PRECIS 1e-12
#define TCONVEC_DISCR_PRECIS 3e-12
#define BOUND_FUNC_HEAT_SOURCE 1001
#define SOLID_PRECIS 1e-7
#define MAX_TC_NONLINEAR_ITER 1000

class Task_convec : public Time_Dependent_Task <Triangular_Mesh>
{
private:
	double PR, GR;
	int PRGRconst; // 1 if no spline needed, 0 if otherwise
	double T_min, THETA; // minimum temp in the system, temp scale
	double T0_solid, T0_fluid, T_source;
	int LOCALIZED_BOTTOM_HEAT_SOURCE;
	double heat_source_length;
	double lambda_fluid;
	double lambda_solid;
	int MELT_SETUP;
	double Ste;
	double T_phase_change; // dimentionless temperature
	int USE_DENSITY_FUNCTION; // 0 for normal linear dependence
	int density_function; // density function number

	virtual void get_conditions (int k_system, char * file_name) override; // gets boundary conditions

	double get_function_value (int k_system, int k_function, double * coordinates, int area) override;

	double function_f (int k_system, double * coordinates, int area) override;
	double function_starting_condition (int k_system, double * coordinates, int area) override;
	double function_FCondition (int k_system, double * coordinates, int area, int boundary) override;
	void get_function_coefficients (int k_system, int function, MathVector * coef);
	void get_density_coefficients (int k_element, int function, MathVector * coef);
	double get_density_T_IP_scaled ();

	void build_system (int k_system) override; // set of functions for the equation

	void print_solutions (int knliter);
	void print_solutions ();
	void get_Theta_edge (int k_system, int k_element, int n1, int n2, int * functions, int boundary, MathVector * Theta) override; 
	virtual double Lambda (int k_element); // returns physical lambda value

	virtual void apply_boundary_conditions (int k_system) override; // applies boundary conditions
	virtual void apply_first_boundary_conditions (int k_system) override; // applies first boundary conditions
	virtual void apply_second_boundary_conditions (int k_system) override; // applies first boundary conditions

	// for faster first boundary conditions
	int K_NON_LINEAR_ITER;
	int ** ffb; // flags for the first boundary condition, 0 - no first condition, anything else is the number of the function
	void set_ffb_flags ();
	int * state; // 0 for fluid, 1 for solid, 2 for mush
	void reset_states ();

	// PR GR as splines
	Spline_1D_Smooth * spline_Pr;
	Spline_1D_Smooth * spline_Gr;
	double get_Pr (int k_element);
	double get_Gr (int k_element);
protected:
	virtual void set_n_systems () override; // sets the amount of equations of the task
public:
	Task_convec ();
	Task_convec (char * file_param);
	Task_convec (const Task_convec & ctt) {};
	~Task_convec ();

	virtual void get_min_max (int k_system, double * min, double * max) override;
	// its own solve function
	virtual bool solve_task (int solver_param[][5]) override;
	// and relax
	virtual double relax (int k_system, double w, const MathVector & prev_non_linear_layer_solution);
	virtual double relax (int k_system, double w_left, double w_right, const MathVector & prev_non_linear_layer_solution);
	virtual double relax (int k_system, const MathVector & prev_non_linear_layer_solution);

	void save_front_pos (char * file_name);
};