#pragma once
#include <math.h>
#include <stdio.h>
#include <conio.h>
#include <string.h>
#include <vector>
#include <map>

#include "Task_Pointer.h"
#include "Mesh_Types.h"
#include "Painter.h"

#define ZERO_task 1e-10
#define IMPLICIT_SIMPLE_ITERATION_MAX_ITER 100
#define DISCR_CHANGE 1e-13
#define NONLINEAR_DISCR 1e-12
#define MAX_NONLINEAR_ITER 100

template <class Type_Mesh>
class Task : public Task_pointer
{
private:
	bool SR; // flag of derivative by SR method
protected:
	// setting section
	virtual void get_conditions (int k_system, char * file_name); // gets boundary conditions

	Painter * painter;

	int symmetrical; // flag of vertical symmetry, 0 - none, 1 - absolute, 2 - with sign change
	bool non_linear; // flag of non-linearity
	// relaxation
	virtual double relax (int k_system, MathVector * prev_solution);

	// actual solution section
	virtual void solution_step (int k_system, int method, int d_type, int depth); // solution step. Might be non-linear
	virtual void solution_step (int k_system, int * solver_param); // solution step. Might be non-linear

	// time section
	int n_time_layers; // 2 for stationary task
	int current_time_layer; // number of current time layer
	double * time_layers; // values of time layers
	int time_sampling; // amount of previous time_layers needed
	virtual void get_time_coefficients (MathVector * c); // get scheme coefficients on current time layers
	virtual void get_time_approx (double t, MathVector * c); // get scheme coefficients on current time layers

	// time_stamps
	bool use_time_mapping;
	std::vector <double> time_stamps;

	// extra variables section
	FILE * log_file;
	bool prepared;

	// mesh section
	Type_Mesh mesh; // a mesh, strictly one
	
	// conditions
	Condition ** conditions;

	// setting section
	virtual void set_time_layers (char * file_name); // sets time layers
	virtual void set_n_systems (); // sets the amount of equations of the task
	virtual bool build_portrait (int k_system); // build matrix portrait
	virtual bool build_portrait (int k_system, Type_Mesh * mesh); // build matrix portrait

	// systems and solutions section
	int n_systems; // amount of differential equations that describe the task
	MathVector ** previous_time_layers_solutions; // values from previous k-layers
	MathVector * non_linear_layers_solutions; // values from previous k-layers
	MathVector * derivative; // saves derivates there, if needed
	virtual void next_time_layer (int k_system); // moves past layers solutions up and puts a new one
	virtual void next_time_layer (int k_system, MathVector * solution); // moves past layers solutions up and puts a new one

	// functions section
	int N_functions; // amount of functions
	virtual int get_function_global_number (int k_element, int k_func) override; // returns global number of k_func local function on k_element
	virtual void get_element_functions (int k_element, int * functions); // returns numbers of functions of the element
	
	std::map<int, double> lambda_val; // lambdas
	std::map<int, double> gamma_val; // gammas

	// system functions section

	// default 0
	virtual void set_starting_conditions (); // sets starting conditions for dynamic tasks
	// default (G + M) * x = B
	virtual void build_system (int k_system); // set of functions for the equation

	// boundary conditions section
	virtual void apply_boundary_conditions (int k_system); // applies boundary conditions
	virtual void apply_first_boundary_conditions (int k_system); // applies first boundary conditions
	virtual void apply_second_boundary_conditions (int k_system); // applies second boundary conditions
	// second
	virtual void get_Theta_edge (int k_system, int k_element, int n1, int n2, int * functions, int boundary, MathVector * Theta); // get vector of theta values
	virtual double function_Theta (int k_system, double * coordinates, int area, int boundary); // returns value for first condition in coordinates
	// first
	virtual void get_FCondition_edge (int k_system, int k_element, int n1, int n2, int * functions, int boundary, MathVector * FCondition); // get vector of q for first condition
	virtual double function_FCondition (int k_system, double * coordinates, int area, int boundary); // returns value for first condition in coordinates

	// derivative section
	// first derivative by using basis functions derivatives
	virtual double get_derivative_first (int k_system, int k_var, double * coordinates) override; // returns value of first derivative by k_var
	// second derivative that also uses numerical derivative, because basis functions for triangles do not have second derivative
	virtual double get_derivative_second (int k_system, int k_var1, int k_var2, double h1, double h2, double * coordinates); // returns value of second derivative  by k_var1 k_var2
	virtual double get_derivative_SR (int k_system, int k_var, double * coordinates); // returns value of first derivative by k_var

	// systems' parameters
	virtual double Lambda (int k_system, int k_function, int area); // returns lambda value, 1.0 by default
	virtual double Gamma (int k_system, int k_function, int area); // returns gamma value, 1.0 by default

	// f functions' section
	virtual void get_local_B (int k_system, int k_element, MathVector * B);
	virtual double function_f (int k_system, double * coordinates, int area); // returns value of function f in coordinates
	
	virtual double get_function_value (int k_system, int k_function, double * coordinates, int area);
	virtual void get_function_coefficients (int k_system, int k_element, int function, MathVector * coef);
public:
	Task (); // constructor 
	Task (const Task & task); // constructor-copy
	virtual ~Task (); // destructor

	void painter_pointer (Painter * p);

	// prepare/solve section
	virtual void prepare (char * mesh_file_name, char * bound_file_name, bool save); // prepare function (sets mesh, sets starting conditions, matrix portraits and stuff like that)
	virtual void prepare (char * mesh_file_name, char * bound_file_name, char * time_file_name, bool save); // prepare function (sets mesh, sets starting conditions, matrix portraits and stuff like that)
	virtual void prepare (char * mesh_file_name, char * bound_file_name[], char * time_file_name, bool save); // prepare function (sets mesh, sets starting conditions, matrix portraits and stuff like that)

	virtual void prepare (char * nodes_file_name, char * elements_file_name, char * bound_file_name); // prepare function (sets mesh, sets starting conditions, matrix portraits and stuff like that)
	virtual void prepare (char * nodes_file_name, char * elements_file_name, char * bound_file_name, char * time_file_name); // prepare function (sets mesh, sets starting conditions, matrix portraits and stuff like that)
	virtual void prepare (char * nodes_file_name, char * elements_file_name, char * bound_file_name[], char * time_file_name); // prepare function (sets mesh, sets starting conditions, matrix portraits and stuff like that)

	virtual void read_materials (char * materials_file_name);
	// output section
	virtual void print_solutions (); // prints current time layer's solutions into files	
	virtual void print_extra_data (); 

	virtual bool solve_task (int method, int d_type, int depth); // solve function. Universal?
	virtual bool solve_task (int * solver_param); // solve function. Universal?
	virtual bool solve_task (int solver_param[][5]);
	virtual void read_solution (char * file_name, int k_system, int k_time_layer);

	virtual void set_symmetry (int Symmetrical);

	// value section
	virtual bool get_solution_in_point (int k_system, double * coordinates, double * value) override; // returns k_system solution on current time layer in point defined by coordinates
	virtual bool get_solution_in_point (int k_system, int k_time_layer, double * coordinates, double * value) override; // returns k_system solution on current time layer in point defined by coordinates
	virtual bool get_non_linear_solution_in_point (int k_system, double * coordinates, double * value); // returns k_system solution on current time layer in point defined by coordinates

	// solution saving section
	virtual void save_solution_in_points (int k_system, char * file_sour, char * file_dest, bool keep_points); // saves k_system solution on current time layer in points from file into another file
	virtual void save_solution_in_points_first_derivative (int k_system, char * file_sour, char * file_dest, bool keep_points, int k_var, double h_der); // saves k_system solution on current time layer in points from file into another file
	virtual void save_solution_in_points_second_derivative (int k_system, char * file_sour, char * file_dest, bool keep_points, int k_var1, int k_var2, double h1, double h2); // saves k_system solution on current time layer in points from file into another file
	virtual void save_slice (int k_system, char * file_sour, char * file_dest, bool keep_points); // saves slice of the solution 
	virtual void save_slice (int k_system, double * boundaries, char * file_sour, char * file_dest, bool keep_points); // saves slice of the solution 
	virtual void save_slice (int k_system, double seek_value,  int k_var, char * file_dest, bool keep_points); // saves slice of the solution 
	virtual void save_slice_weighted (int k_system, double * weights, char * file_sour, char * file_dest, bool keep_points); // saves slice of the solution 
	virtual void save_slice_first_derivative (int k_system, char * file_sour, char * file_dest, bool keep_points); // saves slice of the solution 
	virtual void save_slice_first_derivative (int k_system, int der_var, char * file_sour, char * file_dest, bool keep_points); // saves slice of the solution 
	virtual void save_slice_first_derivative_weighted (int k_system, int der_var, double * weights, char * file_sour, char * file_dest, bool keep_points); // saves slice of the solution 
	virtual void save_slice_first_derivative_absolute (int k_system, int der_var, double * weights, char * file_sour, char * file_dest, bool keep_points); // saves slice of the solution 
	virtual void save_slice_second_derivative (int k_system, char * file_sour, char * file_dest, bool keep_points, int k_var2, double h1, double h2); // saves slice of the solution 
	virtual void save_evenly_2D (int k_system, double * h, char * file_dest, bool keep_points); // save values in points that are calculated from c0 to cN with h step
	virtual void save_evenly_2D_first_derivative (int k_system, int k_var, double * h, char * file_dest, bool keep_points); // save values in points that are calculated from c0 to cN with h step
	virtual void save_evenly_2D_first_derivative_SR (int k_system, int k_var, double * h, char * file_dest, bool keep_points); // save values in points that are calculated from c0 to cN with h step
	
	virtual void fprint_edges_by_nodes (char * file_name);
	virtual void fprint_edges (char * file_name);
	virtual void fprint_functions_by_elements (char * file_name);
	virtual void fprint_edges_by_elements (char * file_name);
	virtual void fprint_elements_by_edges (char * file_name);

	// other useful functions
	virtual void get_mesh_boundaries (double * c0, double * cN) override; // returns mesh boundaries for any kind of mesh
	virtual int get_mesh_dim () override; // returns mesh dimentionality
	Mesh_Prototype * get_mesh_pointer () override;
	virtual double get_Q (int k_system, int time_layer, int k_func) override;

	virtual void get_min_max (int k_system, double * min, double * max) override;

	virtual void read_solution (int k_system, char * file_name);
	virtual void get_full_solution (int k_system, MathVector * v) override;

	virtual int is_symmetrical () override;

	// lambda gamma section
	void get_lambda_values (char * file_name);
	void get_gamma_values (char * file_name);

	// derivative section
	virtual bool get_derivative (int k_system, int k_var, double * coordinates, double *value) override; // returns value of first derivative by k_var
	virtual void calc_derivative (int k_system, int k_var, double * c0, double * cN);
	virtual void calc_gradient (int k_system) override;
	
	// temporary
	void reset_N_functions (int n);

	// refinement section
	void remesh (std::vector <int> elements_to_refine);

	// for painter
	// var = 0..n - derivative, > n - gradient, < 0 - solution value
	virtual bool get_isoline_section (int k_system, int var, int k_element, double value, double * c1, double * c2) override;

	void set_solution (int k_system, int k_time_layer, char * file_name);
};

// functions to change:
// for new 2D task
// simple:
// function_f
// function_FCondition
// get_Theta_edge
// if functions are not 1 for a node:
// get_N_functions
// get_function_global_number
// possibly:
// get_local_B
// apply_first_boundary_conditions
// get_FCondition_edge