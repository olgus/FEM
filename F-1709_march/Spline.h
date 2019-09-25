#pragma once
#include "Prototype_Task.h"

class Spline : public Task <Rectangular_Mesh_Spline>
{
private:
	int n_sol_points; // amount of points that spline smoothes
	double alpha, beta; // parameters of the spline
	double def_value; // default value
	MathVector * solution_values; // f~, values in solution_points
	Matrix * solution_points; // x~, coordinates of solutions points
	MathVector * weigths; // spline weigths

	void get_local_A_matrix (int k_element, Matrix * A_matrix); // returns A matrix for k_element
	void get_local_B (int k_element, MathVector * B_vector); // sets B vector for k_element
	double function_f (int k_system, double * coordinates, int area) override; // returns value of function f in k-node

	void build_system (int k_system) override; // set of functions for the equation
	void apply_boundary_conditions (int k_system) override; //applies boundary conditions	
	void print_solutions () override; // prints current time layer's solutions into files
	void get_conditions (int k_system, char * file_name) override; // gets boundary conditions
public:
	Spline ();
	Spline (const Spline & spline); 
	~Spline ();

	bool get_solution_in_point (int k_system, double * coordinates, double * value) override; // returns k_system solution on current time layer in point defined by coordinates
	bool get_first_derivative (int k_system, int k_var, double * coordinates, double * value); // returns k_var derivative of spline on current time layer in point defined by coordinates
	void prepare (char * mesh_file_name, char * bound_file_name, bool save) override; // call parent one, add point input from files
	void prepare (Mesh_Prototype * parent_mesh, char * spline_for_mesh_file); // point input from files and uniform mesh
	void prepare (double * c0, double * cN, double lvl, char * spline_for_mesh_file); // point input from files and uniform mesh
	void save_slice_first_derivative (int k_system, char * file_sour, char * file_dest, bool keep_points) override; // saves slice of the solution 
	void save_slice_mixed_derivative (int k_system, char * file_sour, char * file_dest, bool keep_points); // saves slice of the solution 
	void default_value (double value);
	void set_alpha_beta (double Alpha, double Beta);
};