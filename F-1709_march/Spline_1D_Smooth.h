#pragma once
#include "Spline_1D_Pointer.h"
#include "Test_1D_Task.h"

#define ZERO_1D_spline_smooth 1e-5

class Spline_1D_Smooth : public Spline_Pointer, public Task_Test_1D <Mesh_1D_Hermitian>
{
private:
	int n_sol_points; // amount of points that spline smoothes
	double alpha, beta; // parameters of the spline
	MathVector * solution_values; // f~, values in solution_points
	MathVector * solution_points; // x~, coordinates of solutions points

	void get_local_A_matrix (int k_element, Matrix * A_matrix); // returns A matrix for k_element
	void get_local_B (int k_element, MathVector * B_vector); // sets B vector for k_element
	double function_f (int k_system, double * coordinates, int area) override; // returns value of function f in k-node

	void build_system (int k_system) override; // set of functions for the equation
	void apply_boundary_conditions (int k_system) override; //applies boundary conditions	

protected:
public:
	Spline_1D_Smooth ();
	Spline_1D_Smooth (const Spline_1D_Smooth & spline);
	~Spline_1D_Smooth ();

	void prepare (char * file_name, double * extra);
	void prepare (std::vector <std::pair<double, double>> points_data, double * extra) override;

	void prepare (char * file_name, double * extra, double n_points_sqrt);
	void prepare (std::vector <std::pair<double, double>> points_data, double * extra, double n_points_sqrt);

	bool get_solution_in_point (double * c, double * value) override;
	void solve_task () override;
	void get_boundaries (double * c0, double * cN) override;

	void prepare (double c0, double cN, int n_x, double Alpha, double Beta, std::vector <std::pair <double, double>> points_data); // prepare mesh and function values
	void prepare (double c0, double cN, int n_x, double Alpha, double Beta, char * file_name); // prepare mesh and function values
};