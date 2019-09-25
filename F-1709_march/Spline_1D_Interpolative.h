#pragma once
#include "Spline_1D_Pointer.h"
#include "Test_1D_Task.h"

class Spline_1D_Interpolative : public Spline_Pointer, public Task_Test_1D <Mesh_1D_Hermitian>
{
private:
	std::vector<double> function_values;

	void renumerate_functions ();
protected:
	
	virtual void build_system (int k_system) override; // set of functions for the equation

	// boundary conditions section
	void apply_boundary_conditions (int k_system) override;
	virtual void apply_first_boundary_conditions (int k_system) override; // applies first boundary conditions

public:
	void prepare (char * file_name);
	
	void prepare (std::vector <std::pair<double, double>> points_data, double * extra) override;
	bool get_solution_in_point (double * c, double * value) override;
	void solve_task () override;
	void get_boundaries (double * c0, double * cN) override;

	Spline_1D_Interpolative ();
	Spline_1D_Interpolative (const Spline_1D_Interpolative & spline);
	~Spline_1D_Interpolative ();
};