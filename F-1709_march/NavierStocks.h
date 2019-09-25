#pragma once
#include "Time_Dependent_Task.h"

#define NS_DISCRPANSCY 1e-7

class NavierStocks : public Time_Dependent_Task <Triangular_Mesh>
{
private:
	// physical parameters
	double phys_T0; // starting temperature
	double phys_T1; // heating temperature
	double PR; // Prandtl number
	double GR; // Grashof number
	double L_source; // length of the source (from 0 to the right)
	double P_step; // step for boundary conditions in W
protected:
	void set_n_systems () override; // sets the amount of equations of the task
	void build_system (int k_system) override; // build system
	void apply_first_boundary_conditions (int k_system) override; // applies first boundary conditions
	void print_solutions () override;
	
	double function_FCondition (int k_system, double * coordinates, int area, int boundary) override; // returns value for first condition in coordinates
	double function_starting_condition (int k_system, double * coordinates, int area) override;

	void solution_step (int k_system, int method, int d_type, int depth) override; // solution step non-linear, inplicit simple iteration
public:
	NavierStocks ();
	NavierStocks (const NavierStocks & nst);
	~NavierStocks ();
};