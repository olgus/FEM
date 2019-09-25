#pragma once
#include "Time_Dependent_Task.h"

#define ZERO_Plumes 1e-10
#define PL_DISCRPANSCY 9.9e-4

template <class TypeMesh>
class Plumes_Prototype : public Time_Dependent_Task <TypeMesh>
{
private:
	// physical parameters
	double phys_T0; // starting temperature
	double phys_T1; // heating temperature
	double PR; // Prandtl number
	double GR; // Grashof number
	double L_source; // length of the source (from 0 to the right)
	double P_step; // step for boundary conditions in W
	int slip; // defines boundary conditions for w

	double therm_solid; // starting temperature
	double therm_fluid; // heating temperature
protected:
	virtual void set_n_systems () override; // sets the amount of equations of the task
	virtual void build_system (int k_system) override; // build system
	virtual void apply_first_boundary_conditions (int k_system) override; // applies first boundary conditions
	virtual void print_solutions () override;
	virtual double function_FCondition (int k_system, double * coordinates, int area, int boundary) override; // returns value for first condition in coordinates
	double function_starting_condition (int k_system, double * coordinates, int area) override;

	virtual void solution_step (int k_system, int method, int d_type, int depth) override; // solution step non-linear, inplicit simple iteration
public:
	Plumes_Prototype ();
	Plumes_Prototype (const Plumes_Prototype & nst);
	~Plumes_Prototype ();
};
