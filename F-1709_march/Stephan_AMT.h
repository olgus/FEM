#pragma once
#include "Time_Dependent_Task.h"

#include "Painter.h" 

#define MAX_NON_LINEAR_ITER_ST 100

// two-dimentional heat task
// supposedely rotated by 30 degrees
// no idea how that is supposed to work

namespace autom_St_task
{
	double Laplace_function (double z); // laplace function for z = x cos + y sin
	double get_ph_ch_position (double t, double beta);
	double get_solution (double t, double x, double beta);
}

template <class TypeMesh>
class Stephan_AMT : public Time_Dependent_Task <TypeMesh>
{
protected:
	Painter * painter;

	double angle; // angle of the phase change position, in dgr
	double beta; // beta coef
	double phase_change_og_position; // z-coordinate of the boundary
	double del_T_smoothing; // interval of T to smooth Dirac's fuction
	double phys_T_phase_change; // phase change temperature
	double T_heat, T_cool; // boundary temperatures
	double L; // phase change heat capacity

	virtual void set_n_systems () override; // sets the amount of equations of the task
	double function_starting_condition (int k_system, double * coordinates, int area) override;

	virtual double function_FCondition (int k_system, double * coordinates, int area, int boundary) override; // returns value for first condition in coordinates
	virtual double function_f (int k_system, double * coordinates, int area) override; // returns value of function f in coordinates
	virtual void print_solutions () override;

public:
	double get_ph_ch_x_position ();
	void painter_pointer (Painter * p);

	Stephan_AMT ();
	Stephan_AMT (char * parameters_file_name);
	~Stephan_AMT ();
};

// task that applies phase change boundary condition directly through smoothing over crossing elements
template <class TypeMesh>
class Stephan_AMT_smoothing : public Stephan_AMT<TypeMesh>
{
private: 
protected:
	virtual void build_system (int k_system) override; // build system
	double Dirac_function (double T, double * coordinates); // smoothed Dirac function
public:

	Stephan_AMT_smoothing ();
	Stephan_AMT_smoothing (char * parameters_file_name);
	Stephan_AMT_smoothing (const Stephan_AMT_smoothing & samts);
	~Stephan_AMT_smoothing ();
};

// task that uses coordinate transformation 
template <class TypeMesh>
class Stephan_AMT_coord_transf : public Stephan_AMT <TypeMesh>
{
};

// task that adapts the mesh by adding a node
//template <class TypeMesh>
class Stephan_AMT_mesh_adapt : public Stephan_AMT <Mesh_1D_L1>
{ 
private:
	Mesh_Prototype ** mesh_pointers_by_time_levels;// meshes for different time levels
	void adapt_mesh (); // adapt mesh
	void reset_areas (); // reset areas

	int L_heat_node;
	bool L_base_node;
public:
	// special solve function
	virtual bool solve_task (int solver_param[][5]) override;
	// special function_starting condition, i assume
	// don't forget to set L_heat_node
	virtual void set_starting_conditions () override;	

	Stephan_AMT_mesh_adapt ();
	~Stephan_AMT_mesh_adapt ();
};