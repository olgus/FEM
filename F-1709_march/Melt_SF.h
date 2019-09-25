#pragma once
#include "Time_Dependent_Task.h"

#include "Painter.h" 

#define MAX_NONLINEAR_ITER_MELT 100
#define NONLINEAR_DISCR_MELT 1e-12
#define DISCR_CHANGE_MELT 1e-13
#define ZERO_MELT 1e-25
#define SOURCE_PRECIS 1e-6
#define INTEGRAL_POINTS 1000

// TODO: build_system, apply_first_boundary_conditions

template <class TypeMesh>
class Task_Melt_SF : public Time_Dependent_Task <TypeMesh>
{
private: 
	// physical parameters
	double PR; // Prandtl number
	double GR; // Grashof number
	double phase_change_og_position; // y-coordinate of the boundary
	double L_source; // length of the source (from 0 to the right)
	double phys_T1; // heating temperature
	double phys_T_phase_change; // phase change temperature
	double del_T_smoothing; // interval of T to smooth Dirac's fuction
	double Ste; // Stephan number
	double phys_T0_F; // fluid starting temperature
	double phys_T0_S; // solid starting temperature
	double source_power; // power of the heat source
	double lambda_fluid; // thermal conductivity of the fluid
	double lambda_solid; // thermal conductivity of the solidus
	double TEMP_T0; // minimum temp in the system
	bool density_inversion; // for water density inversion

	double turn_on_source_time;
	double bottom_temp;

	bool flag_append;
	bool flag_calc_start;

	int melt_rate_counter;
	
	double SCALE_TEMP, SCALE_TIME, SCALE_SIZE, SCALE_SPEED;
protected:
	virtual void set_n_systems () override; // sets the amount of equations of the task
	double function_starting_condition (int k_system, double * coordinates, int area) override;
	virtual void set_starting_conditions () override; // sets starting conditions for dynamic tasks

	// second boundary conditions SHOULD be applied normally, TEST that
	virtual void apply_first_boundary_conditions (int k_system) override; // applies boundary conditions
	virtual double function_FCondition (int k_system, double * coordinates, int area, int boundary) override; // returns value for first condition in coordinates
	virtual double function_f (int k_system, double * coordinates, int area) override; // returns value of function f in coordinates
	double get_lambda (double T);
	virtual void build_system (int k_system) override; // build system
	virtual void print_solutions () override;

	double Dirac_function (double T, double * coordinates); // smoothed Dirac function
	double density_inversion_coef (double T1); // returns coeffitient to multiply dT/dx in w equation
	void next_time_layer (int k_system, MathVector * solution); // moves past layers solutions up and puts a new one

	void print_extra_data () override;

	double get_melt_volume ();
	void get_Qs (double y, double * Qt, double * Qc);
public:
	//bool solve_task (int method, int d_type, int depth) override;
	//bool solve_task (int solver_param[][5]);
	void get_scales ();
	void get_scales (char * scales_file_name);
	void use_density_inversion ();
	void set_starting_conditions_from_task (Task_pointer * task);

	void reset_parameters (double _SCALE_TEMP, double _GR, double _Ste, double _phys_T_phase_change, double _phys_T0_F);
	double get_front_position (double x);

	Task_Melt_SF ();
	Task_Melt_SF (char * parameters_file_name);
	Task_Melt_SF (const Task_Melt_SF & nst);
	~Task_Melt_SF ();
};