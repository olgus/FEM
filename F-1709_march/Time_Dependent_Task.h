#pragma once
#include "Prototype_Task.h"

template <class Type_Mesh>
class Time_Dependent_Task : public Task <Type_Mesh>
{
private:
	// map time
	void map_time (char * file_name);
protected:
	bool start_cond;

	virtual double function_f (int k_system, double * coordinates, int area) override; // returns value of function f in coordinates
	virtual double function_FCondition (int k_system, double * coordinates, int area, int boundary) override; // returns value for first condition in coordinates
																							
	// starting conditions
	virtual void set_starting_conditions () override; // sets starting conditions for dynamic tasks
	virtual double function_starting_condition (int k_system, double * coordinates, int area);
	virtual void get_local_SCB (int k_system, int k_element, MathVector * B); // sets starting conditions vector on the element

	// time section
	void get_time_coefficients (MathVector * c) override; // get scheme coefficients on current time layers
	virtual void get_time_approx (double t, MathVector * c) override; // get scheme coefficients on current time layers
	virtual void set_time_layers (char * file_name) override;  // sets time layers

public:
	Time_Dependent_Task (); // constructor
	Time_Dependent_Task (const Time_Dependent_Task & tdt); // constructor-copy
	virtual ~Time_Dependent_Task (); // destructor

	void set_scheme_order (int Scheme_order); // sets scheme_order 1-3

	// parabolic task
	virtual void build_system (int k_system) override; // set of functions for the equation
	
	using Task<Type_Mesh>::prepare;
	virtual void prepare (char * mesh_file_name, char * bound_file_name, char * time_file_name, char * stamps_file_name, bool save);
	virtual void prepare (char * mesh_file_name, char * bound_file_name[], char * time_file_name, char * stamps_file_name, bool save);
	virtual void prepare (char * nodes_file_name, char * elements_file_name, char * bound_file_name, char * time_file_name, char * stamps_file_name);
	virtual void prepare (char * nodes_file_name, char * elements_file_name, char * bound_file_name[], char * time_file_name, char * stamps_file_name);

	virtual bool start_from_selected_time_layer (double time_value, char * solutions[]);
};
