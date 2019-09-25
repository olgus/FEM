#pragma once
#include <stdio.h>
#include <math.h>
#include "MathVector.h"
#include "Point.h"
#include "Task_Pointer.h"

#define ZERO_integrator 1e-10

class Integrator_section_3D
{
private:
	Task_pointer * task; // FEM-task to integrate 

	double coord0[2]; // left boundaries
	double coordN[2]; // right boundaries
	int sections[2]; // amount of sections for integration on each axis
	int variables[2]; // variables to iterate through

	int fixed_variable; // fixed variable
	double fixed_value; // value of fixed variable

	double inner_integrate (int k_system, int k_func); // integrates on section
public:
	Integrator_section_3D ();
	Integrator_section_3D (const Integrator_section_3D & integrator);
	~Integrator_section_3D ();
	void prepare (Task_pointer * Task, int K_var, double Value, double * Coord0, double * CoordN, int * Sections); // sets parameters of integration and a task
	double function_integrate (int k_system, int k_func, int add_data, double * coordinates); // set a funciton to integrate
	double integrate (int k_system, int k_func); // integrates on section
};

class Integrator_section_2D
{
private:
	Task_pointer * task; // FEM-task to integrate 
	Point_2D points[2]; // points of the circuit
	int sections; // amount of sections for integration
	int var; // variable to derivate
public:
	Integrator_section_2D ();
	Integrator_section_2D (const Integrator_section_2D & integrator);
	~Integrator_section_2D ();

	void prepare (Task_pointer * Task, int Var, double * Coord0, double * CoordN, int Sections); // sets parameters of integration and a task
	
	double function_integrate (int k_system, int k_func, int add_data, double * coordinates); // set a funciton to integrate
	double integrate (int k_system, int k_func); // integrates on section
};
