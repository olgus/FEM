#pragma once
#include "Prototype_Task.h"

// task that takes another task and uniform triangular mesh
// and interpolates the solution (whatever it is) on said mesh
class Task_2D_Mesh_Interpolation : public Task <Triangular_Mesh>
{
private:
	Task_pointer * task; // task to interpolate	
	int task_k_system; // solution to interpolate
	virtual void apply_boundary_conditions (int k_system) override; // no boundary conditions 														  
	virtual double function_f (int k_system, double * coordinates, int area) override;  // Fvector 
public:
	Task_2D_Mesh_Interpolation ();
	~Task_2D_Mesh_Interpolation ();
	void prepare (Task_pointer * Task, int K_system, double * c0, double * cN, int * N_axis); // set uniform mesh and task
};