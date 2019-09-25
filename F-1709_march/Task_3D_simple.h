#pragma once
#include "Prototype_Task.h"

template <class Type_Mesh>
class Task_3D : public Task <Type_Mesh>
{
protected:
	// override application of boundary conditions
	// for now
	// kinda stupid
	// but i don't really need 3D tasks, so whatever

	virtual void apply_first_boundary_conditions (int k_system) override; // applies first boundary conditions
	virtual void apply_second_boundary_conditions (int k_system) override; // applies second boundary conditions
	virtual void get_local_B (int k_system, int k_element, MathVector * B) override;
public:
	Task_3D ();
	Task_3D (const Task_3D & task);
	~Task_3D ();
};