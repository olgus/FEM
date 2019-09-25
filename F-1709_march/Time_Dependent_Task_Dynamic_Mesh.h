#pragma once
#include <vector>
#include <memory>

#include "Time_Dependent_Task.h"
#include "Dynamic_Mesh_Prototype.h"

// time dependent task on a dynamic mesh

template <class Type_Mesh>
class Time_Dependent_Task_Dynamic_Mesh <Type_Mesh> : Time_Dependent_Task <Type_Mesh>
{
private:
protected:
	std::vector <std::unique_ptr<Type_Mesh>> mesh; // keeps meshes for previous time layers, 
public:
	Time_Dependent_Task_Dynamic_Mesh ();
	Time_Dependent_Task_Dynamic_Mesh (const Time_Dependent_Task_Dynamic_Mesh & tdtdm);
	~Time_Dependent_Task_Dynamic_Mesh ();

	using Time_Dependent_Task <Type_Mesh>::get_solution_in_point;
	virtual bool get_solution_in_point (int k_system, int k_time_layer, double * coordinates, double * value) override;  // get a value in a point for specific time layer
	void next_time_layer_mesh (); // add switching mesh_pointer
};