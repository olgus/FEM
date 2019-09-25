#pragma once
#include "Prototype_Task.h"

class Task_Rectangle_Element_3 : public Task <Rectangular_Mesh_3>
{
private:
	double function_f (int k_system, double * coordinates, int area) override;
	double function_FCondition (int k_system, double * coordinates, int area, int boundary) override;
public:
	Task_Rectangle_Element_3 () {};
	Task_Rectangle_Element_3 (const Task_Rectangle_Element_3 & tr3) {};
	~Task_Rectangle_Element_3 () {};
};