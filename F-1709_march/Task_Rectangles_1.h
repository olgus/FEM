#pragma once
#include "Prototype_Task.h"

class Task_Rectangle_Element_1 : public Task <Rectangular_Mesh>
{
private:
	double function_f (int k_system, double * coordinates, int area) override;
	double function_FCondition (int k_system, double * coordinates, int area, int boundary) override; 

public:
	Task_Rectangle_Element_1 () {};
	Task_Rectangle_Element_1 (const Task_Rectangle_Element_1 & tr1) {};
	~Task_Rectangle_Element_1 () {};
};