#include "Task_Rectangles_3.h"

double Task_Rectangle_Element_3::function_f (int k_system, double * coordinates, int area)
{
	double x = coordinates[0];
	double y = coordinates[1];

	return 10.0;
	//return x + y;
}

double Task_Rectangle_Element_3::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	double x = coordinates[0];
	double y = coordinates[1];

	return 10.0;
	//return x + y;
}
