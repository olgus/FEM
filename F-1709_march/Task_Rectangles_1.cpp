#include "Task_Rectangles_1.h"

double Task_Rectangle_Element_1::function_f (int k_system, double * coordinates, int area)
{
	double x = coordinates[0];
	double y = coordinates[1];

	//return 1.0;
	return x + y;
	//return coordinates[0] * coordinates[1];
}

double Task_Rectangle_Element_1::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	double x = coordinates[0];
	double y = coordinates[1];

	//return 1.0;
	return x + y;
}

