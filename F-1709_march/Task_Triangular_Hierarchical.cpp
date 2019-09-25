#include "Task_Triangular_Hierarchical.h"

double Task_Triangular_Hierarchical::function_f (int k_system, double * coordinates, int area)
{
	double x = coordinates[0];
	double y = coordinates[1];

	//return 1.0;
	//return x + y;
	return x * x + y * y - 4.0;
	//return -202.0 + 400.0 * (y - 3 * x * x) + pow (1.0 - x, 2.0) + 100.0 * pow (y - x * x, 2.0);
	//return 3.0 * cos (x + y);
}

double Task_Triangular_Hierarchical::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	double x = coordinates[0];
	double y = coordinates[1];

	//return 1.0;
	//return x + y;
	return x * x + y * y;
	//return pow (1.0 - x, 2.0) + 100.0 * pow (y - x * x, 2.0);
	//return cos (x + y);
}

double Task_Triangular_Hierarchical::function_Theta (int k_system, double * coordinates, int area, int boundary)
{
	return 0.0;
}
