#include "Cubic_Task.h"

double Cubic_Task::function_f (int k_system, double * coordinates, int area)
{
	double x = coordinates[0];
	double y = coordinates[1];
	double z = coordinates[2];

	return x * x + y * y + z * z + x * y * z - 6.0;
}

double Cubic_Task::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	double x = coordinates[0];
	double y = coordinates[1];
	double z = coordinates[2];

	return x * x + y * y + z * z + x * y * z;
}
