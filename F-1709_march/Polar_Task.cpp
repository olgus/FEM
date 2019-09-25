#include "Polar_Task.h"

double Polar_Task::function_f (int k_system, double * coordinates, int area)
{
	double r = coordinates[0];
	double phi = coordinates[1];

	//return 10.0;
	//return r + phi - 1.0 / r;
	//return cos (r * phi) * (2.0 + pow (phi, 2.0) * r) + phi * sin (r * phi) / r;
	//return sin (phi / r) * (1.0 + (1.0 + pow (phi, 2.0)) / pow (r, 4.0)) - cos (phi / r) * phi / pow (r, 3.0);
	return cos (3.0 * r * phi) * (9.0 * pow (phi, 2.0) * (1.0 + 1.0 / pow (r, 2.0))) + sin (3.0 * r * phi) * 3.0 * r * phi;
	//return pow (phi, 2.0) / r - 2.0 / pow (r, 3.0);
}

double Polar_Task::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	double r = coordinates[0];
	double phi = coordinates[1];

	//return 1.0;
	//return r + phi;
	//return cos (r * phi);
	//return sin (phi / r);
	return cos (3.0 * r * phi);
	//return pow (phi, 2.0) / r;
}

double Polar_Task::function_Theta (int k_system, double * coordinates, int area, int boundary)
{
	double r = coordinates[0];
	double phi = coordinates[1];
	double val = 0.0;
	switch (boundary)
	{
	case 1:
		val = -1.0;
		break;
	case 11:
		val = 1.0;
		break;
	case 21:
		val = -1.0;
		break;
	case 31:
		val = 1.0;
		break;
	}
	return val;
}
