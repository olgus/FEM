#include "Task_Mixed_Tr_Mesh.h"

double Task_Mixed_Tr_Mesh::function_f (int k_system, double * coordinates, int area)
{
	double x = coordinates[0];
	double y = coordinates[1];

	//return 1.0;
	return x + y;
	//return x * x + y * y + x + y - 4.0 + x * y;
	//return -202.0 + 400.0 * (y - 3 * x * x) + pow (1.0 - x, 2.0) + 100.0 * pow (y - x * x, 2.0);
	//return cos (x * y) * (1.0 + x * x + y * y);
}

double Task_Mixed_Tr_Mesh::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	double x = coordinates[0];
	double y = coordinates[1];

	//return 1.0;
	return x + y;
	//return x * x + y * y + x + y + x * y;
	//return pow (1.0 - x, 2.0) + 100.0 * pow (y - x * x, 2.0);
	//return cos (x * y);
}

double Task_Mixed_Tr_Mesh::function_Theta (int k_system, double * coordinates, int area, int boundary)
{
	double x = coordinates[0];
	double y = coordinates[1];
	double r = 0.0;

	switch (boundary)
	{
	case 2:
		r = x * sin (x * y);
		break;
	case 3:
		r = -x * sin (x * y);
		break;
	case 0:
		r = y * sin (x * y);
		break;
	case 1:
		r = -y * sin (x * y);
		break;
	}

	r = 0.0;
	return r;
}
