#include "Task_3D_Prism.h"

double Task_3D_Prism::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	double x = coordinates[0];
	double y = coordinates[1];
	double z = coordinates[2];

	//return 1.0;
	//return x + y + z;
	return x * z + y * z;
	//return x * x + y * y + z * z + x * y * z;
}

double Task_3D_Prism::function_f (int k_system, double * coordinates, int area)
{
	double x = coordinates[0];
	double y = coordinates[1];
	double z = coordinates[2];

	//return 1.0;
	//return x + y + z;
	return x * z + y * z;

	//return x * x + y * y + z * z + x * y * z - 6.0;
}
 