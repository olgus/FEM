#include "Cartesian_Test_Task.h"

double Cartesian_Test_Task::function_f (int k_system, double * coordinates, int area)
{
	double x = coordinates[0];
	double y = coordinates[1];

	return x + y;
	//return 2.0;
	//return 2.0 * x - y;
	//return coordinates[0] * coordinates[1];
	//return -2.0 * pow (coordinates[1], 3.0) -
	//	6.0 * pow (coordinates[0], 2.0) * coordinates[1] +
	//	pow (coordinates[0], 2.0) * pow (coordinates[1], 3.0);
	//return -202.0 + 400.0 * (y - 3 * x * x) + pow (1.0 - x, 2.0) + 100.0 * pow (y - x * x, 2.0);
	//return -6.0 * x - 6.0 *y + pow (x, 3) + pow (y, 3);
	//return -4.0 + pow (x, 2) + pow (y, 2);
	//return 3.0 * cos (x + y);

	//return 0.0;
}

double Cartesian_Test_Task::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	double x = coordinates[0];
	double y = coordinates[1];

	return x + y;
	//return 2.0;
	//return 2.0 * x - y;
	//return coordinates[0] * coordinates[1];
	//return pow (coordinates[0], 2.0) * pow (coordinates[1], 3.0);
	//return pow (1.0 - x, 2.0) + 100.0 * pow (y - x * x, 2.0);
	//return pow (x, 3) + pow (y, 3);
	//return pow (x, 2) + y;
	//return cos (x + y);

	//switch (boundary)
	//{
	//case (0):
	//	return 1.0;
	//case (1):
	//	return 0.0;
	//case (2):
	//	return 1.0;
	//case (3):
	//	return 0.0;
	//}

	//return 0.0;
}

double Cartesian_Test_Task::function_Theta (int k_system, double * coordinates, int area, int boundary)
{
	return 1.0;
}

double Cartesian_Test_Task::Lambda (int k_system, int k_function, int area)
{
	double r = 1.0;
	//switch (area)
	//{
	//case (1):
	//	r = 1.0;
	//	break;
	//case (2):
	//	//r = 990.02;
	//	//r = 166.94;
	//	r = 1000.0;
	//	break;
	//case (3):
	//	r = 2.0;
	//	break;
	//default:
	//	r = 0.0;
	//}
	return r;
}

double Cartesian_Test_Task::Gamma (int k_system, int k_function, int area)
{
	return 1.0;
}
