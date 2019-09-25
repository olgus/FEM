#pragma once
#include <math.h>
#include "Task_pointer.h"
#include "Spline_1D_Pointer.h"

#define MAX_ITER_optimization 1000
#define ERROR_optimization 1e-15
#define ERROR_optimization_2 1e-15
#define ERROR_optimization_1D 1e-5
#define ERROR_optimization_1D_H 1e-10
#define ZERO_OPTIMIZATION 1e-15

namespace optimization
{
	double find_minimum_HJ (Task_pointer * task, int k_system, double * point, double * c0, double * cN);
	double find_maximum_HJ (Task_pointer * task, int k_system, double * point, double * c0, double * cN);

	double find_minimum_der_HJ (Task_pointer * task, int k_system, int k_var, double * point, double * c0, double * cN);
	double find_maximum_der_HJ (Task_pointer * task, int k_system, int k_var, double * point, double * c0, double * cN);
	
	double closest_value_point_GS (Task_pointer * task, int k_system, double * c0, double * c1, double val, double * sol_point);
	double closest_der_point_GS (Task_pointer * task, int k_system, int k_var, double * c0, double * c1, double val, double * sol_point);

	double min_Spline_GS (Spline_Pointer * task, double * sol_point);
	double max_Spline_GS (Spline_Pointer * task, double * sol_point);
}