#include "Optimization.h"

namespace optimization
{
	double find_minimum_HJ (Task_pointer * task, int k_system, double * point, double * c0, double * cN)
	{
		int dim = task->get_mesh_dim ();

		// set starting point as the center of the mesh
		double * xcur = new double[dim];
		double * xstep = new double[dim];
		double min = 1e+20;
		for (int i = 0; i < dim; i++)
		{
			xcur[i] = (cN[i] + c0[i]) / 2.0;
			if (cN[i] - c0[i] < min)
				min = cN[i] - c0[i];
		}

		// set stating step as 1/10 of the smallest area side length
		double step = min / 10.0;
		double MIN_step = min / 1e+6;
		bool found = false;
		bool unlucky;
		int failed_steps;
		double value;

		double min_value;
		task->get_solution_in_point (k_system, xcur, &min_value);
		double global_min = min_value;

		for (int i = 0; i < MAX_ITER_optimization && !found; i++)// method iterations
		{
			unlucky = true;
			// while not found new point
			while (unlucky && step > MIN_step)
			{
				failed_steps = 0;
				for (int k = 0; k < dim; k++)
				{
					xstep[k] = xcur[k];
				}
				// iterate through variables
				for (int k = 0; k < dim; k++)
				{
					// make test +step 
					xstep[k] += step;
					// if there is a solution
					if (xstep[k] < cN[k] + ZERO_OPTIMIZATION && xstep[k] > c0[k] - ZERO_OPTIMIZATION)
					{
						if (task->get_solution_in_point (k_system, xstep, &value))
						{
							// if function value is smaller than previous
							if (value < min_value - ERROR_optimization_2)
							{
								// keep the step
								xcur[k] = xstep[k];
								min_value = value;
							}
							else
							{
								// else make test -step 
								xstep[k] = xcur[k];
								xstep[k] -= step;
								// if there is a solution
								if (xstep[k] < cN[k] + ZERO_OPTIMIZATION && xstep[k] > c0[k] - ZERO_OPTIMIZATION)
								{
									if (task->get_solution_in_point (k_system, xstep, &value))
									{
										// if function value is smaller than previous
										if (value < min_value - ERROR_optimization_2)
										{
											// keep the step
											xcur[k] = xstep[k];
											min_value = value;
										}
										// if not both, increase the amount of failed variables' steps
										else
										{
											failed_steps++;
											xstep[k] = xcur[k];
										}
									}
									else
									{
										failed_steps++;
										xstep[k] = xcur[k];
									}
								}
								else
								{
									failed_steps++;
									xstep[k] = xcur[k];
								}
							}
						}
						else
						{
							failed_steps++;
							xstep[k] = xcur[k];
						}
					}
					else
					{
						failed_steps++;
						xstep[k] = xcur[k];
					}
				}
				// if amount of failed steps equals to amount of variables, make the step smaller and start again
				if (failed_steps == dim)
				{
					step /= 2.0;
				}
				// else new point is found
				else
				{
					unlucky = false;
				}
			}

			// if the value hasn't changed or there is no new point, minimum has been found
			task->get_solution_in_point (k_system, xcur, &value);
			if (unlucky || fabs (value - global_min) < ERROR_optimization)
			{
				found = true;
			}
			else
			{
				min_value = value;
				global_min = min_value;
			}
		}

		for (int k = 0; k < dim; k++)
		{
			point[k] = xcur[k];
		}

		delete[] xcur;
		delete[] xstep;


		return min_value;
	}
	double find_maximum_HJ (Task_pointer * task, int k_system, double * point, double * c0, double * cN)
	{
		int dim = task->get_mesh_dim ();

		// set starting point as the center of the mesh
		double * xcur = new double[dim];
		double * xstep = new double[dim];
		double min = 1e+20;
		for (int i = 0; i < dim; i++)
		{
			xcur[i] = (cN[i] + c0[i]) / 2.0;
			if (cN[i] - c0[i] < min)
				min = cN[i] - c0[i];
		}

		// set stating step as 1/10 of the smallest area side length
		double step = min / 10.0;
		double MIN_step = min / 1e+6;
		bool found = false;
		bool unlucky;
		int failed_steps;
		double value;

		double max_value;
		task->get_solution_in_point (k_system, xcur, &max_value);
		double global_max = max_value;
		
		for (int i = 0; i < MAX_ITER_optimization && !found; i++)// method iterations
		{
			unlucky = true;
			// while not found new point
			while (unlucky && step > MIN_step)
			{
				failed_steps = 0;
				for (int k = 0; k < dim; k++)
				{
					xstep[k] = xcur[k];
				}
				// iterate through variables
				for (int k = 0; k < dim; k++)
				{
					// make test +step 
					xstep[k] += step;
					if (xstep[k] < cN[k] + ZERO_OPTIMIZATION && xstep[k] > c0[k] - ZERO_OPTIMIZATION)
					{

						// if there is a solution
						if (task->get_solution_in_point (k_system, xstep, &value))
						{

							// if function value is bigger than previous
							if (value > max_value + ERROR_optimization_2)
							{
								// keep the step
								xcur[k] = xstep[k];
								max_value = value;
							}
							else
							{
								// else make test -step 
								xstep[k] = xcur[k];
								xstep[k] -= step;
								// if there is a solution
								if (xstep[k] < cN[k] + ZERO_OPTIMIZATION && xstep[k] > c0[k] - ZERO_OPTIMIZATION)
								{
									if (task->get_solution_in_point (k_system, xstep, &value))
									{
										// if function value is bigger than previous
										if (value > max_value + ERROR_optimization_2)
										{
											// keep the step
											xcur[k] = xstep[k];
											max_value = value;
										}
										// if not both, increase the amount of failed variables' steps
										else
										{
											failed_steps++;
											xstep[k] = xcur[k];
										}
									}
									else
									{
										failed_steps++;
										xstep[k] = xcur[k];
									}
								}
								else
								{
									failed_steps++;
									xstep[k] = xcur[k];
								}
							}
						}
						else
						{
							failed_steps++;
							xstep[k] = xcur[k];
						}
					}
					else
					{
						failed_steps++;
						xstep[k] = xcur[k];
					}
				}
				// if amount of failed steps equals to amount of variables, make the step smaller and start again
				if (failed_steps == dim)
				{
					step /= 2.0;
				}
				// else new point is found
				else
				{
					unlucky = false;
				}
			}

			// if the value hasn't changed or there is no new point, minimum has been found
			task->get_solution_in_point (k_system, xcur, &value);
			if (unlucky || fabs (value - global_max) < ERROR_optimization)
			{
				found = true;
			}
			else
			{
				max_value = value;
				global_max = max_value;
			}
		}

		for (int k = 0; k < dim; k++)
		{
			point[k] = xcur[k];
		}

		delete[] xcur;
		delete[] xstep;

		return max_value;
	}

	double find_minimum_der_HJ (Task_pointer * task, int k_system, int k_var, double * point, double * c0, double * cN)
	{
		int dim = task->get_mesh_dim ();

		// set starting point as the center of the mesh
		double * xcur = new double[dim];
		double * xstep = new double[dim];
		double min = 1e+20;
		for (int i = 0; i < dim; i++)
		{
			xcur[i] = (cN[i] + c0[i]) / 2.0;
			if (cN[i] - c0[i] < min)
				min = cN[i] - c0[i];
		}

		// set stating step as 1/10 of the smallest area side length
		double step = min / 10.0;
		double MIN_step = min / 1e+6;
		bool found = false;
		bool unlucky;
		int failed_steps;
		double value;

		double min_value;
		task->get_derivative (k_system, k_var, xcur, &min_value);
		double global_min = min_value;

		for (int i = 0; i < MAX_ITER_optimization && !found; i++)// method iterations
		{
			unlucky = true;
			// while not found new point
			while (unlucky && step > MIN_step)
			{
				failed_steps = 0;
				for (int k = 0; k < dim; k++)
				{
					xstep[k] = xcur[k];
				}
				// iterate through variables
				for (int k = 0; k < dim; k++)
				{
					// make test +step 
					xstep[k] += step;
					// if there is a solution
					if (xstep[k] < cN[k] + ZERO_OPTIMIZATION && xstep[k] > c0[k] - ZERO_OPTIMIZATION)
					{
						if (task->get_derivative (k_system, k_var, xstep, &value))
						{
							// if function value is smaller than previous
							if (value < min_value - ERROR_optimization_2)
							{
								// keep the step
								xcur[k] = xstep[k];
								min_value = value;
							}
							else
							{
								// else make test -step 
								xstep[k] = xcur[k];
								xstep[k] -= step;
								// if there is a solution
								if (xstep[k] < cN[k] + ZERO_OPTIMIZATION && xstep[k] > c0[k] - ZERO_OPTIMIZATION)
								{
									if (task->get_derivative (k_system, k_var, xstep, &value))
									{
										// if function value is smaller than previous
										if (value < min_value - ERROR_optimization_2)
										{
											// keep the step
											xcur[k] = xstep[k];
											min_value = value;
										}
										// if not both, increase the amount of failed variables' steps
										else
										{
											failed_steps++;
											xstep[k] = xcur[k];
										}
									}
									else
									{
										failed_steps++;
										xstep[k] = xcur[k];
									}
								}
								else
								{
									failed_steps++;
									xstep[k] = xcur[k];
								}
							}
						}
						else
						{
							failed_steps++;
							xstep[k] = xcur[k];
						}
					}
					else
					{
						failed_steps++;
						xstep[k] = xcur[k];
					}
				}
				// if amount of failed steps equals to amount of variables, make the step smaller and start again
				if (failed_steps == dim)
				{
					step /= 2.0;
				}
				// else new point is found
				else
				{
					unlucky = false;
				}
			}

			// if the value hasn't changed or there is no new point, minimum has been found
			task->get_derivative (k_system, k_var, xcur, &value);
			if (unlucky || fabs (value - global_min) < ERROR_optimization)
			{
				found = true;
			}
			else
			{
				min_value = value;
				global_min = min_value;
			}
		}

		for (int k = 0; k < dim; k++)
		{
			point[k] = xcur[k];
		}

		delete[] xcur;
		delete[] xstep;


		return min_value;
	}
	double find_maximum_der_HJ (Task_pointer * task, int k_system, int k_var, double * point, double * c0, double * cN)
	{
		int dim = task->get_mesh_dim ();

		// set starting point as the center of the mesh
		double * xcur = new double[dim];
		double * xstep = new double[dim];
		double min = 1e+20;
		for (int i = 0; i < dim; i++)
		{
			xcur[i] = (cN[i] + c0[i]) / 2.0;
			if (cN[i] - c0[i] < min)
				min = cN[i] - c0[i];
		}

		// set stating step as 1/10 of the smallest area side length
		double step = min / 10.0;
		double MIN_step = min / 1e+6;
		bool found = false;
		bool unlucky;
		int failed_steps;
		double value;

		double max_value;
		task->get_derivative (k_system, k_var, xcur, &max_value);
		double global_max = max_value;

		for (int i = 0; i < MAX_ITER_optimization && !found; i++)// method iterations
		{
			unlucky = true;
			// while not found new point
			while (unlucky && step > MIN_step)
			{
				failed_steps = 0;
				for (int k = 0; k < dim; k++)
				{
					xstep[k] = xcur[k];
				}
				// iterate through variables
				for (int k = 0; k < dim; k++)
				{
					// make test +step 
					xstep[k] += step;
					if (xstep[k] < cN[k] + ZERO_OPTIMIZATION && xstep[k] > c0[k] - ZERO_OPTIMIZATION)
					{

						// if there is a solution
						if (task->get_derivative (k_system, k_var, xstep, &value))
						{

							// if function value is bigger than previous
							if (value > max_value + ERROR_optimization_2)
							{
								// keep the step
								xcur[k] = xstep[k];
								max_value = value;
							}
							else
							{
								// else make test -step 
								xstep[k] = xcur[k];
								xstep[k] -= step;
								// if there is a solution
								if (xstep[k] < cN[k] + ZERO_OPTIMIZATION && xstep[k] > c0[k] - ZERO_OPTIMIZATION)
								{
									if (task->get_derivative (k_system, k_var, xstep, &value))
									{
										// if function value is bigger than previous
										if (value > max_value + ERROR_optimization_2)
										{
											// keep the step
											xcur[k] = xstep[k];
											max_value = value;
										}
										// if not both, increase the amount of failed variables' steps
										else
										{
											failed_steps++;
											xstep[k] = xcur[k];
										}
									}
									else
									{
										failed_steps++;
										xstep[k] = xcur[k];
									}
								}
								else
								{
									failed_steps++;
									xstep[k] = xcur[k];
								}
							}
						}
						else
						{
							failed_steps++;
							xstep[k] = xcur[k];
						}
					}
					else
					{
						failed_steps++;
						xstep[k] = xcur[k];
					}
				}
				// if amount of failed steps equals to amount of variables, make the step smaller and start again
				if (failed_steps == dim)
				{
					step /= 2.0;
				}
				// else new point is found
				else
				{
					unlucky = false;
				}
			}

			// if the value hasn't changed or there is no new point, minimum has been found
			task->get_derivative (k_system, k_var, xcur, &value);
			if (unlucky || fabs (value - global_max) < ERROR_optimization)
			{
				found = true;
			}
			else
			{
				max_value = value;
				global_max = max_value;
			}
		}

		for (int k = 0; k < dim; k++)
		{
			point[k] = xcur[k];
		}

		delete[] xcur;
		delete[] xstep;

		return max_value;
	}

	double closest_value_point_GS (Task_pointer * task, int k_system, double * c0, double * c1, double val, double * sol_point)
	{
		double min, golden_ratio_coef, f1, f2,
			lenght_prev,    //длина интервала на предыдущем шаге
			lenght_cur,     //длина текущего интервала
			left, right;    //границы текущего интервала
		int k;
		bool l_r_sw;        //выбор из двух интервалов: [left,x2]||[x1,right]
		golden_ratio_coef = (3.0 - pow (5.0, 0.5)) / 2.0;
		double a1, a2;

		if (c0[0] > c1[0])
		{
			left = c0[0];
			c0[0] = c1[0];
			c1[0] = left;
			left = c0[1];
			c0[1] = c1[1];
			c1[1] = left;
		}

		left = 0.0;
		right = 1.0;
		a1 = left + golden_ratio_coef;
		a2 = right - golden_ratio_coef;
		//
		double value;
		double point[2];
		point[0] = a1 * (c1[0] - c0[0]) + c0[0];
		point[1] = a1 * (c1[1] - c0[1]) + c0[1];
		task->get_solution_in_point (k_system, point, &value);
		f1 = fabs (value - val);
		point[0] = a2 * (c1[0] - c0[0]) + c0[0];
		point[1] = a2 * (c1[1] - c0[1]) + c0[1];
		task->get_solution_in_point (k_system, point, &value);
		f2 = fabs (value - val);

		if (f1 < f2)        //выбор промежутка, на котором находится экстремум
			l_r_sw = true;
		else
			l_r_sw = false;

		lenght_prev = 1.0;
		int maxIter = 5000;
		for (k = 0; k < MAX_ITER_optimization && right - left > ERROR_optimization_1D_H; k++)
		{
			if (l_r_sw) //вычисление новой точки в интервале, содержащем экстремум
			{
				a1 = left + golden_ratio_coef * (right - left);
				point[0] = a1 * (c1[0] - c0[0]) + c0[0];
				point[1] = a1 * (c1[1] - c0[1]) + c0[1];
				task->get_solution_in_point (k_system, point, &value);
				f1 = fabs (value - val);
			}
			else
			{
				a2 = right - golden_ratio_coef * (right - left);
				point[0] = a2 * (c1[0] - c0[0]) + c0[0];
				point[1] = a2 * (c1[1] - c0[1]) + c0[1];
				task->get_solution_in_point (k_system, point, &value);
				f2 = fabs (value - val);
			}

			if (f1 < f2)    //сокращение интервала, определение будущего направления движения
			{
				right = a2;
				a2 = a1;
				f2 = f1;
				l_r_sw = true;
			}
			else
			{
				left = a1;
				a1 = a2;
				f1 = f2;
				l_r_sw = false;
			}

			lenght_cur = right - left;
			lenght_prev = lenght_cur;
		}

		min = (left + right) / 2;

		point[0] = min * (c1[0] - c0[0]) + c0[0];
		point[1] = min * (c1[1] - c0[1]) + c0[1];
		task->get_solution_in_point (k_system, point, &value);
		double pr = fabs (value - val);

		for (int i = 0; i < 2; i++)
		{
			sol_point[i] = point[i];
		}
		return pr;
	}
	double closest_der_point_GS (Task_pointer * task, int k_system, int k_var, double * c0, double * c1, double val, double * sol_point)
	{
		double min, golden_ratio_coef, f1, f2,
			lenght_prev,    //длина интервала на предыдущем шаге
			lenght_cur,     //длина текущего интервала
			left, right;    //границы текущего интервала
		int k;
		bool l_r_sw;        //выбор из двух интервалов: [left,x2]||[x1,right]
		golden_ratio_coef = (3.0 - pow (5.0, 0.5)) / 2.0;
		double a1, a2;

		if (c0[0] > c1[0])
		{
			left = c0[0];
			c0[0] = c1[0];
			c1[0] = left;
			left = c0[1];
			c0[1] = c1[1];
			c1[1] = left;
		}

		left = 0.0;
		right = 1.0;
		a1 = left + golden_ratio_coef;
		a2 = right - golden_ratio_coef;
		//
		double value;
		double point[2];
		point[0] = a1 * (c1[0] - c0[0]) + c0[0];
		point[1] = a1 * (c1[1] - c0[1]) + c0[1];
		task->get_derivative (k_system, k_var, point, &value);
		f1 = fabs (value - val);
		point[0] = a2 * (c1[0] - c0[0]) + c0[0];
		point[1] = a2 * (c1[1] - c0[1]) + c0[1];
		task->get_derivative (k_system, k_var, point, &value);
		f2 = fabs (value - val);

		if (f1 < f2)        //выбор промежутка, на котором находится экстремум
			l_r_sw = true;
		else
			l_r_sw = false;

		lenght_prev = 1.0;
		int maxIter = 100;
		for (k = 0; k < maxIter && right - left > ERROR_optimization_1D; k++)
		{
			if (l_r_sw) //вычисление новой точки в интервале, содержащем экстремум
			{
				a1 = left + golden_ratio_coef * (right - left);
				point[0] = a1 * (c1[0] - c0[0]) + c0[0];
				point[1] = a1 * (c1[1] - c0[1]) + c0[1];
				task->get_derivative (k_system, k_var, point, &value);
				f1 = fabs (value - val);
			}
			else
			{
				a2 = right - golden_ratio_coef * (right - left);
				point[0] = a2 * (c1[0] - c0[0]) + c0[0];
				point[1] = a2 * (c1[1] - c0[1]) + c0[1];
				task->get_derivative (k_system, k_var, point, &value);
				f2 = fabs (value - val);
			}

			if (f1 < f2)    //сокращение интервала, определение будущего направления движения
			{
				right = a2;
				a2 = a1;
				f2 = f1;
				l_r_sw = true;
			}
			else
			{
				left = a1;
				a1 = a2;
				f1 = f2;
				l_r_sw = false;
			}

			lenght_cur = right - left;
			lenght_prev = lenght_cur;
		}

		min = (left + right) / 2;

		point[0] = min * (c1[0] - c0[0]) + c0[0];
		point[1] = min * (c1[1] - c0[1]) + c0[1];
		task->get_derivative (k_system, k_var, point, &value);
		double pr = fabs (value - val);

		for (int i = 0; i < 2; i++)
		{
			sol_point[i] = point[i];
		}
		return pr;
	}
	
	double min_Spline_GS (Spline_Pointer * task, double * sol_point)
	{
		double min, golden_ratio_coef, f1, f2,
			lenght_prev,    //длина интервала на предыдущем шаге
			lenght_cur,     //длина текущего интервала
			left, right;    //границы текущего интервала
		int k;
		bool l_r_sw;        //выбор из двух интервалов: [left,x2]||[x1,right]
		golden_ratio_coef = (3.0 - pow (5.0, 0.5)) / 2.0;
		double a1, a2;

		double c0[1], c1[1];
		task->get_boundaries (c0, c1);

		left = 0.0;
		right = 1.0;
		a1 = left + golden_ratio_coef;
		a2 = right - golden_ratio_coef;
		//
		double value;
		double point[1];
		point[0] = a1 * (c1[0] - c0[0]) + c0[0];
		task->get_solution_in_point (point, &value);
		f1 = value;
		point[0] = a2 * (c1[0] - c0[0]) + c0[0];
		task->get_solution_in_point (point, &value);
		f2 = value;

		if (f1 < f2)        //выбор промежутка, на котором находится экстремум
			l_r_sw = true;
		else
			l_r_sw = false;

		lenght_prev = 1.0;
		int maxIter = 100;
		for (k = 0; k < maxIter && right - left > ERROR_optimization_1D; k++)
		{
			if (l_r_sw) //вычисление новой точки в интервале, содержащем экстремум
			{
				a1 = left + golden_ratio_coef * (right - left);
				point[0] = a1 * (c1[0] - c0[0]) + c0[0];
				task->get_solution_in_point (point, &value);
				f1 = value;
			}
			else
			{
				a2 = right - golden_ratio_coef * (right - left);
				point[0] = a2 * (c1[0] - c0[0]) + c0[0];
				task->get_solution_in_point (point, &value);
				f2 = value;
			}

			if (f1 < f2)    //сокращение интервала, определение будущего направления движения
			{
				right = a2;
				a2 = a1;
				f2 = f1;
				l_r_sw = true;
			}
			else
			{
				left = a1;
				a1 = a2;
				f1 = f2;
				l_r_sw = false;
			}

			lenght_cur = right - left;
			lenght_prev = lenght_cur;
		}

		min = (left + right) / 2;

		point[0] = min * (c1[0] - c0[0]) + c0[0];
		task->get_solution_in_point (point, &value);

		double pr = value;
		sol_point[0] = c0[0];
		return pr;
	}
	double max_Spline_GS (Spline_Pointer * task, double * sol_point)
	{
		double max, golden_ratio_coef, f1, f2,
			lenght_prev,    //длина интервала на предыдущем шаге
			lenght_cur,     //длина текущего интервала
			left, right;    //границы текущего интервала
		int k;
		bool l_r_sw;        //выбор из двух интервалов: [left,x2]||[x1,right]
		golden_ratio_coef = (3.0 - pow (5.0, 0.5)) / 2.0;
		double a1, a2;

		double c0[1], c1[1];
		task->get_boundaries (c0, c1);

		left = 0.0;
		right = 1.0;
		a1 = left + golden_ratio_coef;
		a2 = right - golden_ratio_coef;
		//
		double value;
		double point[1];
		point[0] = a1 * (c1[0] - c0[0]) + c0[0];
		task->get_solution_in_point (point, &value);
		f1 = value;
		point[0] = a2 * (c1[0] - c0[0]) + c0[0];
		task->get_solution_in_point (point, &value);
		f2 = value;

		if (f1 > f2)        //выбор промежутка, на котором находится экстремум
			l_r_sw = true;
		else
			l_r_sw = false;

		lenght_prev = 1.0;
		int maxIter = 100;
		for (k = 0; k < maxIter && right - left > ERROR_optimization_1D; k++)
		{
			if (l_r_sw) //вычисление новой точки в интервале, содержащем экстремум
			{
				a1 = left + golden_ratio_coef * (right - left);
				point[0] = a1 * (c1[0] - c0[0]) + c0[0];
				task->get_solution_in_point (point, &value);
				f1 = value;
			}
			else
			{
				a2 = right - golden_ratio_coef * (right - left);
				point[0] = a2 * (c1[0] - c0[0]) + c0[0];
				task->get_solution_in_point (point, &value);
				f2 = value;
			}

			if (f1 > f2)    //сокращение интервала, определение будущего направления движения
			{
				right = a2;
				a2 = a1;
				f2 = f1;
				l_r_sw = true;
			}
			else
			{
				left = a1;
				a1 = a2;
				f1 = f2;
				l_r_sw = false;
			}

			lenght_cur = right - left;
			lenght_prev = lenght_cur;
		}

		max = (left + right) / 2;

		point[0] = max * (c1[0] - c0[0]) + c0[0];
		task->get_solution_in_point (point, &value);

		double pr = value;
		sol_point[0] = c1[0];
		return pr;
	}
}