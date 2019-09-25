#include "Task_Vector_2D.h"

void Task_Vector_2D::get_local_B (int k_system, int k_element, MathVector * B)
{
	// B(i) = integral (basis_function[i] * function_f)
	int n_funcs = mesh_pointer->get_amount_non_zero_functions (k_element);
	int area = mesh_pointer->get_area (k_element);

	double jac = 0.0; // jacobian
	int dim = mesh_pointer->get_dimentionality ();
	int n_integr_points = mesh.elements[k_element]->amount_of_integration_points ();
	double * weigths = new double[n_integr_points];
	double ** points = new double *[n_integr_points];
	for (int i = 0; i < n_integr_points; i++)
	{
		points[i] = new double[dim];
	}

	mesh.elements[k_element]->integration_points (points, weigths, &jac);

	double val;
	double f_value[2], basis_func_value[2];
	double sp;
	for (int i = 0; i < n_funcs; i++)
	{
		// get vector of integrals of scalar product f * function_i
		val = 0.0;
		// get theta values * basis_functions in those points
		for (int k = 0; k < n_integr_points; k++)
		{
			// get basis function value
			mesh.elements[k_element]->get_basis_function_value (i, points[k], basis_func_value);
			// get f value
			function_f (k_system, points[k], area, f_value);
			// get their scalar product
			sp = basis_func_value[0] * f_value[0] + basis_func_value[1] * f_value[1];
			// sum them multiplying by weigths
			val += sp * weigths[k];
		}
		// multiply by jac
		val *= jac;
		B->setElem (i, val);
	}

	delete[] weigths;
	for (int i = 0; i < n_integr_points; i++)
	{
		delete[] points[i]; 
	}
	delete[] points;
}

void Task_Vector_2D::get_FCondition_edge (int k_system, int k_element, int n1, int n2, int * functions, int boundary, MathVector * FCondition)
{
	int n_funcs = FCondition->getSize ();
	int area = mesh_pointer->get_area (k_element);

	double jac = 0.0; // jacobian
	int dim = mesh_pointer->get_dimentionality ();
	int n_integr_points = mesh.elements[k_element]->edge_amount_of_integration_points ();
	double * weigths = new double[n_integr_points];
	double ** points = new double *[n_integr_points];
	for (int i = 0; i < n_integr_points; i++)
	{
		points[i] = new double[dim];
	}

	// get integration points for the edge from the element
	mesh.elements[k_element]->edge_integration_points (*mesh_pointer, n1, n2, points, weigths, &jac);
	// get edge-vector and normal vector for it
	double e[2];
	double n[2];
	mesh.get_node_coordinates (n1, e);
	mesh.get_node_coordinates (n2, n);
	e[0] -= n[0];
	e[1] -= n[1];
	if (fabs (e[1]) > 1e-7)
	{
		n[0] = sqrt (1.0 / (1.0 + pow (e[0] / e[1], 2.0)));
		n[1] = -n[0] * e[0] / e[1];
	}
	else
	{
		n[0] = 0.0;
		n[1] = 1.0;
	}
	if (boundary % 2 == 0)
	{
		n[0] *= -1;
		n[1] *= -1;
	}

	// make matrix 
	Matrix * A = new Matrix (n_funcs, n_funcs);
	MathVector * B = new MathVector (n_funcs);

	double val, sp;
	double value[2], basis_func_value[2];
	double vp1, vp2;
	for (int i = 0; i < n_funcs; i++)
	{
		for (int j = 0; j < n_funcs; j++)
		{
			val = 0.0;
			// get vector of integrals scalar product function_i * function_j
			for (int k = 0; k < n_integr_points; k++)
			{
				// get basis function value
				mesh.elements[k_element]->get_basis_function_value (functions[i], points[k], basis_func_value);
				// get basis function value
				mesh.elements[k_element]->get_basis_function_value (functions[j], points[k], value);
				// get vector products for them
				vp1 = basis_func_value[0] * n[1] - basis_func_value[1] * n[0];
				vp2 = value[0] * n[1] - value[1] * n[0];
				// get the scalar product
				sp = vp1 * vp2;
				// sum them multiplying by weigths
				val += sp * weigths[k];
			}
			// multiply by jac
			val *= jac;
			A->setElem (i, j, val);
		}

		val = 0.0;
		// get Fcondition values * basis_functions in those points
		for (int k = 0; k < n_integr_points; k++)
		{
			// get basis function value
			mesh.elements[k_element]->get_basis_function_value (functions[i], points[k], basis_func_value);
			// get first condition value
			function_FCondition (k_system, points[k], area, boundary, value);
			// get vector products for them
			vp1 = basis_func_value[0] * n[1] - basis_func_value[1] * n[0];
			vp2 = value[0] * n[1] - value[1] * n[0];
			// get the scalar product
			sp = vp1 * vp2;
			// sum them multiplying by weigths
			val += sp * weigths[k];
		}
		// multiply by jac
		val *= jac;
		B->setElem (i, val);
	}
	// solve it
	A->solve (FCondition, *B);
	//FCondition->Print ();
	delete A;
	delete B;
	delete[] weigths;
	for (int i = 0; i < n_integr_points; i++)
	{
		delete[] points[i];
	}
	delete[] points;
}

void Task_Vector_2D::get_Theta_edge (int k_system, int k_element, int n1, int n2, int * functions, int boundary, MathVector * Theta)
{
	// amount of basis functions on the edge is Theta size, since it was established functions earlier
	int n_funcs = Theta->getSize ();
	int area = mesh_pointer->get_area (k_element);

	double jac = 0.0; // jacobian
	int dim = mesh_pointer->get_dimentionality ();
	int n_integr_points = mesh.elements[k_element]->edge_amount_of_integration_points ();
	double * weigths = new double[n_integr_points];
	double ** points = new double *[n_integr_points];
	for (int i = 0; i < n_integr_points; i++)
	{
		points[i] = new double[dim];
	}

	// get integration points for the edge from the element
	mesh.elements[k_element]->edge_integration_points (*mesh_pointer, n1, n2, points, weigths, &jac);
	
	// get normal
	double e[2];
	double n[2];
	mesh.get_node_coordinates (n1, e);
	mesh.get_node_coordinates (n2, n);
	e[0] -= n[0];
	e[1] -= n[1];
	if (fabs (e[1]) > 1e-7)
	{
		n[0] = sqrt (1.0 / (1.0 + pow (e[0] / e[1], 2.0)));
		n[1] = -n[0] * e[0] / e[1];
	}
	else
	{
		n[0] = 0.0;
		n[1] = 1.0;
	}
	if (boundary % 2 == 0)
	{
		n[0] *= -1;
		n[1] *= -1;
	}

	double val;
	double integral_function;
	double basis_function_value[2];
	for (int i = 0; i < n_funcs; i++)
	{
		val = 0.0;
		// get theta * normal * basis_functions in those points
		for (int j = 0; j < n_integr_points; j++)
		{
			// sum them multiplying by weigths
			integral_function = function_Theta (k_system, points[j], area, boundary);
			// get basis function
			mesh.elements[k_element]->get_basis_function_value (functions[i], points[j], basis_function_value);
			integral_function *= (-basis_function_value[0] * n[1] + basis_function_value[1] * n[0]);
			val += integral_function * weigths[j];
		}
		// multiply by jac
		val *= jac;
		Theta->setElem (i, val);
	}

	delete[] weigths;
	for (int i = 0; i < n_integr_points; i++)
	{
		delete[] points[i];
	}
	delete[] points;
}

void Task_Vector_2D::function_f (int k_system, double * coordinates, int area, double * f_value)
{
	double x = coordinates[0];
	double y = coordinates[1];

	//f_value[0] = y;
	//f_value[1] = x;

	f_value[0] = - x - y - 1.0;
	f_value[1] = y + x + 1.0;

	//f_value[0] = x * x + y * y - 1.0;
	//f_value[1] = x * y;

	//f_value[0] = exp (x * y) - 9.0 * x * x * y * y - x * x * exp (x * y);
	//f_value[1] = -x * x * x  * y * y * y + 6.0 * x * y * y * y  + exp (x * y) * (1.0 + x * y);
}

void Task_Vector_2D::function_FCondition (int k_system, double * coordinates, int area, int boundary, double * value)
{
	double x = coordinates[0];
	double y = coordinates[1];

	//value[0] = y;
	//value[1] = x;

	value[0] = - x - y - 1.0;
	value[1] = y + x + 1.0;

	//value[0] = x * x + y * y;
	//value[1] = x * y;

	//value[0] = exp (x * y);
	//value[1] = - x * x * x  * y * y * y;
}

double Task_Vector_2D::Lambda (int k_system, int k_function, int area)
{
	return 1.0;
}

double Task_Vector_2D::Gamma (int k_system, int k_function, int area)
{
	return 1.0;
}

double Task_Vector_2D::function_Theta (int k_system, double * coordinates, int area, int boundary)
{
	// theta is scalar
	// dAy/dx - dAx/dy
	double x = coordinates[0];
	double y = coordinates[1];
	double r = 0.0;

	r = - x * exp (x * y) - 3.0 * x * x * y * y * y;

	return r;
}

bool Task_Vector_2D::get_solution_in_point (int k_system, double * coordinates, double * value)
{
	// find element that contains the point 
	int p_element = mesh_pointer->point_inside (coordinates);
	if (p_element == -1)
	{
		value[0] = 0.0;
		value[1] = 0.0;
		return false;
	}

	// get local functions' values of the element
	int func_amount = mesh_pointer->get_amount_non_zero_functions (p_element);
	int dim = mesh_pointer->get_dimentionality ();

	Matrix * local_func = new Matrix (func_amount, dim);
	mesh_pointer->get_local_function_values (p_element, coordinates, local_func);
	//local_func->Print ();

	// make vector of solutions in those nodes
	MathVector * solution = new MathVector (func_amount);

	int global_function;
	for (int i = 0; i < func_amount; i++)
	{
		global_function = get_function_global_number (p_element, i);
		if (global_function != -1)
		{
			solution->setElem (i, previous_time_layers_solutions[k_system][0].getElem (global_function));
		}
	}

	// multiply them by the solution's values in those nodes
	// make vector of solutions in those nodes
	MathVector * vector_solution = new MathVector (func_amount);

	for (int j = 0; j < dim; j++)
	{
		for (int i = 0; i < func_amount; i++)
		{
			vector_solution->setElem (i, local_func->Elem (i, j));
		}
		value[j] = vector_solution->Scalar_Product (*solution);
	}
	delete local_func;
	delete solution;
	delete vector_solution;

	return true;
}

void Task_Vector_2D::save_evenly_2D (int k_system, double * h, char * file_dest, bool keep_points)
{
	FILE * file = fopen (file_dest, "w");
	int dim = mesh_pointer->get_dimentionality ();
	double * coordinates0 = new double[dim];
	mesh_pointer->get_0_boundaries (coordinates0);
	double * coordinatesN = new double[dim];
	mesh_pointer->get_N_boundaries (coordinatesN);
	double * coordinates = new double[dim];
	bool exists = false;
	double value[2];
	int n_amount = 0;

	// starting point
	for (int i = 0; i < dim; i++)
	{
		coordinates[i] = coordinates0[i];
	}

	for (int i = 0; coordinates[0] < coordinatesN[0]; i++)
	{
		coordinates[0] = coordinates0[0] + i * h[0];
		coordinates[1] = coordinates0[1];
		for (int j = 0; coordinates[1] < coordinatesN[1]; j++)
		{
			coordinates[1] = coordinates0[1] + j * h[1];
			exists = get_solution_in_point (k_system, coordinates, value);
			if (keep_points)
			{
				if (exists)
				{
					for (int k = 0; k < dim; k++)
					{
						fprintf (file, "%.16lf ", coordinates[k]);
					}
					fprintf (file, "%.16lf ", value[0]);
					fprintf (file, "%.16lf\n", value[1]);

					n_amount++;
				}
			}
			else
			{
				fprintf (file, "%.16lf ", value[0]);
				fprintf (file, "%.16lf\n", value[1]);
			}
		}
	}
	fclose (file);

	delete[] coordinates;
	delete[] coordinates0;
	delete[] coordinatesN;

	if (keep_points)
	{
		file = fopen ("Result Files Extra//points_amount.txt", "w");
		fprintf (file, "%i", n_amount);
		fclose (file);
	}

	printf ("\trequested solution is saved\n");
}
