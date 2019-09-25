#include "Integration.h"

Integrator_section_3D::Integrator_section_3D ()
{
	task = NULL;
	sections[0] = sections[1] = 0;
	variables[0] = variables[1] = 0;
	coord0[0] = coord0[1] = 0.0;
	coordN[0] = coordN[1] = 0.0;
	fixed_variable = 0;
	fixed_value = 0.0;
}

Integrator_section_3D::Integrator_section_3D (const Integrator_section_3D & integrator)
{
	task = NULL;
	if (integrator.task != NULL)
	{
		task = integrator.task;
	}

	for (int i = 0; i < 2; i++)
	{
		sections[i] = integrator.sections[i];
		variables[i] = integrator.variables[i];
		coord0[i] = integrator.coord0[i];
		coordN[i] = integrator.coordN[i];
	}
	fixed_variable = integrator.fixed_variable;
	fixed_value = integrator.fixed_value;
}

Integrator_section_3D::~Integrator_section_3D ()
{
}

void Integrator_section_3D::prepare (Task_pointer * Task, int K_var, double Value, double * Coord0, double * CoordN, int * Sections)
{
	task = Task;
	fixed_variable = K_var;
	int counter = 0;
	for (int i = 0; i < 3; i++)
	{
		if (i != fixed_variable)
		{
			variables[counter] = i;
			counter++;
		}
	}
	for (int i = 0; i < 2; i++)
	{
		coord0[i] = Coord0[i];
		coordN[i] = CoordN[i];
		sections[i] = Sections[i];
	}
	fixed_value = Value;
}

double Integrator_section_3D::function_integrate (int k_system, int k_func, int add_data, double * coordinates)
{
	double r;
	//r = coordinates[add_data] / pow (coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1] + coordinates[2] * coordinates[2], 1.5);
	//r /= 4.0;
	//r /= 3.14159265;
	//printf ("%lf ", r);

	switch (k_func)
	{
	case 0:
		// get solution value in coordinates
		task->get_solution_in_point (k_system, coordinates, &r);
		break;
	case 1:
		task->get_derivative (k_system, add_data, coordinates, &r);
		break;
	case 2:
		r = coordinates[add_data] / pow(coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1] + coordinates[2] * coordinates[2], 1.5);
		r /= 4.0;
		r /= 3.14159265;
		break;
	case 3:
		r = task->get_derivative_first (k_system, add_data, coordinates);
		break;
	case 4:
		r = 1.0;
		for (int i = 0; i < 2; i++)
		{
			r *= pow (coordinates[variables[i]], 6.0);
		}
		r += pow (coordinates[variables[0]], 5.0);
		r += pow (coordinates[variables[1]], 4.0);
		r += pow (coordinates[variables[0]], 3.0) * pow (coordinates[variables[1]], 3.0);
		break;
	default:
		r = 0.0;
	}
	//printf ("%lf\n", r);

	return r;
}

double Integrator_section_3D::integrate (int k_system, int k_func)
{
	double r = inner_integrate (k_system, k_func);
	return r;
}

double Integrator_section_3D::inner_integrate (int k_system, int k_func)
{
	// Gauss-4 
	int n_int_nodes = 4;
	// set w(i)
	MathVector * w = new MathVector (n_int_nodes);
	w->setElem (0, (18.0 + pow (30, 0.5)) / 36.0);
	w->setElem (1, (18.0 + pow (30, 0.5)) / 36.0);
	w->setElem (2, (18.0 - pow (30, 0.5)) / 36.0);
	w->setElem (3, (18.0 - pow (30, 0.5)) / 36.0);

	// set master point's cordinates
	MathVector * xi = new MathVector (n_int_nodes);
	xi->setElem (0, -pow (3.0 / 7.0 - 2.0 * pow (6.0 / 5.0, 0.5) / 7.0, 0.5));
	xi->setElem (1, pow (3.0 / 7.0 - 2.0 * pow (6.0 / 5.0, 0.5) / 7.0, 0.5));
	xi->setElem (2, -pow (3.0 / 7.0 + 2.0 * pow (6.0 / 5.0, 0.5) / 7.0, 0.5));
	xi->setElem (3, pow (3.0 / 7.0 + 2.0 * pow (6.0 / 5.0, 0.5) / 7.0, 0.5));

	double r = 0.0;
	double weight;
	double coordinates[3];
	coordinates[fixed_variable] = fixed_value;
	double h[2];
	for (int i = 0; i < 2; i++)
	{
		h[i] = (coordN[i] - coord0[i]) / sections[i];
	}

	double coord0s[2];
	double coordNs[2];
	double value;

	for (int i_sections = 0; i_sections < sections[0]; i_sections++)
	{
		coord0s[0] = coord0[0] + h[0] * i_sections;
		coordNs[0] = coord0[0] + h[0] * (i_sections + 1);
		for (int j_sections = 0; j_sections < sections[1]; j_sections++)
		{
			coord0s[1] = coord0[1] + h[1] * j_sections;
			coordNs[1] = coord0[1] + h[1] * (j_sections + 1);
			for (int i = 0; i < n_int_nodes; i++)
			{
				coordinates[variables[0]] = (coordNs[0] - coord0s[0]) * (xi->getElem (i) + 1.0) / 2.0 + coord0s[0];
				for (int j = 0; j < n_int_nodes; j++)
				{
					weight = w->getElem (i) * w->getElem (j);
					coordinates[variables[1]] = (coordNs[1] - coord0s[1]) * (xi->getElem (j) + 1.0) / 2.0 + coord0s[1];
					value = function_integrate (k_system, k_func, fixed_variable, coordinates);
					r += weight * value;
				}
			}
		}
	}
	double jac = 1.0;
	for (int k = 0; k < 2; k++)
	{
		jac *= coordN[k] - coord0[k];
	}
	jac /= 4.0;
	r *= jac;
	r /= sections[0];
	r /= sections[1];

	return r;
}

Integrator_section_2D::Integrator_section_2D ()
{
	task = NULL;
}

Integrator_section_2D::Integrator_section_2D (const Integrator_section_2D & integrator)
{
	delete task;
	if (integrator.task != NULL)
		task = integrator.task;

	sections = integrator.sections;
	for (int i = 0; i < 2; i++)
	{
		points[i] = integrator.points[i];
	}
}

Integrator_section_2D::~Integrator_section_2D ()
{
}

void Integrator_section_2D::prepare(Task_pointer * Task, int Var, double * Coord0, double * CoordN, int Sections)
{
	sections = Sections;

	points[0].set_point (Coord0);
	points[1].set_point (CoordN);

	task = Task;
	var = Var;
}

double Integrator_section_2D::function_integrate (int k_system, int k_func, int add_data, double * coordinates)
{
	double r;
	switch (k_func)
	{
	case 0:
		// get solution value in coordinates
		task->get_solution_in_point (k_system, coordinates, &r);
		break;
	case 1:
		// get solution value in coordinates
		r = task->get_derivative_first (k_system, add_data, coordinates);
		break;
	case 2:
		// get solution value in coordinates
		task->get_derivative (k_system, add_data, coordinates, &r);
		break;
	case 4:
		r = 0.0;
		for (int i = 0; i < 8; i++)
		{
			r += pow (coordinates[1], i);
		}
		//r += pow (coordinates[0], 5.0);
		//r += pow (coordinates[1], 4.0);
		//r += pow (coordinates[0], 3.0) * pow (coordinates[1], 3.0);
		break;
	case 5:
		r = 1.0;
		break;
	default:
		r = 0.0;
	}
	return r;
}

double Integrator_section_2D::integrate (int k_system, int k_func)
{
	// Gauss-4 
	int n_int_nodes = 4;
	// set w(i)
	MathVector * w = new MathVector (n_int_nodes);
	w->setElem (0, (18.0 + pow (30, 0.5)) / 36.0);
	w->setElem (1, (18.0 + pow (30, 0.5)) / 36.0);
	w->setElem (2, (18.0 - pow (30, 0.5)) / 36.0);
	w->setElem (3, (18.0 - pow (30, 0.5)) / 36.0);

	// set master point's cordinates
	MathVector * xi = new MathVector (n_int_nodes);
	xi->setElem (0, -pow (3.0 / 7.0 - 2.0 * pow (6.0 / 5.0, 0.5) / 7.0, 0.5));
	xi->setElem (1, pow (3.0 / 7.0 - 2.0 * pow (6.0 / 5.0, 0.5) / 7.0, 0.5));
	xi->setElem (2, -pow (3.0 / 7.0 + 2.0 * pow (6.0 / 5.0, 0.5) / 7.0, 0.5));
	xi->setElem (3, pow (3.0 / 7.0 + 2.0 * pow (6.0 / 5.0, 0.5) / 7.0, 0.5));

	double coordinates_d[2];
	double coordinates_1[2];
	double coordinates_2[2];

	double coord_1[2];
	double coord_2[2];
	double r = 0.0;

	points[0].get_point (coordinates_1);
	points[1].get_point (coordinates_2);

	double L = pow (pow (coordinates_2[0] - coordinates_1[0], 2.0) + pow (coordinates_2[1] - coordinates_1[1], 2.0), 0.5);
	double h = L / sections;

	double l;

	coord_1[0] = coordinates_1[0];
	coord_1[1] = coordinates_1[1];
	for (int i = 0; i < sections; i++)
	{
		coord_2[0] = (coordinates_2[0] - coordinates_1[0]) * h * (i + 1) / L + coordinates_1[0];
		coord_2[1] = (coordinates_2[1] - coordinates_1[1]) * h * (i + 1) / L + coordinates_1[1];

		for (int k = 0; k < n_int_nodes; k++)
		{
			// get coordinates through master coordinates
			l = xi->getElem (k) + 1.0;
			coordinates_d[0] = (coord_2[0] - coord_1[0]) * l / 2.0 + coord_1[0];
			coordinates_d[1] = (coord_2[1] - coord_1[1]) * l / 2.0 + coord_1[1];
			// calculate addition to integral value
			r += w->getElem (k) * function_integrate (k_system, k_func, var, coordinates_d);
		}
		coord_1[0] = coord_2[0];
		coord_1[1] = coord_2[1];
	}
	// pseudo jacobian 
	r *= L / 2.0;
	r /= sections;

	delete w;
	delete xi;
	return r;
}
