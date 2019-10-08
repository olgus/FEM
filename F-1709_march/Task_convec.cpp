#include "Task_convec.h"

void Task_convec::get_min_max (int k_system, double * min, double * max)
{
	double mins = 1e+20, maxs = -1e+20, e;
	for (int i = 0; i < mesh_pointer->get_n_nodes (); i++)
	{
		e = previous_time_layers_solutions[k_system][0].getElem (i);
		if (e < mins)
			mins = e;
		if (e > maxs)
			maxs = e;
	}
	*min = mins;
	*max = maxs;
	//printf ("%i %lf %lf\n", k_system, *min, *max);
}

bool Task_convec::solve_task (int solver_param[][5])
{
	{
		if (!prepared)
		{
			printf ("ERROR: missing preparation\n");
			return false;
		}

		// set ffb flags
		set_ffb_flags ();

		// set Pr, Gr splines
		if (PRGRconst == 0)
		{
			double extra[] = {1e-10, 1e-10};
			double step = 1.0;
			spline_Pr = new Spline_1D_Smooth ();
			spline_Pr->prepare ("Source Files//TC//spline_Pr.txt", extra, step);
			spline_Pr->solve_task ();

			spline_Gr = new Spline_1D_Smooth ();
			spline_Gr->prepare ("Source Files//TC//spline_Gr.txt", extra, step);
			spline_Gr->solve_task ();
		}

		// n_systems is set in preparation to 1 by default
		non_linear_layers_solutions = new MathVector[n_systems];
		MathVector * prev_non_linear_layers_solutions = new MathVector[n_systems];

		for (int k_system = 0; k_system < n_systems; k_system++)
		{
			non_linear_layers_solutions[k_system].setSize (eq_matrixes[k_system].Size ()); // default starting solution point is 0
			prev_non_linear_layers_solutions[k_system].setSize (eq_matrixes[k_system].Size ()); // default starting solution point is 0
		}

		double sum_non_linear_dif;
		double sum_non_linear_discr;
		double * non_linear_dif = new double[n_systems];
		double * non_linear_discr = new double[n_systems];
		double w_rel = 1.0;
		double norm;

		// 0 non_linear solution is a start solution
		if (start_cond)
			set_starting_conditions ();

		// states
		state = new int[N_functions];

		int solver_iterations;
		bool found = false;

		// TODO time cycle
		bool stable = false; // indicator of reaching stable condition
		for (current_time_layer = 1; current_time_layer < n_time_layers && !stable; current_time_layer++)
		{
			if (current_time_layer % 2 == 0)
				system ("cls");
			printf ("time layer: %.7lf\n", time_layers[current_time_layer]);
			found = false;
			for (int k_system = 0; k_system < n_systems; k_system++)
			{
				non_linear_layers_solutions[k_system].Copy ((previous_time_layers_solutions[k_system][0]));
			}

			reset_states ();

			for (int k_nl_iter = 0; k_nl_iter < MAX_TC_NONLINEAR_ITER && !found; k_nl_iter++)
			{
				K_NON_LINEAR_ITER = k_nl_iter;
				if ((k_nl_iter + 2) % 15 == 0)
				{
					system ("cls");
					printf ("time layer: %.7lf\n", time_layers[current_time_layer]);
				}

				printf ("k_nl_iter: %i\n", k_nl_iter);
				sum_non_linear_dif = 0.0;
				sum_non_linear_discr = 0.0;

				for (int k_system = 0; k_system < n_systems; k_system++)
				{
					// save prev non linear solution
					prev_non_linear_layers_solutions[k_system].Copy (non_linear_layers_solutions[k_system]);
				}
#pragma omp parallel for num_threads (n_systems) 
				for (int k_system = 0; k_system < n_systems; k_system++)
				{
					// set starting point for the solver
					eq_matrixes[k_system].set_starting_point (&(non_linear_layers_solutions[k_system])); // set last layer solution as x0 for solver
					// build system
					build_system (k_system);
					// apply boundary conditions
					apply_boundary_conditions (k_system); // apply boundary conditions
				}

				// solve it
				//eq_matrixes[0].fprint ("Test//matrix.txt");
				for (int k_system = 0; k_system < n_systems; k_system++)
				{
					solver_iterations = eq_matrixes[k_system].solve (solver_param[k_system]);
				}

#pragma omp parallel for num_threads (n_systems)
				for (int k_system = 0; k_system < n_systems; k_system++)
				{
					// TODO relax
					// for now just copy the solution from the matrix
					//eq_matrixes[k_system].get_solution (&(non_linear_layers_solutions[k_system]));
					////relax (k_system, 1.0, *prev_non_linear_layers_solutions);
					//if (k_system == 0)
					//	non_linear_discr[k_system] = relax (k_system, 1.0, *prev_non_linear_layers_solutions);
					//else
					//	non_linear_discr[k_system] = relax (k_system, 0.3, *prev_non_linear_layers_solutions);
					
					// non_linear_discr[k_system] = relax (k_system, *prev_non_linear_layers_solutions);
				
					if (k_system == 0)
						non_linear_discr[k_system] = relax (k_system, 0.1, 1.0, prev_non_linear_layers_solutions[k_system]);
					else
						non_linear_discr[k_system] = relax (k_system, 0.1, 0.8, prev_non_linear_layers_solutions[k_system]);


					// non_linear_layers_solutions contains entire _new_ solution 
					prev_non_linear_layers_solutions[k_system].Substract (non_linear_layers_solutions[k_system]);
					norm = non_linear_layers_solutions[k_system].Norm ();
					if (norm > 1e-15)
						non_linear_dif[k_system] = prev_non_linear_layers_solutions[k_system].Norm () / norm;
					else
						non_linear_dif[k_system] = prev_non_linear_layers_solutions[k_system].Norm ();

					printf ("\t\tk_system: %i\tdif: %.3e\tnldiscr: %.3e\n", k_system, non_linear_dif[k_system], non_linear_discr[k_system]);
					sum_non_linear_dif += non_linear_dif[k_system];
					sum_non_linear_discr += non_linear_discr[k_system];

					//print_solutions (k_nl_iter);
				}
				if ((k_nl_iter > 0) && (sum_non_linear_dif < TCONVEC_PRECIS) && (sum_non_linear_discr / 3.0 < TCONVEC_DISCR_PRECIS))
					found = true;
			}

			// move time layers
			for (int k_system = 0; k_system < n_systems; k_system++)
			{
				for (int j = time_sampling - 1; j > 0; j--)
				{
					previous_time_layers_solutions[k_system][j] = previous_time_layers_solutions[k_system][j - 1];
				}
				previous_time_layers_solutions[k_system][0] = non_linear_layers_solutions[k_system];
			}
			print_solutions ();
		}

		delete[] non_linear_dif;
		delete[] prev_non_linear_layers_solutions;
		delete[] non_linear_layers_solutions;
		non_linear_layers_solutions = NULL;
		return true;
	}
}

double Task_convec::relax (int k_system, double w, const MathVector & prev_non_linear_layer_solution)
{
	double discr;
	MathVector * sol = new MathVector (eq_matrixes[k_system].Size ()); // matrix solution
	MathVector * mv = new MathVector (eq_matrixes[k_system].Size ());
	
	eq_matrixes[k_system].get_solution (sol);

	non_linear_layers_solutions[k_system].Linear_Combination (*sol, w, prev_non_linear_layer_solution, 1.0 - w);
	build_system (k_system);
	apply_boundary_conditions (k_system);
	eq_matrixes[k_system].mult_A_v (non_linear_layers_solutions[k_system], mv);
	mv->Substract (*(eq_matrixes[k_system].f));

	discr = mv->Norm ();
	
	double f_norm = eq_matrixes[k_system].f->Norm ();
	if (f_norm > 1e-13)
		discr /= f_norm;
	delete mv;
	delete sol;

	return discr;
}

double Task_convec::relax (int k_system, double w_left, double w_right, const MathVector & prev_non_linear_layer_solution)
{
	double w = w_right;
	double f_prev, f;
	double f_norm;
	int i_pos = 0;
	MathVector * sol = new MathVector (eq_matrixes[k_system].Size ()); // matrix solution
	eq_matrixes[k_system].get_solution (sol);
	MathVector * mv = new MathVector (eq_matrixes[k_system].Size ());

	int par_num_threads = 1;
	double * res_omp = NULL;
	if (par_num_threads > 1)
		res_omp = new double[(par_num_threads - 1) * eq_matrixes[k_system].Size ()];

	f_prev = 1e+10;
	bool found = false;
	for (int i = 0; i < 10 && w > w_left && !found; i++)
	{
		non_linear_layers_solutions[k_system].Linear_Combination (*sol, w, prev_non_linear_layer_solution, 1.0 - w);
		build_system (k_system);
		apply_boundary_conditions (k_system);
		eq_matrixes[k_system].mult_A_v (non_linear_layers_solutions[k_system], mv, res_omp, par_num_threads);
		mv->Substract (*(eq_matrixes[k_system].f));
		
		f = mv->Norm ();
		f_norm = eq_matrixes[k_system].f->Norm ();
		if (f_norm > 1e-13)
			f /= f_norm;

		if (f > f_prev)
		{
			i_pos = i - 1;
			found = true;
		}
		else
		{
			f_prev = f;
			w = w / 1.25;
		}
		//printf ("%i", i);
	}
	if (i_pos == -1)
		i_pos = 0;

	w = w_right / pow (1.25, i_pos);
	non_linear_layers_solutions[k_system].Linear_Combination (*sol, w, prev_non_linear_layer_solution, 1.0 - w);
	printf ("w:%.4lf", w);
	delete sol;
	if (par_num_threads > 1)
		delete[] res_omp;
	delete mv;

	return f_prev;
}

double Task_convec::relax (int k_system, const MathVector & prev_non_linear_layer_solution)
{
	double ms;
	clock_t t1, t2;

	double c = (1.0 + sqrt (5.0)) / 2.0; // coef
	MathVector * sol = new MathVector (eq_matrixes[k_system].Size ()); // matrix solution
	MathVector * mv = new MathVector (eq_matrixes[k_system].Size ());
	eq_matrixes[k_system].get_solution (sol);

	// w boundaries
	double w;
	double a = 0.0;
	double b = 1.0;
	// current boundaries
	double x1 = b - (b - a) / c;
	double x2 = a + (b - a) / c;

	double y1, y2; // function values at xs
	double f_norm;

	int par_num_threads = 1;
	double * res_omp = NULL;
	if (par_num_threads > 1)
		res_omp = new double[(par_num_threads - 1) * eq_matrixes[k_system].Size ()];

	// left
	non_linear_layers_solutions[k_system].Linear_Combination (*sol, x1, prev_non_linear_layer_solution, 1.0 - x1);
	build_system (k_system);
	apply_boundary_conditions (k_system);
	eq_matrixes[k_system].mult_A_v (non_linear_layers_solutions[k_system], mv, res_omp, par_num_threads);
	mv->Substract (*(eq_matrixes[k_system].f));
	y1 = mv->Norm ();
	f_norm = eq_matrixes[k_system].f->Norm ();
	if (f_norm > 1e-13)
		y1 /= f_norm;
	//printf ("%lf %.3e\n", x1, y1);
	// right
	non_linear_layers_solutions[k_system].Linear_Combination (*sol, x2, prev_non_linear_layer_solution, 1.0 - x2);
	build_system (k_system);
	apply_boundary_conditions (k_system);
	eq_matrixes[k_system].mult_A_v (non_linear_layers_solutions[k_system], mv, res_omp, par_num_threads);
	mv->Substract (*(eq_matrixes[k_system].f));
	y2 = mv->Norm ();
	f_norm = eq_matrixes[k_system].f->Norm ();
	if (f_norm > 1e-13)
		y2 /= f_norm;
	//printf ("%lf %.3e\n", x2, y2);

	for (int i = 0; (i < 100) && (fabs (b - a) > 2e-2); i++)
	{
		//printf ("%i", i);
		if (y1 >= y2)
		{
			a = x1;
			x1 = x2;
			y1 = y2;

			x2 = a + (b - a) / c;
			non_linear_layers_solutions[k_system].Linear_Combination (*sol, x2, prev_non_linear_layer_solution, 1.0 - x2);

			t1 = clock ();

			build_system (k_system);

			t2 = clock ();
			ms = (double)(t2 - t1);
			printf ("time spend building system: %lf s\n", ms / 1000);

			t1 = clock ();

			apply_boundary_conditions (k_system);

			t2 = clock ();
			ms = (double)(t2 - t1);
			printf ("time spend applying boundary cond: %lf s\n", ms / 1000);

			t1 = clock (); 

			eq_matrixes[k_system].mult_A_v (non_linear_layers_solutions[k_system], mv, res_omp, par_num_threads);

			t2 = clock ();
			ms = (double)(t2 - t1);
			printf ("time spend multiplying: %lf s\n", ms / 1000);

			mv->Substract (*(eq_matrixes[k_system].f));
			y2 = mv->Norm ();
			f_norm = eq_matrixes[k_system].f->Norm ();
			if (f_norm > 1e-13)
				y2 /= f_norm;
			//printf ("%lf %.3e\n", x2, y2);
		}
		else
		{
			b = x2;
			x2 = x1;
			y2 = y1;

			x1 = b - (b - a) / c;
			non_linear_layers_solutions[k_system].Linear_Combination (*sol, x1, prev_non_linear_layer_solution, 1.0 - x1);
			t1 = clock ();

			build_system (k_system);

			t2 = clock ();
			ms = (double)(t2 - t1);
			printf ("time spend building system: %lf s\n", ms / 1000);

			t1 = clock ();

			apply_boundary_conditions (k_system);

			t2 = clock ();
			ms = (double)(t2 - t1);
			printf ("time spend applying boundary cond: %lf s\n", ms / 1000);

			t1 = clock ();

			eq_matrixes[k_system].mult_A_v (non_linear_layers_solutions[k_system], mv, res_omp, par_num_threads);

			t2 = clock ();
			ms = (double)(t2 - t1);
			printf ("time spend multiplying: %lf s\n", ms / 1000);

			mv->Substract (*(eq_matrixes[k_system].f));
			y1 = mv->Norm ();
			f_norm = eq_matrixes[k_system].f->Norm ();
			if (f_norm > 1e-13)
				y1 /= f_norm;
			//printf ("%lf %.3e\n", x1, y1);
		}
	}


	// test purposes
	//if (K_NON_LINEAR_ITER > 9)
	//{
	//	char file_name[128];
	//	sprintf (file_name, "Test//s%i_t_%.5lf_nonlineariter%i.txt", k_system, time_layers[current_time_layer], K_NON_LINEAR_ITER);
	//	FILE * file = fopen (file_name, "w");
	//	double x;
	//	for (int i = 0; i < 101; i++)
	//	{
	//		printf ("%i ", i);
	//		x = i * 0.01;
	//		non_linear_layers_solutions[k_system].Linear_Combination (*sol, x, prev_non_linear_layer_solution, 1.0 - x);
	//		build_system (k_system);
	//		apply_boundary_conditions (k_system);
	//		eq_matrixes[k_system].mult_A_v (non_linear_layers_solutions[k_system], mv);
	//		mv->Substract (*(eq_matrixes[k_system].f));
	//		y1 = mv->Norm ();
	//		fprintf (file, "%lf %.5e\n", x, y1);
	//	}
	//	fclose (file);
	//	printf ("\n");
	//}


	w = (a + b) / 2.0;
	if (w > 0.99)
		w = 1.0;
	printf ("w:%.5lf ", w);
	non_linear_layers_solutions[k_system].Linear_Combination (*sol, w, prev_non_linear_layer_solution, 1.0 - w);

	delete sol;
	if (par_num_threads > 1)
		delete[] res_omp;
	delete mv;

	return y1;
}

void Task_convec::save_front_pos (char * file_name)
{
	double localMin, localMax;
	double point[2];
	int nodes[3];
	double v;
	double c1[2], c2[2];
	std::vector <std::pair<double, double>> fr_points;
	std::vector <std::pair<double, double>> fr_points_un;

	for (int k_element = 0, k_end = mesh_pointer->get_n_elements (); k_element < k_end; k_element++)
	{
		mesh_pointer->get_base_nodes (k_element, nodes);
		localMin = 1e+20;
		localMax = -1e+20;

		for (int k = 0; k < 3; k++)
		{
			mesh_pointer->get_node_coordinates (nodes[k], point);
			get_solution_in_point (0, point, &v);

			if (v > localMax)
				localMax = v;
			if (v < localMin)
				localMin = v;
		}

		if ((localMin - 1e-10 < T_phase_change) && (T_phase_change < localMax + 1e-10))
		{
			get_isoline_section (0, -1, k_element, T_phase_change, c1, c2);
			fr_points.push_back (std::pair<double, double> (c1[0], c1[1]));
			fr_points.push_back (std::pair<double, double>(c2[0], c2[1]));
		}
	}

	// sort points in X increasing order
	std::sort (fr_points.begin (), fr_points.end (), [](const std::pair<double, double> & a, const std::pair<double, double> & b)
	{
		return a.first < b.first;
	});

	{
		int i = 0;
		while (i < fr_points.size () - 1)
		{
			fr_points_un.push_back (fr_points[i]);
			if (fabs (fr_points[i].first - fr_points[i + 1].first) < 1e-7)
			{
				i += 2;
			}
			else
			{
				i += 1;
			}
		}
		if (i == (int)fr_points.size () - 1)
			fr_points_un.push_back (fr_points[(int)(fr_points.size ()) - 1]);
	}
	
	FILE * file = fopen (file_name, "w");
	for (int i = 0; i < fr_points_un.size (); i++)
	{
		fprintf (file, "%.15lf\t%.15lf\n", fr_points_un[i].first, fr_points_un[i].second);
	}
	fclose (file);
}

void Task_convec::get_conditions (int k_system, char * file_name)
{
	int t, f;
	switch (k_system)
	{
	case 1:
	{
		delete[] conditions[k_system];
		conditions[k_system] = new Condition[mesh_pointer->get_dimentionality () * 2 + 1];

		FILE * file = fopen (file_name, "r");
		for (int i = 0, i_end = mesh_pointer->get_dimentionality () * 2 + 1; i < i_end; i++)
		{
			fscanf (file, "%i %i", &t, &f);
			conditions[k_system][i].set_Type (t);
			conditions[k_system][i].set_function (f);
		}
		fclose (file);

		break;
	}
	default:
	{
		FILE * file = fopen (file_name, "r");
		for (int i = 0, i_end = mesh_pointer->get_dimentionality () * 2; i < i_end; i++)
		{
			fscanf (file, "%i %i", &t, &f);
			conditions[k_system][i].set_Type (t);
			conditions[k_system][i].set_function (f);
		}
		fclose (file);
	}
	}
}

double Task_convec::get_function_value (int k_system, int k_function, double * coordinates, int area)
{
	double r = 0.0;
	switch (k_function)
	{
	case 0:
		// f function
		r = function_f (k_system, coordinates, area);
		break;
	}
	return r;
}

double Task_convec::function_f (int k_system, double * coordinates, int area)
{
	double r = 0.0;
	return r;
}

double Task_convec::function_starting_condition (int k_system, double * coordinates, int area)
{
	double r = 0.0;
	switch (k_system)
	{
	case 0:
	{
		// switch area
		switch (area)
		{
		case 1:
			r = T0_fluid;
			break;
		case 2:
			r = T0_solid;
			break;
		}
		break;
	}
	case 1:
	case 2:
		r = 0.0;
		break;
	}
	return r;
}

double Task_convec::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	double r = 0.0;

	// boundary == function number
	switch (k_system)
	{
	case 0:
	{
		switch (boundary)
		{
		case 2:
			r = T_source;
			break;
		case 3:
			r = T0_solid;
			break;
		case BOUND_FUNC_HEAT_SOURCE:
			r = T_source;
			break;
		}
		break;
	}
	case 1:
	{
		double step = 1e-2;
		double point[] = { 0.0, 0.0 };
		double c0[2], cn[2];
		mesh_pointer->get_0_boundaries (c0);
		mesh_pointer->get_N_boundaries (cn);
		switch (boundary)
		{
		case 1:
		{
			point[0] = c0[0] + step;
			point[1] = coordinates[1];
			break;
		}
		case 2:
		{
			point[0] = cn[0] - step;
			point[1] = coordinates[1];
			break;
		}
		case 3:
		{
			point[0] = coordinates[0];
			point[1] = c0[1] + step;
			break;
		}
		case 4:
		{
			point[0] = coordinates[0];
			point[1] = cn[1] - step;
			break;
		}
		case 5:
		{
			// anti-symmetry
			return 0.0;
			break;
		}
		}

		double value = 0.0;
		get_solution_in_point (2, point, &value);
		r = -2.0 * value / pow (step, 2.0);
		//r = 0.0;
		break;
	}
	case 2:
	{
		r = 0.0;
		break;
	}
	}
	return r;
}

void Task_convec::get_function_coefficients (int k_system, int function, MathVector * coef)
{
	double coord[2];
	int area;
	for (int i = 0; i < mesh_pointer->get_n_nodes (); i++)
	{
		mesh_pointer->get_node_coordinates (i, coord);
		area = mesh_pointer->get_area (coord);
		coef->setElem (i, get_function_value (k_system, function, coord, area));
	}
}

void Task_convec::get_density_coefficients (int k_element, int function, MathVector * coef)
{
	// for water
	
	{
		int iF;
		int n_functions = mesh_pointer->get_amount_non_zero_functions (0);
		double T0 = 4.029325;
		double q = 1.894816;
		double T, T_NON_SCALED;
		double DENS_SCALE = 999.972;
		double dens;

		for (int i = 0; i < n_functions; i++)
		{
			iF = get_function_global_number (k_element, i);
			T = previous_time_layers_solutions[0][0].getElem (iF);
			T_NON_SCALED = T_min + T * THETA;
			dens = pow (abs (T_NON_SCALED - T0), q);
			dens /= THETA;
			//dens /= DENS_SCALE; // ????
			coef->setElem (i, dens);
		}
	}
}

double Task_convec::get_density_T_IP_scaled ()
{
	// for water
	double T0 = 4.029325;
	double T0_SCALED = (T0 - T_min) / THETA;
	return T0_SCALED;
}

void Task_convec::build_system (int k_system)
{
	eq_matrixes[k_system].Clear ();
	int n_functions = mesh_pointer->get_amount_non_zero_functions (0);
	int iF, jF, kF;

	switch (k_system)
	{
	case 0:
	{
		Matrix * G;
		Matrix * M;
		MathVector * DDx;
		MathVector * DDy;
		MathVector * Time_coef;
		MathVector * T_prev;
		MathVector * Psi;

		G = new Matrix (n_functions, n_functions);
		M = new Matrix (n_functions, n_functions);
		DDx = new MathVector (n_functions * n_functions * n_functions);
		DDy = new MathVector (n_functions * n_functions * n_functions);
		T_prev = new MathVector (n_functions);
		Psi = new MathVector (n_functions);

		Time_coef = new MathVector (time_sampling);
		get_time_coefficients (Time_coef);

		double lambda;
		double melt_coef;
		double Pr;

		for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
		{
			// get G matrix
			mesh_pointer->get_G_local_matrix (k_element, G);
			// get M matrix
			mesh_pointer->get_M_local_matrix (k_element, M);
			// build Psi
			for (int i = 0; i < n_functions; i++)
			{
				Psi->setElem (i, non_linear_layers_solutions[2].getElem (get_function_global_number (k_element, i)));
			}
			// get DDs
			mesh.elements[k_element]->get_DDF (1, 0, DDx);
			mesh.elements[k_element]->get_DDF (0, 1, DDy);

			// get element's lambda
			lambda = Lambda (k_element);

			// melt_coef
			if (MELT_SETUP == 1)
			{
				melt_coef = 1.0;
				
				double min_T = 1e+20, max_T = -1e+20;
				double T_val;

				for (int i = 0; i < n_functions; i++)
				{
					T_val = previous_time_layers_solutions[0][0].getElem (get_function_global_number (k_element, i));
					if (T_val > max_T)
						max_T = T_val;
					if (T_val < min_T)
						min_T = T_val;
				}
				if ((min_T - SOLID_PRECIS < T_phase_change) && (T_phase_change < max_T + SOLID_PRECIS))
				{
					melt_coef += 1.0 / Ste;
				}
			}
			else
				melt_coef = 1.0;

			Pr = get_Pr (k_element);

			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (k_element, i);
				if ((iF != -1) && (!ffb[k_system][iF]))
				{
					for (int j = 0; j < n_functions; j++)
					{
						jF = get_function_global_number (k_element, j);

						if (jF != -1)
						{
							// add G
							eq_matrixes[k_system].add_to_entry (iF, jF, G->Elem (i, j) * lambda / (lambda_fluid * Pr));

							// add M
							eq_matrixes[k_system].add_to_entry (iF, jF, melt_coef * Time_coef->getElem (0) * M->Elem (i, j));

							// add and substract partial derivatives
							for (int k = 0; k < n_functions; k++)
							{
								eq_matrixes[k_system].add_to_entry (iF, jF, melt_coef * DDx->getElem (i * n_functions * n_functions + j * n_functions + k) * Psi->getElem (k));
								eq_matrixes[k_system].add_to_entry (iF, jF, melt_coef * (-DDy->getElem (i * n_functions * n_functions + j * n_functions + k)) * Psi->getElem (k));
							}
						}
					}
				}
			}

			// add prev time layers data
			{
				// go by amount of time layers that count
				for (int t = 1; t < time_sampling; t++)
				{
					// get previous time layers solutions
					for (int i = 0; i < n_functions; i++)
					{
						iF = get_function_global_number (k_element, i);
						if (iF != -1)
						{
							T_prev->setElem (i, previous_time_layers_solutions[k_system][t - 1].getElem (iF));
						}
					}
					// go by element's functions
					for (int i = 0; i < n_functions; i++)
					{
						iF = get_function_global_number (k_element, i);
						if ((iF != -1) && (!ffb[k_system][iF]))
						{
							for (int j = 0; j < n_functions; j++)
							{
								eq_matrixes[k_system].add_to_f_entry (iF, melt_coef * T_prev->getElem (j) * Time_coef->getElem (t) * M->Elem (i, j));
							}
						}
					}
				}
			}
		}

		delete G;
		delete M;
		delete DDx;
		delete DDy;
		delete Time_coef;
		delete T_prev;
		delete Psi;
		break;
	}
	case 1:
	{
		Matrix * G;
		Matrix * M;
		Matrix * D;
		MathVector * DDx;
		MathVector * DDy;
		MathVector * Time_coef;
		MathVector * W_prev;
		MathVector * T;
		MathVector * Psi;

		G = new Matrix (n_functions, n_functions);
		M = new Matrix (n_functions, n_functions);
		D = new Matrix (n_functions, n_functions);
		DDx = new MathVector (n_functions * n_functions * n_functions);
		DDy = new MathVector (n_functions * n_functions * n_functions);
		W_prev = new MathVector (n_functions);
		Psi = new MathVector (n_functions);
		T = new MathVector (n_functions);

		Time_coef = new MathVector (time_sampling);
		get_time_coefficients (Time_coef);

		double Gr;

		for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
		{
			// get G matrix
			mesh_pointer->get_G_local_matrix (k_element, G);
			// get M matrix
			mesh_pointer->get_M_local_matrix (k_element, M);
			// get bf * dbf/dx
			mesh.elements[k_element]->get_D (0, D);

			// build Psi
			for (int i = 0; i < n_functions; i++)
			{
				Psi->setElem (i, non_linear_layers_solutions[2].getElem (get_function_global_number (k_element, i)));
			}
			// build T
			if (USE_DENSITY_FUNCTION == 0)
			{
				for (int i = 0; i < n_functions; i++)
				{
					T->setElem (i, non_linear_layers_solutions[0].getElem (get_function_global_number (k_element, i)));
				}
			}
			// or density
			if (USE_DENSITY_FUNCTION == 1)
				get_density_coefficients (k_element, density_function, T);

			// get DDs
			mesh.elements[k_element]->get_DDF (1, 0, DDx);
			mesh.elements[k_element]->get_DDF (0, 1, DDy);

			Gr = get_Gr (k_element);

			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (k_element, i);
				if ((iF != -1) && (!ffb[k_system][iF]) && (state[iF] == 0))
				{
					for (int j = 0; j < n_functions; j++)
					{
						jF = get_function_global_number (k_element, j);

						if (jF != -1)
						{
							// add G
							eq_matrixes[k_system].add_to_entry (iF, jF, G->Elem (i, j));

							// add M
							eq_matrixes[k_system].add_to_entry (iF, jF, Time_coef->getElem (0) * M->Elem (i, j));

							// add D (dT/dx) into f
							eq_matrixes[k_system].add_to_f_entry (iF, D->Elem (i, j) * T->getElem (j) * Gr);

							// add and substract partial derivatives
							for (int k = 0; k < n_functions; k++)
							{
								eq_matrixes[k_system].add_to_entry (iF, jF, DDx->getElem (i * n_functions * n_functions + j * n_functions + k) * Psi->getElem (k));
								eq_matrixes[k_system].add_to_entry (iF, jF, -DDy->getElem (i * n_functions * n_functions + j * n_functions + k) * Psi->getElem (k));
							}
						}
					}
				}
			}

			// add prev time layers data
			{
				// go by amount of time layers that count
				for (int t = 1; t < time_sampling; t++)
				{
					// get previous time layers solutions
					for (int i = 0; i < n_functions; i++)
					{
						iF = get_function_global_number (k_element, i);
						if (iF != -1)
						{
							W_prev->setElem (i, previous_time_layers_solutions[k_system][t - 1].getElem (iF));
						}
					}
					// go by element's functions
					for (int i = 0; i < n_functions; i++)
					{
						iF = get_function_global_number (k_element, i);
						if ((iF != -1) && (!ffb[k_system][iF]) && (state[iF] == 0))
						{
							for (int j = 0; j < n_functions; j++)
							{
								eq_matrixes[k_system].add_to_f_entry (iF, W_prev->getElem (j) * Time_coef->getElem (t) * M->Elem (i, j));
							}
						}
					}
				}
			}
		}

		delete G;
		delete M;
		delete D;
		delete DDx;
		delete DDy;
		delete Time_coef;
		delete W_prev;
		delete Psi;
		delete T;
		break;
	}
	case 2:
	{
		Matrix * G;
		Matrix * M;
		MathVector * W;

		G = new Matrix (n_functions, n_functions);
		M = new Matrix (n_functions, n_functions);
		W = new MathVector (n_functions);

		for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
		{
			// get G matrix
			mesh_pointer->get_G_local_matrix (k_element, G);
			// get M matrix
			mesh_pointer->get_M_local_matrix (k_element, M);
			// get W vector
			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (k_element, i);
				if (iF != -1)
				{
					W->setElem (i, non_linear_layers_solutions[1].getElem (iF));
				}
			}

			// go by element's functions
			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (k_element, i);
				if ((iF != -1) && (!ffb[k_system][iF]) && (state[iF] == 0))
				{
					for (int j = 0; j < n_functions; j++)
					{
						jF = get_function_global_number (k_element, j);

						if (jF != -1)
						{
							// add G into respective places (by def_nodes)
							eq_matrixes[k_system].add_to_entry (iF, jF, G->Elem (i, j));

							// add W * M into f
							eq_matrixes[k_system].add_to_f_entry (iF, W->getElem (j) * M->Elem (i, j));
						}
					}
				}
			}

		}

		delete G;
		delete M;
		delete W;
		break;
	}
	}
}

void Task_convec::print_solutions (int knliter)
{
	// picture for every non-linear solution 
	std::vector<double> local_time_stamps;
	bool save_solution = false;

	if (use_time_mapping)
	{
		// save solution in current time
		save_solution = true;
		local_time_stamps.push_back (time_layers[current_time_layer]);
	}
	else
	{
		// if there is some scheme
		// save which times fall between last one and current
		for (size_t i = 0, i_end = time_stamps.size (); i < i_end; i++)
		{
			if (time_layers[current_time_layer - 1] + 1e-7 < time_stamps[i] && time_stamps[i] < time_layers[current_time_layer] + 1e-7)
			{
				local_time_stamps.push_back (time_stamps[i]);
			}
		}
		if (local_time_stamps.size () != 0)
		{
			save_solution = true;
		}
	}
	if (save_solution)
	{
		FILE * file;
		MathVector * c = new MathVector (time_sampling);

		char name[64];
		for (int k_system = 0; k_system < n_systems; k_system++)
		{
			// fprint solution
			for (size_t i = 0, i_end = local_time_stamps.size (); i < i_end; i++)
			{
				get_time_approx (local_time_stamps[i], c);
				sprintf (name, "Result//TC//s%i_t_%.5lf_nliter_%i.txt", k_system, local_time_stamps[i], knliter);
				file = fopen (name, "w");
				non_linear_layers_solutions[k_system].FPrint (file);
				fclose (file);
			}
		}
		delete c;
	}

}

void Task_convec::print_solutions ()
{
	std::vector<double> local_time_stamps;
	bool save_solution = false;

	if (use_time_mapping)
	{
		// save solution in current time
		save_solution = true;
		local_time_stamps.push_back (time_layers[current_time_layer]);
	}
	else
	{
		// if there is some scheme
		// save which times fall between last one and current
		for (size_t i = 0, i_end = time_stamps.size (); i < i_end; i++)
		{
			if (time_layers[current_time_layer - 1] + 1e-7 < time_stamps[i] && time_stamps[i] < time_layers[current_time_layer] + 1e-7)
			{
				local_time_stamps.push_back (time_stamps[i]);
			}
		}
		if (local_time_stamps.size () != 0)
		{
			save_solution = true;
		}
	}
	if (save_solution)
	{
		FILE * file;
		MathVector * c = new MathVector (time_sampling);

		char name[64];
		for (int k_system = 0; k_system < n_systems; k_system++)
		{
			// fprint solution
			//if (current_time_layer % 100 == 0)
			{
				for (size_t i = 0, i_end = local_time_stamps.size (); i < i_end; i++)
				{
					get_time_approx (local_time_stamps[i], c);
					printf ("\tsaving solutions into files, time layer:\t%.13lf\n", local_time_stamps[i]);
					sprintf (name, "Result//TC//s%i_t_%.5lf.txt", k_system, local_time_stamps[i]);
					file = fopen (name, "w");
					double val;
					for (int j = 0, j_end = previous_time_layers_solutions[k_system][0].getSize (); j < j_end; j++)
					{
						val = 0;
						for (int t = 0; t < time_sampling; t++)
						{
							val += c->getElem (t) * previous_time_layers_solutions[k_system][t].getElem (j);
						}
						fprintf (file, "%.13lf\n", val);
					}
					fclose (file);
				}
			}

			char pic_name[128];
			wchar_t wtemp[128];

			// draw solution
			for (size_t i = 0, i_end = local_time_stamps.size (); i < i_end; i++)
			{
				switch (k_system)
				{
				case 0:
					//set_symmetry (1);
					painter->set_min_max (T0_solid, T_source);
					painter->draw_contour_lines (k_system, 10);
					painter->set_legend_full_iso ();

					if (MELT_SETUP == 1)
						painter->add_isoline (T_phase_change);
					if (USE_DENSITY_FUNCTION)
						painter->add_isoline (get_density_T_IP_scaled());
					break;
				case 1:
					painter->draw_contour_lines (k_system, 20);
					//set_symmetry (2);
					break;
				case 2:
					painter->draw_contour_lines (k_system, 20);
					//set_symmetry (2);
					break;
				}
				painter->set_max_resolution (900);
				painter->set_axis_names ("X", "Y");
				//painter->add_area_dividers ();

				painter->set_iso_line_width (1.0);
				painter->set_axes_scale (10, 10);

				painter->draw_field (k_system, COLOR_SCALE_RED_YELLOW);

				sprintf (name, "Pictures//TC//s%i_t_%.5lf.png", k_system, local_time_stamps[i]);

				mbstowcs (wtemp, name, strlen (name) + 1);
				std::wstring w_name = wtemp;

				painter->draw_to_file (w_name);
				painter->reset ();
				set_symmetry (0);
			}
		}
		delete c;

		// TEST slices
		//for (int k_system = 0; k_system < n_systems; k_system++)
		//{
		//	for (size_t i = 0, i_end = local_time_stamps.size (); i < i_end; i++)
		//	{
		//		{
		//			char name[64];
		//			sprintf (name, "Result//TC//svert%i_t_%.5lf.txt", k_system, local_time_stamps[i]);
		//			save_slice (k_system, "Source Files//TC//slice_vert.txt", name, false);
		//		}
		//		{
		//			char name[64];
		//			sprintf (name, "Result//TC//shor%i_t_%.5lf.txt", k_system, local_time_stamps[i]);
		//			save_slice (k_system, "Source Files//TC//slice_hor.txt", name, false);
		//		}
		//	}
		//}
	}
}

void Task_convec::get_Theta_edge (int k_system, int k_element, int n1, int n2, int * functions, int boundary, MathVector * Theta)
{
	// only 0 values and only for k_system = 0;
	for (int i = 0; i < Theta->getSize (); i++)
		Theta->setElem (i, 0.0);
}

double Task_convec::Lambda (int k_element)
{
	double l = 0.0;

	// non-melt
	if (MELT_SETUP == 0)
	{
		int area = mesh_pointer->get_area (k_element);
		switch (area)
		{
		case 1:
			l = lambda_fluid;
			break;
		case 2:
			l = lambda_solid;
			break;
		}
	}

	if (MELT_SETUP == 1)
	{ 
		// average on the element
		int n_functions = mesh_pointer->get_amount_non_zero_functions (k_element);
		int iF;

		for (int i = 0; i < n_functions; i++)
		{
			iF = get_function_global_number (k_element, i);
			if (iF != -1)
			{
				if (previous_time_layers_solutions[0][0].getElem (iF) > T_phase_change - SOLID_PRECIS)
					l += lambda_fluid;
				else
					l += lambda_solid;
			}
		}

		l /= n_functions;
	}

	return l;
}

void Task_convec::apply_boundary_conditions (int k_system)
{
	apply_second_boundary_conditions (k_system);
	apply_first_boundary_conditions (k_system);
}

void Task_convec::apply_first_boundary_conditions (int k_system)
{
	// apply first boundary conditions
	{
		double coord[2];

		// set diagonal element if ffb = true
		for (int iF = 0; iF < N_functions; iF++)
		{
			if (ffb[k_system][iF])
			{
				mesh_pointer->get_node_coordinates (iF, coord);
				eq_matrixes[k_system].set_entry (iF, iF, 1.0);
				eq_matrixes[k_system].set_f_entry (iF, function_FCondition (k_system, coord, 0, ffb[k_system][iF]));
			}
		}
	}

	// add extra stuff
	switch (k_system)
	{
	case 0:
	{
		if (LOCALIZED_BOTTOM_HEAT_SOURCE)
		{
			int boundary;
			int iF;
			MathVector * FCondition;
			int * functions;
			int f_amount;
			int k_element;
			int n1, n2;
			double coord[2];
			int nodes_within;

			bool * replaced = new bool[N_functions];
			for (int i = 0; i < N_functions; i++)
			{
				replaced[i] = false;
			}

			// first condition by edges
			for (int i_edge = 0, i_end = mesh_pointer->edges->get_n_entries (); i_edge < i_end; i_edge++)
			{
				// get edges's nodes
				mesh_pointer->edges->get_edge_nodes (i_edge, &n1, &n2);
				// get boundary number for those nodes
				boundary = mesh_pointer->function_boundary_edge (n1, n2);
				if (boundary == 2) // bottom boundary
				{
					// check that nodes are within the length
					nodes_within = 0;
					mesh_pointer->get_node_coordinates (n1, coord);
					if (coord[0] < heat_source_length + 1e-13)
						nodes_within++;
					mesh_pointer->get_node_coordinates (n2, coord);
					if (coord[0] < heat_source_length + 1e-13)
						nodes_within++;

					if (nodes_within == 2)
					{
						// find element with that edge
						k_element = mesh_pointer->belonging_element (n1, n2);
						// get amount of values to put into f (from element)
						f_amount = mesh_pointer->get_amount_second_condition (k_element);
						FCondition = new MathVector (f_amount);

						functions = new int[f_amount];
						// get functions that are not 0 on the edge
						mesh_pointer->get_edge_functions (k_element, n1, n2, functions);

						// get first condition from element
						get_FCondition_edge (k_system, k_element, n1, n2, functions, BOUND_FUNC_HEAT_SOURCE, FCondition);

						// add it to f in places by functions
						for (int i = 0; i < f_amount; i++)
						{
							iF = get_function_global_number (k_element, functions[i]);
							if (iF != -1)
							{
								if (!replaced[iF])
								{
									replaced[iF] = true;
									eq_matrixes[k_system].clear_row (iF); // replace row
									eq_matrixes[k_system].set_entry (iF, iF, 1.0);
									eq_matrixes[k_system].set_f_entry (iF, FCondition->getElem (i));
								}
							}
						}

						delete FCondition;
						delete[] functions;
					}
				}
			}
			delete[] replaced;
		}
		break;
	}
	case 1:
	{
		double coord[2];

		for (int iF = 0; iF < N_functions; iF++)
		{
			mesh_pointer->get_node_coordinates (iF, coord);
			if (state[iF] == 1)
			{
				eq_matrixes[k_system].set_entry (iF, iF, 1.0);
				eq_matrixes[k_system].set_f_entry (iF, 0.0);
			}

			if (state[iF] == 2)
			{
				eq_matrixes[k_system].set_entry (iF, iF, function_FCondition (1, coord, 0, conditions[k_system][4].Type ()));
				eq_matrixes[k_system].set_f_entry (iF, 0.0);
			}
		}

	}
	case 2:
	{
		double coord[2];

		for (int iF = 0; iF < N_functions; iF++)
		{
			mesh_pointer->get_node_coordinates (iF, coord);
			if (state[iF] > 0)
			{
				eq_matrixes[k_system].set_entry (iF, iF, 1.0);
				eq_matrixes[k_system].set_f_entry (iF, 0.0);
			}
		}
		break;

	}
	}
}

void Task_convec::apply_second_boundary_conditions (int k_system)
{
}

void Task_convec::set_ffb_flags ()
{
	ffb = new int *[n_systems];
	for (int i = 0; i < n_systems; i++)
	{
		ffb[i] = new int[N_functions];
		for (int k = 0; k < N_functions; k++)
		{
			ffb[i][k] = 0;
		}
	} 

	int boundary;
	int iF;
	MathVector * FCondition;
	int * functions;
	int f_amount;
	int k_element;
	int n1, n2;
	
	for (int k_system = 0; k_system < n_systems; k_system++)
	{
		// first condition by edges
		for (int i_edge = 0, i_end = mesh_pointer->edges->get_n_entries (); i_edge < i_end; i_edge++)
		{
			// get edges's nodes
			mesh_pointer->edges->get_edge_nodes (i_edge, &n1, &n2);
			// get boundary number for those nodes
			boundary = mesh_pointer->function_boundary_edge (n1, n2);
			if (boundary != -1)
			{
				if (conditions[k_system][boundary].Type () == 1)
				{
					// find element with that edge
					k_element = mesh_pointer->belonging_element (n1, n2);
					// get amount of values to put into f (from element)
					f_amount = mesh_pointer->get_amount_second_condition (k_element);

					FCondition = new MathVector (f_amount);
					functions = new int[f_amount];
					// get functions that are not 0 on the edge
					mesh_pointer->get_edge_functions (k_element, n1, n2, functions);

					for (int i = 0; i < f_amount; i++)
					{
						iF = get_function_global_number (k_element, functions[i]);
						ffb[k_system][iF] = conditions[k_system][boundary].Function();
					}

					delete FCondition;
					delete[] functions;
				}
			}
		}
	}
}

void Task_convec::reset_states ()
{
	int area;
	int zone;
	int iF;
	int n_functions;

	for (int iF = 0; iF < N_functions; iF++)
		state[iF] = 0;

	for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
	{
		zone = 0; // fluid
		area = mesh_pointer->get_area (k_element);
		n_functions = mesh_pointer->get_amount_non_zero_functions (k_element);

		// for non-melt
		if (MELT_SETUP == 0)
		{
			if (area == 2)
				zone = 1; // solid

			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (k_element, i);
				state[iF] = zone;
			}
		}
		if (MELT_SETUP == 1)
		{
			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (k_element, i);
				if (iF != -1)
				{
					if ((T_phase_change - SOLID_PRECIS < previous_time_layers_solutions[0][0].getElem (iF)) &&
						(previous_time_layers_solutions[0][0].getElem (iF) < T_phase_change + SOLID_PRECIS))
						zone = 2; // mush

					if (previous_time_layers_solutions[0][0].getElem (iF) < T_phase_change - SOLID_PRECIS)
						zone = 1; // solid

					if (state[iF] < zone)
						state[iF] = zone;
				}
			}
		}
	}
}

double Task_convec::get_Pr (int k_element)
{
	if (PRGRconst)
		return PR;
	double Pr = 0.0;
	double v, c[1];
	int i, i_end;

	for (i = 0, i_end = mesh_pointer->get_amount_non_zero_functions (k_element); i < i_end; i++)
	{
		v = 0.0;
		c[0] = previous_time_layers_solutions[0][0].getElem (mesh_pointer->get_function_global_number(k_element, i))
			* THETA + T_min;
		spline_Pr->get_solution_in_point (c, &v);
		Pr += v;
	}

	Pr /= (double)i_end;
	return Pr;
}

double Task_convec::get_Gr (int k_element)
{
	if (PRGRconst)
		return GR;

	double Gr = 0.0;
	double v, c[1];
	int i, i_end;

	for (i = 0, i_end = mesh_pointer->get_amount_non_zero_functions (k_element); i < i_end; i++)
	{
		v = 0.0;
		c[0] = previous_time_layers_solutions[0][0].getElem (mesh_pointer->get_function_global_number (k_element, i))
			* THETA + T_min;
		spline_Gr->get_solution_in_point (c, &v);
		Gr += v;
	}

	Gr /= (double)i_end;
	return Gr;
}

void Task_convec::set_n_systems ()
{
	n_systems = 3;
	eq_matrixes = new compressed_matrix[n_systems];
	printf ("systems settled:\t%i\n", n_systems);
	fprintf (log_file, "systems settled:\t%i\n", n_systems);
}

Task_convec::Task_convec ()
{
	ffb = NULL;
	spline_Pr = NULL;
	spline_Gr = NULL;
	state = NULL;
}

Task_convec::Task_convec (char * file_param)
{
	FILE * file = fopen (file_param, "r");
	fscanf (file, "%i", &PRGRconst);
	fscanf (file, "%lf", &PR); // PR
	fscanf (file, "%lf", &GR); // GR
	fscanf (file, "%lf", &T_min); // minimum temp in the system
	fscanf (file, "%lf", &THETA); // temp scale

	fscanf (file, "%lf", &T0_solid); // T0_solid
	fscanf (file, "%lf", &T0_fluid); // T0_fluid
	fscanf (file, "%lf", &T_source); // T_source, temperature of the heat source

	fscanf (file, "%i", &LOCALIZED_BOTTOM_HEAT_SOURCE); // 1 for localized bottom heat source
	if (LOCALIZED_BOTTOM_HEAT_SOURCE)
		fscanf (file, "%lf", &heat_source_length);

	fscanf (file, "%lf", &lambda_fluid); 
	fscanf (file, "%lf", &lambda_solid);

	fscanf (file, "%i", &MELT_SETUP); // 1 for including melting processes
	if (MELT_SETUP)
	{
		fscanf (file, "%lf", &Ste);
		fscanf (file, "%lf", &T_phase_change);
	}
	fscanf (file, "%i", &USE_DENSITY_FUNCTION); // 1 for non-linear density
	if (USE_DENSITY_FUNCTION)
	{
		fscanf (file, "%i", &density_function); 
	}

	fclose (file);

	ffb = NULL;
	spline_Pr = NULL;
	spline_Gr = NULL;
	state = NULL;
}

Task_convec::~Task_convec ()
{
	if (ffb != NULL)
	{
		for (int i = 0; i < n_systems; i++)
			delete[] ffb[i];
		delete[] ffb;
	}

	if (spline_Pr != NULL)
		delete spline_Pr;
	if (spline_Gr != NULL)
		delete spline_Gr;
	if (state != NULL)
		delete[] state;
}
