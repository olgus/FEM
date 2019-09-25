#include "Find_theta.h"
#include "Free_meshing.h"

#include <omp.h>

void Task_Find_Theta::solve (char * file_front_name, Painter * painter)
{
	double theta_av = 15.0;
	double GR_no_theta, Ste_no_theta, PC_no_theta, phys_T0_F_no_Theta;

	GR_no_theta = 570.491136260751;
	Ste_no_theta = 0.00492537313432836;
	PC_no_theta = 22.0 - 21.0;
	phys_T0_F_no_Theta = 23.0 - 21.0;

	double GR, Ste, PC, phys_T0_F;
	double theta = theta_av;
	double theta_prev;
	double L; // new direction
	double alpha; // regularization
	double r_prev, r, s, i;
	double beta; // relaxation?

	// get front position
	int n_receivers;
	FILE * file_front = fopen (file_front_name, "r");
	std::vector<receiver_data> front_pos;
	{
		double x, y;
		while (!feof (file_front))
		{

			fscanf (file_front, "%lf %lf", &x, &y);
			front_pos.push_back ({ x, y });
		}
		front_pos.pop_back ();
		n_receivers = (int)front_pos.size ();
		printf ("Got %i receivers\n", n_receivers);
	}
	fclose (file_front);

	MathVector * R = new MathVector (n_receivers);
	MathVector * F = new MathVector (n_receivers);
	//// same mesh for the task
	//Free_mesh mesh;
	//mesh.build_Mesh_wo_renumerate ("Source Files//Melt_SF//Mesh_free.txt", "Source Files//Melt_SF//density.txt", "Source Files//Melt_SF//order.txt", "MeshData//Nodes.txt", "MeshData//Triangles.txt");

	FILE * file_RT = fopen ("Result//Melt_SF//RT//rt.txt", "a+");
	fprintf (file_RT, "iter r i s alpha theta GR Ste PC phys_T0_F\n");
	fclose (file_RT);

	double h = 0.005;
	r = 1.0;
	r_prev = 0.0;
	alpha = 0.0;

	Task_Melt_SF<Triangular_Mesh> * melt;
	Task_Melt_SF<Triangular_Mesh> * melt_d;

	// 0 iter
	bool found = false;
	for (int iter = 0; iter < FIND_THETA_MAX_ITER && !found; iter++)
	{
		beta = 1.0;
		theta_prev = theta;
		while (beta > BETA_PRECIS)
		{
			// calc new parameters
			GR = GR_no_theta * theta;
			Ste = Ste_no_theta * theta;
			PC = PC_no_theta / theta;
			phys_T0_F = phys_T0_F_no_Theta / theta;

			printf ("Parameters: %lf %lf %lf %lf\n", GR, Ste, PC, phys_T0_F);
#pragma omp parallel sections
			{
				// solve task with current theta
#pragma omp section
				{
					melt = new Task_Melt_SF<Triangular_Mesh> ("Source Files//Melt_SF//Parameters.txt");
					char * bound[] = { "Source Files//Melt_SF//Conditions_T.txt",
						"Source Files//Melt_SF//Conditions_P.txt",
						"Source Files//Melt_SF//Conditions_W.txt" };
					melt->prepare ("MeshData//Nodes.txt", "MeshData//Triangles.txt", bound, "Source Files//Melt_SF//time_mapping.txt", "Source Files//Melt_SF//time_stamps.txt");
					melt->get_scales ("Source Files//Melt_SF//Scales.txt");
					melt->reset_parameters (theta, GR, Ste, PC, phys_T0_F);
					int param[][5] = { { SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 1, 4, SOLVER_MKL_NO },
					{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 1, 1, SOLVER_MKL_NO },
					{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 1, 4, SOLVER_MKL_NO } };
					melt->solve_task (param);
				}

				// solve task for derivatives
#pragma omp section
				{
					melt_d = new  Task_Melt_SF<Triangular_Mesh> ("Source Files//Melt_SF//Parameters.txt");
					double GR_der = GR_no_theta * (theta + h * theta);
					double Ste_der = Ste_no_theta * (theta + h * theta);
					double PC_der = PC_no_theta / (theta + h * theta);
					double phys_T0_F_der = phys_T0_F_no_Theta / (theta + h * theta);
					char * bound[] = { "Source Files//Melt_SF//Conditions_T.txt",
						"Source Files//Melt_SF//Conditions_P.txt",
						"Source Files//Melt_SF//Conditions_W.txt" };
					melt_d->prepare ("MeshData//Nodes.txt", "MeshData//Triangles.txt", bound, "Source Files//Melt_SF//time_mapping.txt", "Source Files//Melt_SF//time_stamps.txt");
					melt_d->get_scales ("Source Files//Melt_SF//Scales.txt");
					melt_d->reset_parameters (theta + h * theta, GR_der, Ste_der, PC_der, phys_T0_F_der);
					int param[][5] = { { SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 1, 4, SOLVER_MKL_NO },
					{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 1, 1, SOLVER_MKL_NO },
					{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 1, 4, SOLVER_MKL_NO } };
					melt_d->solve_task (param);
				}
			}
			// save pics and full front position for future references
			{
				// pic
				char name[128];
				wchar_t wtemp[128];

				sprintf (name, "Pictures//Melt_SF//RT//iter_%i_theta_%lf.png", iter, theta);
				mbstowcs (wtemp, name, strlen (name) + 1);
				std::wstring w_name = wtemp;

				Task_pointer * tp = melt;
				painter->set_max_resolution (800);
				painter->set_task (tp);
				painter->set_min_max (0.0, 1.0);
				painter->draw_field (0, COLOR_SCALE_RED_YELLOW);
				painter->draw_contour_lines (0, 20);
				painter->set_iso_line_width (1.0);
				painter->add_isoline (PC);
				painter->draw_to_file (w_name);
				painter->reset ();

				// front
				char front_name[128];
				sprintf (front_name, "Result//Melt_SF//RT//iter_%i_theta_%lf.txt", iter, theta);
				FILE * file_front_loc = fopen (front_name, "w");
				double fp;
				for (int i_r = 0; i_r < n_receivers; i_r++)
				{
					fp = melt->get_front_position (front_pos[i_r].x);
					fprintf (file_front_loc, "%lf\n", fp);
				}
				fclose (file_front_loc);
			}
			{
				// pic
				char name[128];
				wchar_t wtemp[128];

				sprintf (name, "Pictures//Melt_SF//RT//iter_%i_theta_%lf.png", iter, theta + h * theta);
				mbstowcs (wtemp, name, strlen (name) + 1);
				std::wstring w_name = wtemp;

				Task_pointer * tp = melt_d;
				painter->set_max_resolution (800);
				painter->set_task (tp);
				painter->set_min_max (0.0, 1.0);
				painter->draw_field (0, COLOR_SCALE_RED_YELLOW);
				painter->draw_contour_lines (0, 20);
				painter->set_iso_line_width (1.0);
				painter->add_isoline (PC_no_theta / (theta + h * theta));
				painter->draw_to_file (w_name);
				painter->reset ();

				// front
				char front_name[128];
				sprintf (front_name, "Result//Melt_SF//RT//iter_%i_theta_%lf.txt", iter, theta + h * theta);
				FILE * file_front_loc = fopen (front_name, "w");
				double fp;
				for (int i_r = 0; i_r < n_receivers; i_r++)
				{
					fp = melt->get_front_position (front_pos[i_r].x);
					fprintf (file_front_loc, "%lf\n", fp);
				}
				fclose (file_front_loc);
			}
			// calc r
			{
				double fp;
				r = i = 0.0;
				for (int i_r = 0; i_r < n_receivers; i_r++)
				{
					fp = melt->get_front_position (front_pos[i_r].x);
					i += pow (fp - front_pos[i_r].y, 2.0);
					R->setElem (i_r, fp - front_pos[i_r].y);
				}
				s = pow (theta - theta_av, 2.0);
				r = i;
			}
			// calc F
			{
				double fp, fp_der;
				for (int i_r = 0; i_r < n_receivers; i_r++)
				{
					fp = melt->get_front_position (front_pos[i_r].x);
					fp_der = melt_d->get_front_position (front_pos[i_r].x);
					F->setElem (i_r, (fp_der - fp) / (h * theta));
				}
			}

			// save current result
			FILE * file_RT = fopen ("Result//Melt_SF//RT//rt.txt", "a+");
			fprintf (file_RT, "%i %e %e %e %e %.16lf %lf %lf %lf %lf\n", iter, r, i, s, alpha, theta, GR, Ste, PC, phys_T0_F);
			printf ("%i %e %e %e %e %.16lf %lf %lf %lf %lf\n", iter, r, i, s, alpha, theta, GR, Ste, PC, phys_T0_F);
			for (int i_r = 0; i_r < n_receivers; i_r++)
			{
				fprintf (file_RT, "%.16lf %.16lf\n", R->getElem (i_r), F->getElem (i_r));
			}
			fclose (file_RT);

			// calc details for the next iteration
			if (fabs (r - r_prev) < (FIND_THETA_PRECIS / 10.0) || r < FIND_THETA_PRECIS)
				found = true;
			if (!found)
			{
				r_prev = r;
				// calc new direction and alpha
				if (iter > 0)
				{
					alpha = 1e-4 * i / s;
				}

				L = F->Scalar_Product (*R) + alpha * (theta - theta_av);
				// new theta
				theta = theta_prev - beta * L / (F->Scalar_Product (*F) + alpha);

				FILE * file_next = fopen ("Result//Melt_SF//RT//next.txt", "a+");
				fprintf (file_next, "next iter: %e %lf\n", alpha, theta);
				fclose (file_next);
			}

			delete melt, melt_d;
		}
	}
}

void Task_Find_Theta::solve_ (char * file_front_name, Painter * painter)
{
	// theta boundaries
	double theta_min = 6.9;
	double theta_max = 12.1;
	double theta_av = (theta_max + theta_min) / 2.0;
	double theta = theta_av;

	// parameters
	double GR_no_theta, Ste_no_theta, PC_no_theta, phys_T0_F_no_Theta;
	GR_no_theta = 570.491136260751;
	Ste_no_theta = 0.00492537313432836;
	PC_no_theta = 22.0 - 21.0;
	phys_T0_F_no_Theta = 23.0 - 21.0;

	// get front position
	int n_receivers;
	FILE * file_front = fopen (file_front_name, "r");
	std::vector<receiver_data> front_pos;
	{
		double x, y;
		while (!feof (file_front))
		{

			fscanf (file_front, "%lf %lf", &x, &y);
			front_pos.push_back ({ x, y });
		}
		front_pos.pop_back ();
		n_receivers = (int)front_pos.size ();
		printf ("Got %i receivers\n", n_receivers);
	}
	fclose (file_front);

	// tasks
	Task_Melt_SF<Triangular_Mesh> * melt;
	Task_Melt_SF<Triangular_Mesh> * melt_d;
	Task_Melt_SF<Triangular_Mesh> * melt_d_b;

	//
	MathVector * R = new MathVector (n_receivers); // R = residual
	MathVector * D = new MathVector (n_receivers); // D = derivative
	MathVector * W = new MathVector (n_receivers); // W = weights

	{
		double w;
		FILE * file_weights = fopen ("Source Files//Melt_SF//weights.txt", "r");
		for (int i_r = 0; i_r < n_receivers; i_r++)
		{
			fscanf (file_weights, "%lf", &w);
			W->setElem (i_r, w);
		}
		fclose (file_weights);
	}

	bool found = false;
	double h = 0.025;
	double del_theta;
	double func_av;
	double alpha = 0.0;
	double r;

	{
		FILE * file_RT = fopen ("Result//Melt_SF//RT//rt.txt", "a+");
		fprintf (file_RT, "iter theta r del_theta func_av func_sum alpha\n");
		fclose (file_RT);
	}

	{
		FILE * file_beta = fopen ("Result//Melt_SF//RT//beta.txt", "a+");
		fprintf (file_beta, "iter beta theta r\n");
		fclose (file_beta);
	}


	// 0 - iter
	// calc new parameters
	{
		double GR = GR_no_theta * theta;
		double Ste = Ste_no_theta * theta;
		double PC = PC_no_theta / theta;
		double phys_T0_F = phys_T0_F_no_Theta / theta;

		melt = new Task_Melt_SF<Triangular_Mesh> ("Source Files//Melt_SF//Parameters.txt");
		char * bound[] = { "Source Files//Melt_SF//Conditions_T.txt",
			"Source Files//Melt_SF//Conditions_P.txt",
			"Source Files//Melt_SF//Conditions_W.txt" };
		melt->prepare ("MeshData//Nodes.txt", "MeshData//Triangles.txt", bound, "Source Files//Melt_SF//time_mapping.txt", "Source Files//Melt_SF//time_stamps.txt");
		melt->get_scales ("Source Files//Melt_SF//Scales.txt");
		melt->reset_parameters (theta, GR, Ste, PC, phys_T0_F);
		int param[][5] = { { SOLVER_METHOD_GMRES, SOLVER_DECOMP_TYPE_LU, 30, 4, SOLVER_MKL_NO },
		{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 4, 4, SOLVER_MKL_NO },
		{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 4, 4, SOLVER_MKL_NO } };
		melt->solve_task (param);
	}

	for (int iter = 0; iter < FIND_THETA_MAX_ITER && !found; iter++)
	{
#pragma omp parallel sections
		{
			// solve task with current theta
			// solve task for derivatives
#pragma omp section
			{
				melt_d = new  Task_Melt_SF<Triangular_Mesh> ("Source Files//Melt_SF//Parameters.txt");
				double GR_der = GR_no_theta * (theta + h * theta);
				double Ste_der = Ste_no_theta * (theta + h * theta);
				double PC_der = PC_no_theta / (theta + h * theta);
				double phys_T0_F_der = phys_T0_F_no_Theta / (theta + h * theta);
				char * bound[] = { "Source Files//Melt_SF//Conditions_T.txt",
					"Source Files//Melt_SF//Conditions_P.txt",
					"Source Files//Melt_SF//Conditions_W.txt" };
				melt_d->prepare ("MeshData//Nodes.txt", "MeshData//Triangles.txt", bound, "Source Files//Melt_SF//time_mapping.txt", "Source Files//Melt_SF//time_stamps.txt");
				melt_d->get_scales ("Source Files//Melt_SF//Scales.txt");
				melt_d->reset_parameters (theta + h * theta, GR_der, Ste_der, PC_der, phys_T0_F_der);
				int param[][5] = { { SOLVER_METHOD_GMRES, SOLVER_DECOMP_TYPE_LU, 30, 3, SOLVER_MKL_NO },
				{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 2, 3, SOLVER_MKL_NO },
				{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 2, 3, SOLVER_MKL_NO } };
				melt_d->solve_task (param);
			}
			// solve tasks for derivatives
#pragma omp section
			{
				melt_d_b = new  Task_Melt_SF<Triangular_Mesh> ("Source Files//Melt_SF//Parameters.txt");
				double GR_der = GR_no_theta * (theta - h * theta);
				double Ste_der = Ste_no_theta * (theta - h * theta);
				double PC_der = PC_no_theta / (theta - h * theta);
				double phys_T0_F_der = phys_T0_F_no_Theta / (theta - h * theta);
				char * bound[] = { "Source Files//Melt_SF//Conditions_T.txt",
					"Source Files//Melt_SF//Conditions_P.txt",
					"Source Files//Melt_SF//Conditions_W.txt" };
				melt_d_b->prepare ("MeshData//Nodes.txt", "MeshData//Triangles.txt", bound, "Source Files//Melt_SF//time_mapping.txt", "Source Files//Melt_SF//time_stamps.txt");
				melt_d_b->get_scales ("Source Files//Melt_SF//Scales.txt");
				melt_d_b->reset_parameters (theta - h * theta, GR_der, Ste_der, PC_der, phys_T0_F_der);
				int param[][5] = { { SOLVER_METHOD_GMRES, SOLVER_DECOMP_TYPE_LU, 30, 3, SOLVER_MKL_NO },
				{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 2, 3, SOLVER_MKL_NO },
				{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 2, 3, SOLVER_MKL_NO } };
				melt_d_b->solve_task (param);
			}
		}
		// save pics and full front position for future references
		{
			// pic
			char name[128];
			wchar_t wtemp[128];

			sprintf (name, "Pictures//Melt_SF//RT//iter_%i_theta_%lf.png", iter, theta);
			mbstowcs (wtemp, name, strlen (name) + 1);
			std::wstring w_name = wtemp;

			Task_pointer * tp = melt;
			painter->set_max_resolution (800);
			painter->set_task (tp);
			painter->set_min_max (0.0, 1.0);
			painter->draw_field (0, COLOR_SCALE_RED_YELLOW);
			painter->draw_contour_lines (0, 20);
			painter->set_iso_line_width (1.0);
			painter->add_isoline (PC_no_theta / theta);
			painter->draw_to_file (w_name);
			painter->reset ();

			// front
			char front_name[128];
			sprintf (front_name, "Result//Melt_SF//RT//iter_%i_theta_%lf.txt", iter, theta);
			FILE * file_front_loc = fopen (front_name, "w");
			double fp;
			for (int i_r = 0; i_r < n_receivers; i_r++)
			{
				fp = melt->get_front_position (front_pos[i_r].x);
				fprintf (file_front_loc, "%.16lf\n", fp);
			}
			fclose (file_front_loc);
		}
		{
			// pic
			char name[128];
			wchar_t wtemp[128];

			sprintf (name, "Pictures//Melt_SF//RT//iter_%i_theta_%lf.png", iter, theta + h * theta);
			mbstowcs (wtemp, name, strlen (name) + 1);
			std::wstring w_name = wtemp;
			
			Task_pointer * tp = melt_d;
			painter->set_max_resolution (800);
			painter->set_task (tp);
			painter->set_min_max (0.0, 1.0);
			painter->draw_field (0, COLOR_SCALE_RED_YELLOW);
			painter->draw_contour_lines (0, 20);
			painter->set_iso_line_width (1.0);
			painter->add_isoline (PC_no_theta / (theta + h * theta));
			painter->draw_to_file (w_name);
			painter->reset ();

			// front
			char front_name[128];
			sprintf (front_name, "Result//Melt_SF//RT//iter_%i_theta_%lf.txt", iter, theta + h * theta);
			FILE * file_front_loc = fopen (front_name, "w");
			double fp;
			for (int i_r = 0; i_r < n_receivers; i_r++)
			{
				fp = melt_d->get_front_position (front_pos[i_r].x);
				fprintf (file_front_loc, "%.16lf\n", fp);
			}
			fclose (file_front_loc);
		}
		{
			// pic
			char name[128];
			wchar_t wtemp[128];

			sprintf (name, "Pictures//Melt_SF//RT//iter_%i_theta_%lf.png", iter, theta - h * theta);
			mbstowcs (wtemp, name, strlen (name) + 1);
			std::wstring w_name = wtemp;

			Task_pointer * tp = melt_d_b;
			painter->set_max_resolution (800);
			painter->set_task (tp);
			painter->set_min_max (0.0, 1.0);
			painter->draw_field (0, COLOR_SCALE_RED_YELLOW);
			painter->draw_contour_lines (0, 20);
			painter->set_iso_line_width (1.0);
			painter->add_isoline (PC_no_theta / (theta - h * theta));
			painter->draw_to_file (w_name);
			painter->reset ();

			// front
			char front_name[128];
			sprintf (front_name, "Result//Melt_SF//RT//iter_%i_theta_%lf.txt", iter, theta - h * theta);
			FILE * file_front_loc = fopen (front_name, "w");
			double fp;
			for (int i_r = 0; i_r < n_receivers; i_r++)
			{
				fp = melt_d_b->get_front_position (front_pos[i_r].x);
				fprintf (file_front_loc, "%.16lf\n", fp);
			}
			fclose (file_front_loc);
		}
		// get del_theta
		{
			double fp, fp_der, fp_der_b;
			func_av = r = 0.0;
			for (int i_r = 0; i_r < n_receivers; i_r++)
			{
				fp = melt->get_front_position (front_pos[i_r].x);
				fp_der = melt_d->get_front_position (front_pos[i_r].x);
				fp_der_b = melt_d_b->get_front_position (front_pos[i_r].x);
				// R
				R->setElem (i_r, W->getElem (i_r) * (fp - front_pos[i_r].y));
				// D
				D->setElem (i_r, W->getElem (i_r) * (fp_der - fp_der_b) / (2.0 * theta * h));
			}
			// get residual
			r = R->Scalar_Product (*R);
			func_av = pow (theta - theta_av, 2.0);

			// alpha
			if (fabs (theta - theta_av) > 1e-7)
				alpha = 1e-3 * r / func_av;
			del_theta = -(D->Scalar_Product (*R) + alpha * (theta - theta_av)) /
				(D->Scalar_Product (*D) + alpha);
		}
		delete melt_d, melt_d_b;

		{
			FILE * file_RT = fopen ("Result//Melt_SF//RT//rt.txt", "a+");
			fprintf (file_RT, "%i %.16lf %e %e %e %e %.16lf\n", iter, theta, del_theta, r, func_av, r + func_av, alpha);
			fclose (file_RT);
		}

		if (r < FIND_THETA_PRECIS)
			break;

		// beta cycle
		{
			double beta = 1.0;
			bool next_try = true;
			double residual = 0.0;
			int direction = 0;
			int degree = 1;
			while (next_try && direction < 2)
			{
				// new theta
				double theta_try = theta + beta * del_theta;
				// check the boundaries
				if (theta_min < theta_try && theta_try < theta_max)
				{
					// solve task for it
					{
						delete melt;
						double GR = GR_no_theta * theta_try;
						double Ste = Ste_no_theta * theta_try;
						double PC = PC_no_theta / theta_try;
						double phys_T0_F = phys_T0_F_no_Theta / theta_try;

						melt = new Task_Melt_SF<Triangular_Mesh> ("Source Files//Melt_SF//Parameters.txt");
						char * bound[] = { "Source Files//Melt_SF//Conditions_T.txt",
							"Source Files//Melt_SF//Conditions_P.txt",
							"Source Files//Melt_SF//Conditions_W.txt" };
						melt->prepare ("MeshData//Nodes.txt", "MeshData//Triangles.txt", bound, "Source Files//Melt_SF//time_mapping.txt", "Source Files//Melt_SF//time_stamps.txt");
						melt->get_scales ("Source Files//Melt_SF//Scales.txt");
						melt->reset_parameters (theta_try, GR, Ste, PC, phys_T0_F);
						int param[][5] = { { SOLVER_METHOD_GMRES, SOLVER_DECOMP_TYPE_LU, 30, 4, SOLVER_MKL_NO },
						{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 4, 4, SOLVER_MKL_NO },
						{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 4, 4, SOLVER_MKL_NO } };
						melt->solve_task (param);
					}

					// get main residual 
					double fp;
					residual = 0.0;
					for (int i_r = 0; i_r < n_receivers; i_r++)
					{
						fp = melt->get_front_position (front_pos[i_r].x);
						residual += pow (W->getElem (i_r) *(fp - front_pos[i_r].y), 2.0);
					}
				}
				else
					residual = 1e+10;
				{
					FILE * file_beta = fopen ("Result//Melt_SF//RT//beta.txt", "a+");
					fprintf (file_beta, "%i %.16lf %.16lf %e\n", iter, beta, theta_try, residual);
					fclose (file_beta);
				}
				// if residual is bigger, then need to make beta smaller
				if (residual > r)
				{
					beta = beta / pow (2.0, degree);
					degree++;
				}
				else
				{
					next_try = false;
				}
				// if beta is too small, no next_try
				if (beta < BETA_PRECIS)
				{
					if (direction == 0)
					{
						direction++;
						del_theta = -del_theta; // try in other direction
						beta = 1.0;
						degree = 1;
					}
				}
			}

			// we got our beta
			// if beta is too small, stagnation, no new solution possible
			if (beta < BETA_PRECIS)
			{
				found = true;
			}

			// if residual is small, then we got our solution
			if (residual < FIND_THETA_PRECIS)
			{
				found = true;
			}
			// otherwise go to the next cycle iteration
			theta = theta + beta * del_theta;
		}
	}
}

void Task_Find_Theta::solve_ (char * file_front_name, Painter * painter, double alpha)
{
	// theta boundaries
	double theta_min = 6.9;
	double theta_max = 12.1;
	double theta_av = (theta_max + theta_min) / 2.0;
	double theta = theta_av;

	// parameters
	double GR_no_theta, Ste_no_theta, PC_no_theta, phys_T0_F_no_Theta;
	GR_no_theta = 570.491136260751;
	Ste_no_theta = 0.00492537313432836;
	PC_no_theta = 22.0 - 21.0;
	phys_T0_F_no_Theta = 23.0 - 21.0;

	// get front position
	int n_receivers;
	FILE * file_front = fopen (file_front_name, "r");
	std::vector<receiver_data> front_pos;
	{
		double x, y;
		while (!feof (file_front))
		{

			fscanf (file_front, "%lf %lf", &x, &y);
			front_pos.push_back ({ x, y });
		}
		front_pos.pop_back ();
		n_receivers = (int)front_pos.size ();
		printf ("Got %i receivers\n", n_receivers);
	}
	fclose (file_front);

	// tasks
	Task_Melt_SF<Triangular_Mesh> * melt;
	Task_Melt_SF<Triangular_Mesh> * melt_d;
	Task_Melt_SF<Triangular_Mesh> * melt_d_b;

	//
	MathVector * R = new MathVector (n_receivers); // R = residual
	MathVector * D = new MathVector (n_receivers); // D = derivative
	MathVector * W = new MathVector (n_receivers); // W = weights

	{
		double w;
		FILE * file_weights = fopen ("Source Files//Melt_SF//weights.txt", "r");
		for (int i_r = 0; i_r < n_receivers; i_r++)
		{
			fscanf (file_weights, "%lf", &w);
			W->setElem (i_r, w);
		}
		fclose (file_weights);
	}

	bool found = false;
	double h = 0.025;
	double del_theta;
	double func_av;
	double r;

	{
		char name[128];
		sprintf (name, "Result//Melt_SF//RT//alpha_%e_rt.txt", alpha);
		FILE * file_RT = fopen (name, "a+");
		fprintf (file_RT, "iter theta r func_av func_sum alpha\n");
		fclose (file_RT);
	}

	{
		char name[128];
		sprintf (name, "Result//Melt_SF//RT//alpha_%e_beta.txt", alpha);
		FILE * file_beta = fopen (name, "a+");
		fprintf (file_beta, "iter beta theta r\n");
		fclose (file_beta);
	}


	// 0 - iter
	// calc new parameters
	{
		double GR = GR_no_theta * theta;
		double Ste = Ste_no_theta * theta;
		double PC = PC_no_theta / theta;
		double phys_T0_F = phys_T0_F_no_Theta / theta;

		melt = new Task_Melt_SF<Triangular_Mesh> ("Source Files//Melt_SF//Parameters.txt");
		char * bound[] = { "Source Files//Melt_SF//Conditions_T.txt",
			"Source Files//Melt_SF//Conditions_P.txt",
			"Source Files//Melt_SF//Conditions_W.txt" };
		melt->prepare ("MeshData//Nodes.txt", "MeshData//Triangles.txt", bound, "Source Files//Melt_SF//time_mapping.txt", "Source Files//Melt_SF//time_stamps.txt");
		melt->get_scales ("Source Files//Melt_SF//Scales.txt");
		melt->reset_parameters (theta, GR, Ste, PC, phys_T0_F);
		int param[][5] = { { SOLVER_METHOD_GMRES, SOLVER_DECOMP_TYPE_LU, 30, 4, SOLVER_MKL_NO },
		{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 4, 4, SOLVER_MKL_NO },
		{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 4, 4, SOLVER_MKL_NO } };
		melt->solve_task (param);
	}

	for (int iter = 0; iter < FIND_THETA_MAX_ITER && !found; iter++)
	{

		//// for time testing	
		//clock_t t1, t2;
		//t1 = clock ();
		//// for time testing	

#pragma omp parallel sections
		{
			// solve task with current theta
			// solve task for derivatives
#pragma omp section
			{
				melt_d = new  Task_Melt_SF<Triangular_Mesh> ("Source Files//Melt_SF//Parameters.txt");
				double GR_der = GR_no_theta * (theta + h * theta);
				double Ste_der = Ste_no_theta * (theta + h * theta);
				double PC_der = PC_no_theta / (theta + h * theta);
				double phys_T0_F_der = phys_T0_F_no_Theta / (theta + h * theta);
				char * bound[] = { "Source Files//Melt_SF//Conditions_T.txt",
					"Source Files//Melt_SF//Conditions_P.txt",
					"Source Files//Melt_SF//Conditions_W.txt" };
				melt_d->prepare ("MeshData//Nodes.txt", "MeshData//Triangles.txt", bound, "Source Files//Melt_SF//time_mapping.txt", "Source Files//Melt_SF//time_stamps.txt");
				melt_d->get_scales ("Source Files//Melt_SF//Scales.txt");
				melt_d->reset_parameters (theta + h * theta, GR_der, Ste_der, PC_der, phys_T0_F_der);
				int param[][5] = { { SOLVER_METHOD_GMRES, SOLVER_DECOMP_TYPE_LU, 30, 3, SOLVER_MKL_NO },
				{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 2, 3, SOLVER_MKL_NO },
				{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 2, 3, SOLVER_MKL_NO } };
				melt_d->solve_task (param);
			}
			// solve tasks for derivatives
#pragma omp section
			{
				melt_d_b = new  Task_Melt_SF<Triangular_Mesh> ("Source Files//Melt_SF//Parameters.txt");
				double GR_der = GR_no_theta * (theta - h * theta);
				double Ste_der = Ste_no_theta * (theta - h * theta);
				double PC_der = PC_no_theta / (theta - h * theta);
				double phys_T0_F_der = phys_T0_F_no_Theta / (theta - h * theta);
				char * bound[] = { "Source Files//Melt_SF//Conditions_T.txt",
					"Source Files//Melt_SF//Conditions_P.txt",
					"Source Files//Melt_SF//Conditions_W.txt" };
				melt_d_b->prepare ("MeshData//Nodes.txt", "MeshData//Triangles.txt", bound, "Source Files//Melt_SF//time_mapping.txt", "Source Files//Melt_SF//time_stamps.txt");
				melt_d_b->get_scales ("Source Files//Melt_SF//Scales.txt");
				melt_d_b->reset_parameters (theta - h * theta, GR_der, Ste_der, PC_der, phys_T0_F_der);
				int param[][5] = { { SOLVER_METHOD_GMRES, SOLVER_DECOMP_TYPE_LU, 30, 3, SOLVER_MKL_NO },
				{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 2, 3, SOLVER_MKL_NO },
				{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 2, 3, SOLVER_MKL_NO } };
				melt_d_b->solve_task (param);
			}
		}

		//// for time testing	
		//t2 = clock ();
		//printf ("%lf seconds\n", (double)(t2 - t1) / 1000.0);
		//break;
		//// for time testing	

		// save pics and full front position for future references
		{
			// pic
			char name[128];
			wchar_t wtemp[128];

			sprintf (name, "Pictures//Melt_SF//RT//alpha_%e_iter_%i_theta_%lf.png", alpha, iter, theta);
			mbstowcs (wtemp, name, strlen (name) + 1);
			std::wstring w_name = wtemp;

			Task_pointer * tp = melt;
			painter->set_max_resolution (800);
			painter->set_task (tp);
			painter->set_min_max (0.0, 1.0);
			painter->draw_field (0, COLOR_SCALE_RED_YELLOW);
			painter->draw_contour_lines (0, 20);
			painter->set_iso_line_width (1.0);
			painter->add_isoline (PC_no_theta / theta);
			painter->draw_to_file (w_name);
			painter->reset ();

			// front
			char front_name[128];
			sprintf (front_name, "Result//Melt_SF//RT//alpha_%e_iter_%i_theta_%lf.txt", alpha, iter, theta);
			FILE * file_front_loc = fopen (front_name, "w");
			double fp;
			for (int i_r = 0; i_r < n_receivers; i_r++)
			{
				fp = melt->get_front_position (front_pos[i_r].x);
				fprintf (file_front_loc, "%.16lf\n", fp);
			}
			fclose (file_front_loc);
		}
		{
			// pic
			char name[128];
			wchar_t wtemp[128];

			sprintf (name, "Pictures//Melt_SF//RT//alpha_%e_iter_%i_theta_%lf.png", alpha, iter, theta + h * theta);
			mbstowcs (wtemp, name, strlen (name) + 1);
			std::wstring w_name = wtemp;

			Task_pointer * tp = melt_d;
			painter->set_max_resolution (800);
			painter->set_task (tp);
			painter->set_min_max (0.0, 1.0);
			painter->draw_field (0, COLOR_SCALE_RED_YELLOW);
			painter->draw_contour_lines (0, 20);
			painter->set_iso_line_width (1.0);
			painter->add_isoline (PC_no_theta / (theta + h * theta));
			painter->draw_to_file (w_name);
			painter->reset ();

			// front
			char front_name[128];
			sprintf (front_name, "Result//Melt_SF//RT//alpha_%e_iter_%i_theta_%lf.txt", alpha, iter, theta + h * theta);
			FILE * file_front_loc = fopen (front_name, "w");
			double fp;
			for (int i_r = 0; i_r < n_receivers; i_r++)
			{
				fp = melt_d->get_front_position (front_pos[i_r].x);
				fprintf (file_front_loc, "%.16lf\n", fp);
			}
			fclose (file_front_loc);
		}
		{
			// pic
			char name[128];
			wchar_t wtemp[128];

			sprintf (name, "Pictures//Melt_SF//RT//alpha_%e_iter_%i_theta_%lf.png", alpha, iter, theta - h * theta);
			mbstowcs (wtemp, name, strlen (name) + 1);
			std::wstring w_name = wtemp;

			Task_pointer * tp = melt_d_b;
			painter->set_max_resolution (800);
			painter->set_task (tp);
			painter->set_min_max (0.0, 1.0);
			painter->draw_field (0, COLOR_SCALE_RED_YELLOW);
			painter->draw_contour_lines (0, 20);
			painter->set_iso_line_width (1.0);
			painter->add_isoline (PC_no_theta / (theta - h * theta));
			painter->draw_to_file (w_name);
			painter->reset ();

			// front
			char front_name[128];
			sprintf (front_name, "Result//Melt_SF//RT//alpha_%e_iter_%i_theta_%lf.txt", alpha, iter, theta - h * theta);
			FILE * file_front_loc = fopen (front_name, "w");
			double fp;
			for (int i_r = 0; i_r < n_receivers; i_r++)
			{
				fp = melt_d_b->get_front_position (front_pos[i_r].x);
				fprintf (file_front_loc, "%.16lf\n", fp);
			}
			fclose (file_front_loc);
		}
		// get del_theta
		{
			double fp, fp_der, fp_der_b;
			func_av = r = 0.0;
			for (int i_r = 0; i_r < n_receivers; i_r++)
			{
				fp = melt->get_front_position (front_pos[i_r].x);
				fp_der = melt_d->get_front_position (front_pos[i_r].x);
				fp_der_b = melt_d_b->get_front_position (front_pos[i_r].x);
				// R
				R->setElem (i_r, W->getElem (i_r) * (fp - front_pos[i_r].y));
				// D
				D->setElem (i_r, W->getElem (i_r) * (fp_der - fp_der_b) / (2.0 * theta * h));
			}
			// get residual
			r = R->Scalar_Product (*R);
			func_av = pow (theta - theta_av, 2.0);

			// alpha
			del_theta = -(D->Scalar_Product (*R) + alpha * (theta - theta_av)) /
				(D->Scalar_Product (*D) + alpha);
		}
		delete melt_d, melt_d_b;

		{
			char name[128];
			sprintf (name, "Result//Melt_SF//RT//alpha_%e_rt.txt", alpha);
			FILE * file_RT = fopen (name, "a+");
			fprintf (file_RT, "%i %.16lf %e %e %e %.16lf\n", iter, theta, r, func_av, r + func_av, alpha);
			fclose (file_RT);
		}

		if (r < FIND_THETA_PRECIS)
			break;

		// beta cycle
		{
			double beta = 1.0;
			bool next_try = true;
			double residual = 0.0;
			int direction = 0;
			int degree = 1;
			while (next_try && direction < 2)
			{
				// new theta
				double theta_try = theta + beta * del_theta;
				// check the boundaries
				if (theta_min < theta_try && theta_try < theta_max)
				{
					// solve task for it
					{
						delete melt;
						double GR = GR_no_theta * theta_try;
						double Ste = Ste_no_theta * theta_try;
						double PC = PC_no_theta / theta_try;
						double phys_T0_F = phys_T0_F_no_Theta / theta_try;

						melt = new Task_Melt_SF<Triangular_Mesh> ("Source Files//Melt_SF//Parameters.txt");
						char * bound[] = { "Source Files//Melt_SF//Conditions_T.txt",
							"Source Files//Melt_SF//Conditions_P.txt",
							"Source Files//Melt_SF//Conditions_W.txt" };
						melt->prepare ("MeshData//Nodes.txt", "MeshData//Triangles.txt", bound, "Source Files//Melt_SF//time_mapping.txt", "Source Files//Melt_SF//time_stamps.txt");
						melt->get_scales ("Source Files//Melt_SF//Scales.txt");
						melt->reset_parameters (theta_try, GR, Ste, PC, phys_T0_F);
						int param[][5] = { { SOLVER_METHOD_GMRES, SOLVER_DECOMP_TYPE_LU, 30, 4, SOLVER_MKL_NO },
						{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 4, 4, SOLVER_MKL_NO },
						{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 4, 4, SOLVER_MKL_NO } };
						melt->solve_task (param);
					}

					// get main residual 
					double fp;
					residual = 0.0;
					for (int i_r = 0; i_r < n_receivers; i_r++)
					{
						fp = melt->get_front_position (front_pos[i_r].x);
						residual += pow (W->getElem (i_r) *(fp - front_pos[i_r].y), 2.0);
					}
				}
				else
					residual = 1e+10;
				{
					char name[128];
					sprintf (name, "Result//Melt_SF//RT//alpha_%e_beta.txt", alpha);
					FILE * file_beta = fopen (name, "a+");
					fprintf (file_beta, "%i %.16lf %.16lf %e\n", iter, beta, theta_try, residual);
					fclose (file_beta);
				}
				// if residual is bigger, then need to make beta smaller
				if (residual > r)
				{
					beta = beta / pow (2.0, degree);
					degree++;
				}
				else
				{
					next_try = false;
				}
				// if beta is too small, no next_try
				if (beta < BETA_PRECIS)
				{
					next_try = false;
					//if (direction == 0)
					//{
					//	del_theta = -del_theta; // try in other direction
					//	beta = 1.0;
					//	degree = 1;
					//}
					//direction++;
				}
			}

			// we got our beta
			// if beta is too small, stagnation, no new solution possible
			if (beta < BETA_PRECIS)
			{
				found = true;
			}

			// if residual is small, then we got our solution
			if (residual < FIND_THETA_PRECIS)
			{
				found = true;
			}
			// otherwise go to the next cycle iteration
			theta = theta + beta * del_theta;
		}
	}

	delete R;
	delete W;
	delete D;
	delete melt;
}
