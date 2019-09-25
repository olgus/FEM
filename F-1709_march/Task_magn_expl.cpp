#include "Task_magn_expl.h"

void Task_magn_expl::initialize ()
{
	sP->Initialize_const (0.0);
	Gamma->Initialize_const (gamma_start);
}

void Task_magn_expl::build_system ()
{
	// black out
	sA->Zero ();
	sB->Zero ();
	 
	int n_elements = mesh.get_n_elements ();
	// go by receivers
	for (int i_receiver = 0; i_receiver < n_received; i_receiver++)
	{
		// go by elements
		for (int k_element = 0; k_element < n_elements; k_element++)
		{
			// go by elements
			for (int m_element = 0; m_element < n_elements; m_element++)
			{
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						sA->addToElem (k_element * 3 + i, m_element * 3 + j, 
							L->Elem (i_receiver * n_elements + k_element, i) * L->Elem (i_receiver * n_elements + m_element, j));
					}
				}
			}

			for (int i = 0; i < 3; i++)
			{
				sB->addToElem (k_element * 3 + i, L->Elem (i_receiver * n_elements + k_element, i) * b_received[i_receiver][3]);

			}
		}
	}
}

void Task_magn_expl::output ()
{
	// write P vector into file
	FILE * file = fopen ("Result//magn//magn_expl_p.txt", "w");
	// go by vector or results
	for (int i = 0, i_end = mesh.get_n_elements (); i < i_end; i++)
	{
		for (int j = i * 3; j < (i + 1) * 3; j++)
		{
			fprintf (file, "%.16lf ", sP->getElem (j));
		}
		fprintf (file, "\n");
	}
	fclose (file);
}

void Task_magn_expl::regularize ()
{
	MathVector * local_A = new MathVector (3);
	MathVector * local_D = new MathVector (3);

	// go by elements
	for (int k_element = 0, k_element_end = mesh.get_n_elements (); k_element < k_element_end; k_element++)
	{
		local_D->Zero ();
		// go by elements
		for (int m_element = 0, m_element_end = mesh.get_n_elements (); m_element < m_element_end; m_element++)
		{
			// check neighbouring state
			if (mesh.check_neighbouring_state (k_element, m_element))
			{
				local_A->Zero ();
				for (int i = 0; i < 3; i++)
				{
					local_A->addToElem (i, Gamma->getElem (k_element * 3 + i));
					local_A->addToElem (i, Gamma->getElem (m_element * 3 + i));
				}
				// save additions for diagonal block
				for (int i = 0; i < 3; i++)
				{
					local_D->addToElem (i, Gamma->getElem (k_element * 3 + i));
					local_D->addToElem (i, Gamma->getElem (m_element * 3 + i));
				}
				// add local_A
				local_A->MultiplyByNumber (-1.0);
				for (int i = 0; i < 3; i++)
				{
					sA->addToElem (k_element * 3 + i, m_element * 3 + i, local_A->getElem (i));
				}
			}
		}
		// add diagonal block
		for (int i = 0; i < 3; i++)
		{
			sA->addToElem (k_element * 3 + i, k_element * 3 + i, local_D->getElem (i));
		}
	}

	delete local_A;
	delete local_D;
}

int Task_magn_expl::adapt ()
{
	double p_sum[3];
	int changed_gammas = 0;
	// go by elements
	for (int k_element = 0, k_element_end = mesh.get_n_elements (); k_element < k_element_end; k_element++)
	{
		for (int i = 0; i < 3; i++)
		{
			p_sum[i] = 0.0;
		}
		// go by elements
		for (int m_element = 0, m_element_end = mesh.get_n_elements (); m_element < m_element_end; m_element++)
		{
			// check neighbouring state
			if (mesh.check_neighbouring_state (k_element, m_element))
			{
				// go by 3 p components
				for (int i = 0; i < 3; i++)
				{
					// check if p absolute is bigger than 0
					if (fabs (sP->getElem (k_element * 3 + i)) > P_ZERO)
					{
						// calc sum of differences in neighbouring cells
						{
							p_sum[i] += (sP->getElem (k_element * 3 + i) - sP->getElem (m_element * 3 + i));
						}
					}
				}
			}
		}
		// go by 3 p components
		for (int i = 0; i < 3; i++)
		{
			// if sum of differences is bigger than desired, mult that gamma by 10
			if (p_sum[i] > DES_DIF)
			{
				Gamma->multiplyElem (k_element * 3 + i, 1e+1);
				changed_gammas++;
				printf ("%i\t", k_element * 3 + i);
			}
		}
	}
	return changed_gammas;
}

double Task_magn_expl::calc_functional ()
{
	double f = 0.0;
	double By = 0.0;
	double p_dif = 0.0;
	// go by receivers
	for (int i = 0; i < n_received; i++)
	{
		// substract calculated B from received
		By = receive (i, true);
		f += pow (By - b_received[i][3], 2.0);
	}
	// go by elements
	for (int k_element = 0, k_element_end = mesh.get_n_elements (); k_element < k_element_end; k_element++)
	{
		// go by elements
		for (int m_element = 0, m_element_end = mesh.get_n_elements (); m_element < m_element_end; m_element++)
		{
			// check neighbouring state
			if (mesh.check_neighbouring_state (k_element, m_element))
			{
				for (int i = 0; i < 3; i++)
				{
					f += Gamma->getElem (k_element * 3 + i) * pow (sP->getElem (k_element * 3 + i) - sP->getElem (m_element * 3 + i), 2.0);
				}
			}
		}
	}
	return f;
}

void Task_magn_expl::initialize_l ()
{
	double coef = 1.0 / (4.0 * PI);

	double p_coef[3];
	double r;
	double MG[3];
	int n_integration_points = mesh.amount_of_integration_points (0);
	double ** integration_points = new double *[n_integration_points];
	for (int i = 0; i < n_integration_points; i++)
	{
		integration_points[i] = new double[3];
	}
	double * weights = new double[n_integration_points];
	double jac = 0.0;

	L->Zero ();
	int n_elements = mesh.get_n_elements ();
	// calc L's
	{
		// go by receivers
		for (int i_receiver = 0; i_receiver < n_received; i_receiver++)
		{
			// go by cells
			for (int i_element = 0; i_element < n_elements; i_element++)
			{
				mesh.integration_points (i_element, integration_points, weights, &jac);
				// go by integration poitns
				for (int i = 0; i < n_integration_points; i++)
				{
					// move M into G and get r (M, G)
					r = 0.0;
					for (int j = 0; j < 3; j++)
					{
						MG[j] = b_received[i_receiver][j] - integration_points[i][j];
						r += pow (MG[j], 2.0);
					}
					r = 1.0 / sqrt (r);
					// p components
					p_coef[0] = 3.0 * MG[0] * MG[1] * pow (r, 2.0);
					p_coef[1] = 3.0 * MG[1] * MG[1] * pow (r, 2.0) - 1.0;
					p_coef[2] = 3.0 * MG[2] * MG[1] * pow (r, 2.0);

					for (int k = 0; k < 3; k++)
					{
						L->addToElem (i_receiver * n_elements + i_element, k, p_coef[k] * weights[i] * pow (r, 3.0) * coef * jac);
					}
				}
			}
		}
	}

	// cleanup
	{
		for (int i = 0; i < n_integration_points; i++)
		{
			delete[] integration_points[i];
		}
		delete[] integration_points;
		delete[] weights;
	}

	
	//double p_coef[3];
	//double r;
	//double MG[3];
	//double center[3];
	//L->Zero ();
	//int n_elements = mesh.get_n_elements ();
	//double coef = 1.0 / (4.0 * PI);
	//double jac;
	//// calc L's
	//{
	//	// go by receivers
	//	for (int i_receiver = 0; i_receiver < n_received; i_receiver++)
	//	{
	//		// go by cells
	//		for (int i_element = 0; i_element < n_elements; i_element++)
	//		{
	//			jac = mesh.elements[i_element]->get_geometrical_area ();
	//			mesh.get_mass_center (i_element, center);
	//			{
	//				// move M into G and get r (M, G)
	//				r = 0.0;
	//				for (int j = 0; j < 3; j++)
	//				{
	//					MG[j] = b_received[i_receiver][j] - center[j];
	//					r += pow (MG[j], 2.0);
	//				}
	//				r = 1.0 / sqrt (r);
	//				// p components
	//				p_coef[0] = 3.0 * MG[0] * MG[0] * pow (r, 2.0) -1.0;
	//				p_coef[1] = 3.0 * MG[1] * MG[0] * pow (r, 2.0);
	//				p_coef[2] = 3.0 * MG[2] * MG[0] * pow (r, 2.0);
	//				for (int k = 0; k < 3; k++)
	//				{
	//					L->addToElem (i_receiver * n_elements + i_element, k, p_coef[k] * coef * pow (r, 3.0) * jac);
	//				}
	//			}
	//		}
	//	}
	//}

	L->FPrint ("Result//magn//L.txt");
}

Task_magn_expl::Task_magn_expl ()
{
	n_areas = 0;
	p = NULL;

	n_received = 0;
	b_received = NULL;

	rec_z_levels = NULL;
	rec_x_y_levels = NULL;
	rec_y_x_levels = NULL;

	sA = NULL;
	sB = NULL;
	sP = NULL;
	Gamma = NULL;
	L = NULL;

	gamma_start = 1e-15;

	FUNCTIONAL_THRESHOLD_PERCENTILE = 0.1;
	DES_DIF = 10.0;
	P_ZERO = 1e-2;
}

Task_magn_expl::Task_magn_expl (const Task_magn_expl & task)
{
}

Task_magn_expl::~Task_magn_expl ()
{
	if (p != NULL)
	{
		for (int i = 0; i < n_areas; i++)
		{
			delete[] p[i];
		}
		delete[] p;
	}
	p = NULL;

	if (b_received != NULL)
	{
		for (int i = 0; i < n_areas; i++)
		{
			delete[] b_received[i];
		}
		delete[] b_received;
	}
	b_received = NULL;

	if (rec_z_levels != NULL)
		delete[] rec_z_levels;
	rec_z_levels = NULL;
	if (rec_x_y_levels != NULL)
		delete[] rec_x_y_levels;
	rec_x_y_levels = NULL;
	if (rec_y_x_levels != NULL)
		delete[] rec_y_x_levels;
	rec_y_x_levels = NULL;
	if (sA != NULL)
		delete sA;
	sA = NULL;
	if (sB != NULL)
		delete sB;
	sB = NULL;
	if (sP != NULL)
		delete sP;
	sP = NULL;
	if (Gamma != NULL)
		delete Gamma;
	Gamma = NULL;
	if (L != NULL)
		delete L;
	L = NULL;
}

void Task_magn_expl::get_receivers_data (char * file_receivers_data)
{
	// read points and receivers values in them
	FILE * file = fopen (file_receivers_data, "r");

	for (int i = 0; i < n_received; i++)
	{
		fscanf (file, "%lf", &b_received[i][3]);
	}
	fclose (file);
}

void Task_magn_expl::save_field (char * file_result, bool calc_p)
{
	FILE * file = fopen (file_result, "w");
	FILE * file_extra = fopen ("Result//magn//receivers.txt", "w");
	double By;
	int n_elements = mesh.get_n_elements ();
	// go by receivers
	for (int i = 0; i < n_received; i++)
	{
		By = receive (i, calc_p);
		fprintf (file, "%.16lf\n", By);
		for (int j = 0; j < 3; j++)
		{
			fprintf (file_extra, "%.16lf ", b_received[i][j]);
		}
		fprintf (file_extra, "\n");
	}
	fclose (file);
	fclose (file_extra);
}

double Task_magn_expl::receive (int i_receiver, bool calc)
{
	double By = 0.0;
	int area;
	double p_local;
	// go by cells
	int n_elements = mesh.get_n_elements ();
	for (int i_element = 0; i_element < n_elements; i_element++)
	{
		area = mesh.get_area (i_element);
		// multiply respected L's 
		for (int k = 0; k < 3; k++)
		{
			if (calc)
			{
				p_local = sP->getElem (i_element * 3 + k);
			}
			else
			{
				p_local = p[area][k];
			}
			By += L->Elem (i_receiver * n_elements + i_element, k) * p_local;
		}
	}
	return By;
}

double Task_magn_expl::solve ()
{
	// set memory
	if (sA == NULL)
		sA = new Matrix (mesh.get_n_elements () * 3);
	if (sB == NULL)
		sB = new MathVector (mesh.get_n_elements () * 3);
	if (sP == NULL)
		sP = new MathVector (mesh.get_n_elements () * 3);
	if (Gamma == NULL)
		Gamma = new MathVector (mesh.get_n_elements () * 3);	

	// initial state
	initialize ();
	// build matrix
	build_system ();
	// create a copy
	Matrix * sA_copy = new Matrix (*sA);
	MathVector * sB_copy = new MathVector (*sB);
	double f, f_start, THRESHOLD;
	f = 1e+20;
	// log
	FILE * log = fopen ("Result//magn//magn_log.txt", "w");
	//sA->FPrint (log);
	//sB->FPrint (log);
	// first step	
	regularize ();
	sA->solve (sP, *sB);
	f_start = calc_functional ();
	THRESHOLD = FUNCTIONAL_THRESHOLD_PERCENTILE * f_start;
	int cg = adapt ();
	fprintf (log, "%e\n", f_start);
	printf ("%e\n", f_start);

	bool stop = false;
	for (int i = 0; i < MAX_STEPS && !stop; i++)
	{
		// get rid of previous regularization
		sA->Copy (*sA_copy);
		sB->Copy (*sB_copy);
		// regularization
		regularize ();
		// solve
		sA->solve (sP, *sB);
		// get functional
		f = calc_functional ();
		fprintf (log, "%i\t%e\n", cg, f);
		printf ("\n%i\t%e\t%i\n", i, f, cg);
		// increase needed gammas
		cg = adapt ();
		// check if functional has increased less than allowed
		if ((f - f_start) > THRESHOLD)
			stop = true;
		if (cg == 0)
			stop = true;
	}
	fclose (log);
	// save best result separately
	output ();

	delete sA_copy;
	delete sB_copy;
	return f;
}

double Task_magn_expl::get_P (int k_element, int k_component)
{
	return sP->getElem (k_element * 3 + k_component);
}

double Task_magn_expl::get_P_true (int k_element, int k_component)
{
	int area = mesh.get_area (k_element);
	return p[area][k_component];
}

void Task_magn_expl::set_parameters (double Gamma_start, double threshold, double des_dif, double p_zero)
{
	gamma_start = Gamma_start;
	FUNCTIONAL_THRESHOLD_PERCENTILE = threshold;
	DES_DIF = des_dif;
	P_ZERO = p_zero;
}

void Task_magn_expl::prepare_field (char * mesh_file, char * file_areas, char * file_area_data)
{
	// make mesh
	int n;
	mesh.build_Mesh (mesh_file, true, &n);
	mesh.reset_areas (file_areas);
	mesh.output ();

	// read data about areas
	FILE * file = fopen (file_area_data, "r");
	fscanf (file, "%i", &n_areas);
	n_areas++;

	// areas go from 0 to n_areas - 1
	if (p == NULL)
	{
		p = new double *[n_areas];
		for (int i = 0; i < n_areas; i++)
		{
			p[i] = new double[3];
		}
	}
	p[0][0] = 0.0;
	p[0][1] = 0.0;
	p[0][2] = 0.0;
	for (int i = 1; i < n_areas; i++)
	{
		fscanf (file, "%lf %lf %lf", &p[i][0], &p[i][1], &p[i][2]);
	}
	fclose (file);
	
	// set Ls
	if (L == NULL)
		L = new Matrix (n_received * mesh.get_n_elements (), 3);
	initialize_l ();
}

void Task_magn_expl::set_receivers_data (char * file_receivers_data)
{
	FILE * file = fopen (file_receivers_data, "r");

	// amount of z levels
	fscanf (file, "%i", &rec_n_z_levels);
	// their values
	rec_z_levels = new double[rec_n_z_levels];
	for (int i = 0; i < rec_n_z_levels; i++)
	{
		fscanf (file, "%lf", &rec_z_levels[i]);
	}

	// amount of x receiver levels
	// boundaries 
	fscanf (file, "%i", &rec_n_x_levels);
	if (rec_n_x_levels > 0)
	{
	fscanf (file, "%lf %lf", &rec_x0, &rec_xN);
	fscanf (file, "%i", &rec_n_x_by_level);

	rec_x_y_levels = new double[rec_n_x_levels];
	// y components of the receivers
	for (int i = 0; i < rec_n_x_levels; i++)
	{
		fscanf (file, "%lf", &rec_x_y_levels[i]);
	}
	}
	// amount of y receiver levels
	// boundaries 
	fscanf (file, "%i", &rec_n_y_levels);
	if (rec_n_y_levels > 0)
	{
		fscanf (file, "%lf %lf", &rec_y0, &rec_yN);
		fscanf (file, "%i", &rec_n_y_by_level);

		rec_y_x_levels = new double[rec_n_y_levels];
		// x components of the receivers
		for (int i = 0; i < rec_n_y_levels; i++)
		{
			fscanf (file, "%lf", &rec_y_x_levels[i]);
		}
	}
	fclose (file);

	n_received = rec_n_z_levels * (rec_n_x_levels * rec_n_x_by_level + rec_n_y_levels * rec_n_y_by_level);

	if (b_received == NULL)
	{
		b_received = new double *[n_received];
		for (int i = 0; i < n_received; i++)
		{
			b_received[i] = new double[4]; // point and y value
		}
	}

	// calc receivers positions
	int rec_cur = 0;
	// go by z levels
	double x, y, z;
	for (int k = 0; k < rec_n_z_levels; k++)
	{
		z = rec_z_levels[k];
		// do x levels
		{
			// calc x step 
			double x_step = (rec_xN - rec_x0) / (rec_n_x_by_level - 1);
			// go by amount of levels
			for (int j = 0; j < rec_n_x_levels; j++)
			{
				// get y position
				y = rec_x_y_levels[j];
				// go by the amount on the level
				for (int i = 0; i < rec_n_x_by_level; i++)
				{
					x = rec_x0 + x_step * i;
					b_received[rec_cur][0] = x;
					b_received[rec_cur][1] = y;
					b_received[rec_cur][2] = z;
					rec_cur++;
				}
			}
		}
		// do y levels
		{
			// calc y step 
			double y_step = (rec_yN - rec_y0) / (rec_n_y_by_level - 1);
			// go by amount of levels
			for (int i = 0; i < rec_n_y_levels; i++)
			{
				// get y position
				x = rec_y_x_levels[i];
				// go by the amount on the level
				for (int j = 0; j < rec_n_y_by_level; j++)
				{
					y = rec_y0 + y_step * j;
					b_received[rec_cur][0] = x;
					b_received[rec_cur][1] = y;
					b_received[rec_cur][2] = z;
					rec_cur++;
				}
			}
		}
	}
	if (n_received != rec_cur)
		printf ("ERROR save_field\n");
}
 