#include "Prototype_Task.h"

template class Task <Mesh<Node<Point_Prototype>, Element>>;
template class Task <Mesh<Node<Point_2D>, Element>>;
template class Task <Mesh<Node<Point_3D>, Element>>;
template class Task <Mesh<Node<Point_2D_Polar>, Element>>;
template class Task <Mesh<Node_2D, Element>>;
template class Task <Mesh<Node_2D_Polar, Element>>;
template class Task <Mesh<Node_1D, Element>>;
template class Task <Mesh<Node_3D, Element>>;

template class Task <Mesh_1D_Hier>;
template class Task <Mesh_1D_Hermitian>;
template class Task <Mesh_1D_L1>;
template class Task <Triangular_Mesh>;
template class Task <Triangle_Vector_Mesh>;
template class Task <Triangular_Polar_Mesh>;
template class Task <Triangular_Mesh_Hier>;
template class Task <Rectangular_Mesh>;
template class Task <Rectangular_Mesh_3>;
template class Task <Rectangular_Mesh_Spline>;
template class Task <Prismatic_Mesh>;
template class Task <Cubic_Mesh>;
template class Task <Mixed_Triangular_Mesh>;

Painter * painter;

template<class Type_Mesh>
double Task<Type_Mesh>::relax (int k_system, MathVector * prev_solution)
{
	{
	//	if (k_system == 1)
	//	{
	//		FILE * file = fopen ("w.txt", "w");
	//		double w = 0.0;
	//		MathVector * sol = new MathVector (eq_matrixes[k_system].Size ());
	//		MathVector * discr = new MathVector (eq_matrixes[k_system].Size ());

	//		for (int i = 0; i < 101; i++)
	//		{
	//			w = i * 0.01;
	//			eq_matrixes[k_system].get_solution (sol); // get a(k)
	//			sol->Linear_Combination (w, *prev_solution, 1.0 - w);
	//			eq_matrixes[k_system].mult_A_v (*sol, discr);
	//			discr->Substract (*(eq_matrixes[k_system].f));
	//			fprintf (file, "%lf\t%e\n", w, discr->Norm ());
	//		}
	//		fclose (file);
	//		delete sol;
	//		delete discr;
	//	}
	}

	double dif = 0.0;
	double w1, w2, golden_ratio_coef, f1, f2,
		lenght_prev,	//длина интервала на предыдущем шаге
		lenght_cur,		//длина текущего интервала
		left, right;	//границы текущего интервала

	int i;
	bool l_r_sw;		//выбор из двух интервалов: [left,x2]||[x1,right]
	golden_ratio_coef = (3.0 - powf (5.0, 0.5)) / 2.0; // коэффициент золотого сечения

													   // границы поиска 
	left = 0.0;
	right = 1.0;

	int par_num_threads = 4;
	double * res_omp = NULL;
	if (par_num_threads > 1)
		res_omp = new double[eq_matrixes[k_system].Size() * (par_num_threads - 1)];
	
	// отпределение начальных точек поиска
	w1 = left + golden_ratio_coef * (right - left);
	w2 = right - golden_ratio_coef * (right - left);

	MathVector * v1 = new MathVector (eq_matrixes[k_system].Size ());
	MathVector * v2 = new MathVector (eq_matrixes[k_system].Size ());

	//расчет невязок f1, f2
	// fi = || A(wi * a(k) - (1.0 - wi) * a(k-1)) - f ||
	eq_matrixes[k_system].get_solution (v1); // get a(k)
	//v1->MultiplyByNumber (w1);
	v2->Linear_Combination (0.0, *prev_solution, 1.0 - w1);
	//v2->Copy (*prev_solution); // get a(k - 1)
	//v2->MultiplyByNumber (1.0 - w1);
	//v1->Add (*v2); 
	// make (wa(k) + (1-w)a(k-1))
	v1->Linear_Combination (w1, *v2, 1.0);
	eq_matrixes[k_system].mult_A_v (*v1, v2, res_omp, par_num_threads);
	v2->Substract (*(eq_matrixes[k_system].f));
	f1 = v2->Norm (); // discr

	eq_matrixes[k_system].get_solution (v1); // get a(k)
	//v1->MultiplyByNumber (w2);
	//v2->Copy (*prev_solution); // get a(k - 1)
	//v2->MultiplyByNumber (1.0 - w2);
	v2->Linear_Combination (0.0, *prev_solution, 1.0 - w2);
	v1->Linear_Combination (w2, *v2, 1.0);
	//v1->Add (*v2); // make (wa(k) + (1-w)a(k-1))
	eq_matrixes[k_system].mult_A_v (*v1, v2, res_omp, par_num_threads);
	v2->Substract (*(eq_matrixes[k_system].f));
	f2 = v2->Norm (); // discr

	if (f1 < f2)		//выбор промежутка, на котором находится экстремум
		l_r_sw = true;
	else
		l_r_sw = false;
	lenght_prev = right - left;

	// точность поиска w
	double tolerance = 1e-5;
	// максимум итераций
	int N = 100;

	for (i = 0; i < N && right - left > tolerance; i++)
	{
		if (l_r_sw)	//вычисление новой точки в интервале, содержащем экстремум
		{
			w1 = left + golden_ratio_coef * (right - left);


			eq_matrixes[k_system].get_solution (v1); // get a(k)
			v2->Linear_Combination (0.0, *prev_solution, 1.0 - w1);
			v1->Linear_Combination (w1, *v2, 1.0);
			eq_matrixes[k_system].mult_A_v (*v1, v2, res_omp, par_num_threads);
			v2->Substract (*(eq_matrixes[k_system].f));
			f1 = v2->Norm (); // discr
		}
		else
		{
			w2 = right - golden_ratio_coef * (right - left);


			eq_matrixes[k_system].get_solution (v1); // get a(k)
			v2->Linear_Combination (0.0, *prev_solution, 1.0 - w2);
			v1->Linear_Combination (w2, *v2, 1.0);
			eq_matrixes[k_system].mult_A_v (*v1, v2, res_omp, par_num_threads);
			v2->Substract (*(eq_matrixes[k_system].f));
			f2 = v2->Norm (); // discr
		}

		if (f1 <= f2)	//сокращение интервала, определение будущего направления движения
		{
			right = w2;
			w2 = w1;
			f2 = f1;
			l_r_sw = true;
		}
		else
		{
			left = w1;
			w1 = w2;
			f1 = f2;
			l_r_sw = false;
		}

		lenght_cur = right - left;
		lenght_prev = lenght_cur;
	}
	// сохр решения
	double w = (left + right) / 2;

	if (w < 6e-5)
		w = 0.0;
	if (w > 0.9999)
		w = 1.0;
	eq_matrixes[k_system].get_solution (v1); // get a(k)
	v1->MultiplyByNumber (w);
	v2->Copy (*prev_solution); // get a(k - 1)
	v2->MultiplyByNumber (1.0 - w);
	v1->Add (*v2); // make (wa(k) + (1-w)a(k-1))
	prev_solution->Copy (*v1);
	eq_matrixes[k_system].mult_A_v (*v1, v2);
	v2->Substract (*(eq_matrixes[k_system].f));
	dif = v2->Norm () / eq_matrixes[k_system].f->Norm (); // discr
	
	//// emergency discr
	//eq_matrixes[k_system].get_solution (v1); // get a(k)
	//double solution_dif = v1->Norm ();
	//if (solution_dif > 1e-15)
	//	solution_dif = 1.0 / solution_dif;
	//else
	//	solution_dif = 1.0;
	//v1->Substract (*prev_solution);
	//solution_dif *= v1->Norm ();

	//// if the solution has changed only by 0.5%, replace the discrepancy
	//if (solution_dif < 0.005 && solution_dif < dif)
	//	dif = solution_dif;

	printf ("\tsystem %i, w: %lf, discr: %e\n", k_system, w, dif);
	fprintf (log_file, "\tw: %lf, discr: %e\n", w, dif);

	delete v1;
	delete v2;

	if (res_omp != NULL)
		delete[] res_omp;
	return dif;
}

template<class Type_Mesh>
bool Task<Type_Mesh>::build_portrait (int k_system)
{
	std::vector<int> listbeg;

	for (int i = 0; i < N_functions; i++)
	{
		listbeg.push_back (-1);
	}
	int listsize = -1;
	std::vector<int> list1;
	std::vector<int> list2;

	int Link1;
	int n_elements = mesh_pointer->get_n_elements ();
	int n_functions;
	int link1, link2; // global numbers of linked functions
	int iaddr;
	// cycle through elements
	for (int i_element = 0; i_element < n_elements; i_element++)
	{
		// get amount of global functions that are not 0 on the element
		n_functions = mesh_pointer->get_amount_non_zero_functions (i_element);
		// cycle through them
		for (int i_function = 0; i_function < n_functions; i_function++)
		{
			// get number of global function for i_function local function
			Link1 = get_function_global_number (i_element, i_function);
			if (Link1 != -1)
			{
				// for following local functions on the element
				for (int j_function = i_function + 1; j_function < n_functions; j_function++)
				{
					link1 = Link1;
					// get global number of j_function
					link2 = get_function_global_number (i_element, j_function);
					if (link2 != -1)
					{
						// sort link1 and link2, so link2 > link1 
						if (link1 > link2)
						{
							link1 = link2;
							link2 = Link1;
						}
						// get place where list for link2 starts
						iaddr = listbeg[link2];
						// if the list hasn't existed
						if (iaddr == -1)
						{
							listsize++;
							listbeg[link2] = listsize;
							list1.push_back (link1);
							list2.push_back (-1);
						}
						else
						{
							// list exists 
							while (list1[iaddr] < link1 && list2[iaddr] > -1) // look for link1 in it
							{
								iaddr = list2[iaddr];
							}
							// next link is bigger than link1
							if (list1[iaddr] > link1)
							{
								listsize++;
								list1.push_back (list1[iaddr]);
								list2.push_back (list2[iaddr]);
								list1[iaddr] = link1;
								list2[iaddr] = listsize;
							}
							else  // end of list
							{
								if (list1[iaddr] < link1) // to skip already existing link
								{
									listsize++;
									list2[iaddr] = listsize;
									list1.push_back (link1);
									list2.push_back (-1);
								}
							}
						}
					}
				}
			}
		}
	}
	listsize++;
	// make portrait
	// set ig, jg
	int * ig = new int[N_functions + 1];
	int * jg = new int[listsize];

	int i = 0;
	ig[0] = 0;
	// go by global functions
	for (i = 0; i < N_functions; i++)
	{
		ig[i + 1] = ig[i];
		iaddr = listbeg[i];
		// till the end of links for current function
		while (iaddr != -1)
		{
			// save the link
			jg[ig[i + 1]] = list1[iaddr];
			// increase ig
			ig[i + 1]++;
			// go to the next link
			iaddr = list2[iaddr];
		}
	}

	// set matrix
	eq_matrixes[k_system].set_matrix_size (N_functions, listsize);
	eq_matrixes[k_system].set_ig_jg (ig, jg);

	// code for test purposes
	//eq_matrixes[k_system].print ();

	eq_matrixes[k_system].Clear ();
	//eq_matrixes[k_system].fprint ();
	bool result = true;
	if (ig[N_functions] != listsize)
		result = false;
	// clear memory
	delete[] ig;
	delete[] jg;

	if (result)
		printf ("\tmatrix portrait has been build successfully.\n\tit's fully automatic and can't go wrong, unlike you, human.\n");
	else
		printf ("ERROR: something went wrong in matrix portrait building.\ni bet you did it on purpose\n");

	return result;
}

template<class Type_Mesh>
bool Task<Type_Mesh>::build_portrait (int k_system, Type_Mesh * mesh)
{
	std::vector<int> listbeg;

	for (int i = 0; i < N_functions; i++)
	{
		listbeg.push_back (-1);
	}
	int listsize = -1;
	std::vector<int> list1;
	std::vector<int> list2;

	int Link1;
	int n_elements = mesh->get_n_elements ();
	int n_functions;
	int link1, link2; // global numbers of linked functions
	int iaddr;
	// cycle through elements
	for (int i_element = 0; i_element < n_elements; i_element++)
	{
		// get amount of global functions that are not 0 on the element
		n_functions = mesh->get_amount_non_zero_functions (i_element);
		// cycle through them
		for (int i_function = 0; i_function < n_functions; i_function++)
		{
			// get number of global function for i_function local function
			Link1 = mesh->get_function_global_number (i_element, i_function);
			if (Link1 != -1)
			{
				// for following local functions on the element
				for (int j_function = i_function + 1; j_function < n_functions; j_function++)
				{
					link1 = Link1;
					// get global number of j_function
					link2 = mesh->get_function_global_number (i_element, j_function);
					if (link2 != -1)
					{
						// sort link1 and link2, so link2 > link1 
						if (link1 > link2)
						{
							link1 = link2;
							link2 = Link1;
						}
						// get place where list for link2 starts
						iaddr = listbeg[link2];
						// if the list hasn't existed
						if (iaddr == -1)
						{
							listsize++;
							listbeg[link2] = listsize;
							list1.push_back (link1);
							list2.push_back (-1);
						}
						else
						{
							// list exists 
							while (list1[iaddr] < link1 && list2[iaddr] > -1) // look for link1 in it
							{
								iaddr = list2[iaddr];
							}
							// next link is bigger than link1
							if (list1[iaddr] > link1)
							{
								listsize++;
								list1.push_back (list1[iaddr]);
								list2.push_back (list2[iaddr]);
								list1[iaddr] = link1;
								list2[iaddr] = listsize;
							}
							else  // end of list
							{
								if (list1[iaddr] < link1) // to skip already existing link
								{
									listsize++;
									list2[iaddr] = listsize;
									list1.push_back (link1);
									list2.push_back (-1);
								}
							}
						}
					}
				}
			}
		}
	}
	listsize++;
	// make portrait
	// set ig, jg
	int * ig = new int[N_functions + 1];
	int * jg = new int[listsize];

	int i = 0;
	ig[0] = 0;
	// go by global functions
	for (i = 0; i < N_functions; i++)
	{
		ig[i + 1] = ig[i];
		iaddr = listbeg[i];
		// till the end of links for current function
		while (iaddr != -1)
		{
			// save the link
			jg[ig[i + 1]] = list1[iaddr];
			// increase ig
			ig[i + 1]++;
			// go to the next link
			iaddr = list2[iaddr];
		}
	}

	// set matrix
	eq_matrixes[k_system].set_matrix_size (N_functions, listsize);
	eq_matrixes[k_system].set_ig_jg (ig, jg);

	// code for test purposes
	//eq_matrixes[k_system].print ();

	eq_matrixes[k_system].Clear ();
	//eq_matrixes[k_system].fprint ();
	bool result = true;
	if (ig[N_functions] != listsize)
		result = false;
	// clear memory
	delete[] ig;
	delete[] jg;

	if (result)
		printf ("\tmatrix portrait has been build successfully.\n\tit's fully automatic and can't go wrong, unlike you, human.\n");
	else
		printf ("ERROR: something went wrong in matrix portrait building.\ni bet you did it on purpose\n");

	return result;
}

template<class Type_Mesh>
void Task<Type_Mesh>::get_conditions (int k_system, char * file_name)
{
	int t, f;
	switch (k_system)
	{
	case 100:
		break;
	default:
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

template <class Type_Mesh>
void Task<Type_Mesh>::prepare (char * mesh_file_name, char * bound_file_name, bool save)
{
	printf ("\tpreparing time layers and elements\n");

	log_file = fopen ("Result Files Extra//log.txt", "w");
	// build a mesh
	mesh_pointer = &mesh;
	prepared = mesh_pointer->build_Mesh (mesh_file_name, save, &N_functions);

	// set amount of equations
	set_n_systems ();
	// set time layers
	set_time_layers ("");
	// prepares elements
	for (int i = 0, i_end = mesh_pointer->get_n_elements (); i < i_end; i++)
	{
		mesh_pointer->prepare_element (i);
	}

	// set conditions
	conditions = new Condition *[n_systems];
	for (int i = 0; i < n_systems; i++)
	{
		conditions[i] = new Condition[mesh_pointer->get_dimentionality () * 2];
		get_conditions (i, bound_file_name);
	}

	use_time_mapping = true;

	// build matrices' portraits
	for (int k_system = 0; k_system < n_systems; k_system++)
		if (!build_portrait (k_system))
		{
			printf ("ERROR: portrait\n");
			prepared = false;
		}

	// if portrait building was successful, set starting conditions
	if (prepared)
		set_starting_conditions ();
}

template<class Type_Mesh>
void Task<Type_Mesh>::prepare (char * mesh_file_name, char * bound_file_name, char * time_file_name, bool save)
{
	printf ("\tpreparing time layers and elements\n");

	log_file = fopen ("Result Files Extra//log.txt", "w");
	// build a mesh
	mesh_pointer = &mesh;
	prepared = mesh_pointer->build_Mesh (mesh_file_name, save, &N_functions);

	// set amount of equations
	set_n_systems ();
	// set time layers
	set_time_layers (time_file_name);
	// prepares elements
	for (int i = 0, i_end = mesh_pointer->get_n_elements (); i < i_end; i++)
	{
		mesh_pointer->prepare_element (i);
	}

	// set conditions
	conditions = new Condition *[n_systems];
	for (int i = 0; i < n_systems; i++)
	{
		conditions[i] = new Condition[mesh_pointer->get_dimentionality () * 2];
		get_conditions (i, bound_file_name);
	}

	use_time_mapping = true;

	// build matrices' portraits
	for (int k_system = 0; k_system < n_systems; k_system++)
		if (!build_portrait (k_system))
		{
			printf ("ERROR: portrait\n");
			prepared = false;
		}

	// if portrait building was successful, set starting conditions
	if (prepared)
		set_starting_conditions ();
}

template<class Type_Mesh>
void Task<Type_Mesh>::prepare (char * mesh_file_name, char * bound_file_name[], char * time_file_name, bool save)
{
	printf ("\tpreparing time layers and elements\n");

	log_file = fopen ("Result Files Extra//log.txt", "w");
	// build a mesh
	mesh_pointer = &mesh;
	prepared = mesh_pointer->build_Mesh (mesh_file_name, save, &N_functions);

	// set amount of equations
	set_n_systems ();
	// set time layers
	set_time_layers (time_file_name);
	// prepares elements
	for (int i = 0, i_end = mesh_pointer->get_n_elements (); i < i_end; i++)
	{
		mesh_pointer->prepare_element (i);
	}

	// set conditions
	conditions = new Condition *[n_systems];
	for (int i = 0; i < n_systems; i++)
	{
		conditions[i] = new Condition[mesh_pointer->get_dimentionality () * 2];
		get_conditions (i, bound_file_name[i]);
	}

	use_time_mapping = true;

	// build matrices' portraits
	for (int k_system = 0; k_system < n_systems; k_system++)
		if (!build_portrait (k_system))
		{
			printf ("ERROR: portrait\n");
			prepared = false;
		}

	// if portrait building was successful, set starting conditions
	if (prepared)
		set_starting_conditions ();
}

template<class Type_Mesh>
void Task<Type_Mesh>::prepare (char * nodes_file_name, char * elements_file_name, char * bound_file_name)
{
	printf ("\tpreparing time layers and elements\n");

	log_file = fopen ("Result Files Extra//log.txt", "w");
	// build a mesh
	mesh_pointer = &mesh;
	prepared = mesh_pointer->build_Mesh (nodes_file_name, elements_file_name, &N_functions);

	// set amount of equations
	set_n_systems ();
	// set time layers
	set_time_layers ("");
	// prepares elements
	for (int i = 0, i_end = mesh_pointer->get_n_elements (); i < i_end; i++)
	{
		mesh_pointer->prepare_element (i);
	}

	// set conditions
	conditions = new Condition *[n_systems];
	for (int i = 0; i < n_systems; i++)
	{
		conditions[i] = new Condition[mesh_pointer->get_dimentionality () * 2];
		get_conditions (i, bound_file_name);
	}

	use_time_mapping = true;

	// build matrices' portraits
	for (int k_system = 0; k_system < n_systems; k_system++)
		if (!build_portrait (k_system))
		{
			printf ("ERROR: portrait\n");
			prepared = false;
		}

	// if portrait building was successful, set starting conditions
	if (prepared)
		set_starting_conditions ();
}

template<class Type_Mesh>
void Task<Type_Mesh>::prepare (char * nodes_file_name, char * elements_file_name, char * bound_file_name, char * time_file_name)
{
	printf ("\tpreparing time layers and elements\n");

	log_file = fopen ("Result Files Extra//log.txt", "w");
	// build a mesh
	mesh_pointer = &mesh;
	prepared = mesh_pointer->build_Mesh (nodes_file_name, elements_file_name, &N_functions);

	// set amount of equations
	set_n_systems ();
	// set time layers
	set_time_layers (time_file_name);
	// prepares elements
	for (int i = 0, i_end = mesh_pointer->get_n_elements (); i < i_end; i++)
	{
		mesh_pointer->prepare_element (i);
	}

	// set conditions
	conditions = new Condition *[n_systems];
	for (int i = 0; i < n_systems; i++)
	{
		conditions[i] = new Condition[mesh_pointer->get_dimentionality () * 2];
		get_conditions (i, bound_file_name);
	}

	use_time_mapping = true;

	// build matrices' portraits
	for (int k_system = 0; k_system < n_systems; k_system++)
		if (!build_portrait (k_system))
		{
			printf ("ERROR: portrait\n");
			prepared = false;
		}

	// if portrait building was successful, set starting conditions
	if (prepared)
		set_starting_conditions ();
}

template<class Type_Mesh>
void Task<Type_Mesh>::prepare (char * nodes_file_name, char * elements_file_name, char * bound_file_name[], char * time_file_name)
{
	printf ("\tpreparing time layers and elements\n");

	log_file = fopen ("Result Files Extra//log.txt", "w");
	// build a mesh
	mesh_pointer = &mesh;
	prepared = mesh_pointer->build_Mesh (nodes_file_name, elements_file_name, &N_functions);

	// set amount of equations
	set_n_systems ();
	// set time layers
	set_time_layers (time_file_name);
	// prepares elements
	for (int i = 0, i_end = mesh_pointer->get_n_elements (); i < i_end; i++)
	{
		mesh_pointer->prepare_element (i);
	}

	// set conditions
	conditions = new Condition *[n_systems];
	for (int i = 0; i < n_systems; i++)
	{
		conditions[i] = new Condition[mesh_pointer->get_dimentionality () * 2];
		get_conditions (i, bound_file_name[i]);
	}

	use_time_mapping = true;

	// build matrices' portraits
	for (int k_system = 0; k_system < n_systems; k_system++)
		if (!build_portrait (k_system))
		{
			printf ("ERROR: portrait\n");
			prepared = false;
		}

	// if portrait building was successful, set starting conditions
	if (prepared)
		set_starting_conditions ();

}

template<class Type_Mesh>
void Task<Type_Mesh>::read_materials (char * materials_file_name)
{
	FILE * file = fopen (materials_file_name, "r");

	lambda_val.clear ();
	gamma_val.clear ();

	int n;
	int id;
	double l, g;
	fscanf (file, "%i", &n);
	for (int i = 0; i < n; i++)
	{
		fscanf (file, "%i %lf %lf", &id, &l, &g);
		lambda_val.insert (std::pair<int, double>(id, l));
		gamma_val.insert (std::pair<int, double>(id, g));
	}
	fclose (file);
}

template<class Type_Mesh>
void Task<Type_Mesh>::set_time_layers (char * file_name)
{
	n_time_layers = 2; // default value
	time_layers = new double[n_time_layers];
	time_layers[0] = 0.0; // "starting condition"
	time_layers[1] = 0.0; // first-level solution

	time_sampling = 1; // default value
	previous_time_layers_solutions = new MathVector *[n_systems]; // solutions are saved for each system
	for (int k_systems = 0; k_systems < n_systems; k_systems++)
	{
		previous_time_layers_solutions[k_systems] = new MathVector[time_sampling];
	}

	printf ("time layers:\t%i\ntime sampling:\t%i \n", n_time_layers, time_sampling);
	fprintf (log_file, "time layers:\t%i\ntime sampling:\t%i \n", n_time_layers, time_sampling);
}

template<class Type_Mesh>
void Task<Type_Mesh>::set_n_systems ()
{
	n_systems = 1; // default value
	eq_matrixes = new compressed_matrix[n_systems];
	printf ("systems settled:\t%i\n", n_systems);
	fprintf (log_file, "systems settled:\t%i\n", n_systems);
}

template<class Type_Mesh>
void Task<Type_Mesh>::set_starting_conditions ()
{
	for (int k_system = 0; k_system < n_systems; k_system++)
	{
		for (int i = 0; i < time_sampling; i++)
		{
			previous_time_layers_solutions[k_system][i].setSize (eq_matrixes[k_system].Size ()); // default starting solution point is 0
		}
	}

	// make space for possible derivatives
	int dim = mesh_pointer->get_dimentionality ();
	derivative = new MathVector[dim]; // solutions are saved for each system
	for (int i = 0; i < dim; i++)
		derivative[i].setSize (eq_matrixes[0].Size ());

	printf ("\tstarting conditions are settled\n");
}

template<class Type_Mesh>
void Task<Type_Mesh>::solution_step (int k_system, int method, int d_type, int depth)
{
	int solver_iterations;

	// for linear system
	eq_matrixes[k_system].Clear ();
	eq_matrixes[k_system].set_starting_point (&(previous_time_layers_solutions[k_system][0])); // set last layer solution as x0 for solver
	build_system (k_system);  // reset system for current layer
	apply_boundary_conditions (k_system); //apply boundary conditions
	//eq_matrixes[k_system].fprint ();

	int solver_param[] = { method , d_type, depth, 1, SOLVER_MKL_NO };
	solver_iterations = eq_matrixes[k_system].solve (solver_param);

	printf ("system: %i, time layer:\t%lf, iter:\t%i\n", k_system, time_layers[current_time_layer], solver_iterations);
	fprintf (log_file, "system: %i, time layer: %lf, iter: %i\n", k_system, time_layers[current_time_layer], solver_iterations);
	next_time_layer (k_system); // save solution, move up time layer
}

template<class Type_Mesh>
void Task<Type_Mesh>::solution_step (int k_system, int * solver_param)
{
	int solver_iterations;

	// for linear system
	eq_matrixes[k_system].Clear ();
	eq_matrixes[k_system].set_starting_point (&(previous_time_layers_solutions[k_system][0])); // set last layer solution as x0 for solver
	build_system (k_system);  // reset system for current layer
	apply_boundary_conditions (k_system); //apply boundary conditions
										  //eq_matrixes[k_system].fprint ();

	solver_iterations = eq_matrixes[k_system].solve (solver_param);

	printf ("system: %i, time layer:\t%lf, iter:\t%i\n", k_system, time_layers[current_time_layer], solver_iterations);
	fprintf (log_file, "system: %i, time layer: %lf, iter: %i\n", k_system, time_layers[current_time_layer], solver_iterations);
	next_time_layer (k_system); // save solution, move up time layer
}

template<class Type_Mesh>
void Task<Type_Mesh>::next_time_layer (int k_system)
{
	// TEST
	for (int j = time_sampling - 1; j > 0; j--)
	{
		previous_time_layers_solutions[k_system][j] = previous_time_layers_solutions[k_system][j - 1];
	}

	eq_matrixes[k_system].get_solution (&(previous_time_layers_solutions[k_system][0]));
}

template<class Type_Mesh>
void Task<Type_Mesh>::next_time_layer (int k_system, MathVector * solution)
{
	for (int j = time_sampling - 1; j > 0; j--)
	{
		previous_time_layers_solutions[k_system][j].Copy ((previous_time_layers_solutions[k_system][j - 1]));
	}

	previous_time_layers_solutions[k_system][0].Copy (*solution);
}

template<class Type_Mesh>
void Task<Type_Mesh>::print_solutions ()
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
		for (int k = 0; k < n_systems; k++)
		{
			for (size_t i = 0, i_end = local_time_stamps.size (); i < i_end; i++)
			{
				get_time_approx (local_time_stamps[i], c);

				printf ("\tsaving solutions into files, time layer:\t%.16lf\n", local_time_stamps[i]);
				sprintf (name, "Result//s%i_t_%0.5lf.txt", k, local_time_stamps[i]);
				file = fopen (name, "w");
				//previous_time_layers_solutions[k][0].FPrint (file);

				double val;
				for (int j = 0, j_end = previous_time_layers_solutions[k][0].getSize (); j < j_end; j++)
				{
					val = 0;
					for (int t = 0; t < time_sampling; t++)
					{
						val += c->getElem (t) * previous_time_layers_solutions[k][t].getElem (j);
					}
					fprintf (file, "%.16lf\n", val);
				}
				fclose (file);

				if (false)
				{
					wchar_t pic_name[128];
					{
						swprintf (pic_name, 128, L"Pictures//s%i_t_%lf_iso.png", k, local_time_stamps[i]);
						std::wstring pic_name_iso (pic_name);
						//painter->draw_isovalues (k, pic_name_iso);

						swprintf (pic_name, 128, L"Pictures//s%i_t_%lf_field.png", k, local_time_stamps[i]);
						std::wstring pic_name_field (pic_name);
						//painter->draw_field (k, pic_name_field);
					}
					
					{
						swprintf (pic_name, 128, L"Pictures//s%i_t_%lf.png", k, local_time_stamps[i]);
						std::wstring pic_name_field (pic_name);
						//painter->draw_field_and_isovalues (k, pic_name_field);
					}
				}
				fclose (file);
			}
		}
		delete c;
	}
}

template<class Type_Mesh>
void Task<Type_Mesh>::print_extra_data ()
{
}

template<class Type_Mesh>
void Task<Type_Mesh>::get_time_coefficients (MathVector * c)
{
	c->setElem (0, 1.0);
	c->setElem (1, 0.0);
}

template<class Type_Mesh>
void Task<Type_Mesh>::get_time_approx (double t, MathVector * c)
{
	c->setElem (0, 1.0);
}

template<class Type_Mesh>
bool Task<Type_Mesh>::get_solution_in_point (int k_system, double * coordinates, double * value)
{
	double * c = new double[mesh_pointer->get_dimentionality ()];
	for (int i = 0; i < mesh_pointer->get_dimentionality (); i++)
		c[i] = coordinates[i];

	bool sym_point = false;
	// find element that contains the point 
	int p_element = mesh_pointer->point_inside (coordinates);
	if (p_element == -1 && symmetrical != 0)
	{	
		c[0] = -c[0];
		p_element = mesh_pointer->point_inside (c);
		sym_point = true;
	}
	if (p_element == -1)
	{
		*value = 0.0;
		return false;
	}

	// get local functions' values of the element
	int func_amount = mesh_pointer->get_amount_non_zero_functions (p_element);
	MathVector * local_func = new MathVector (func_amount);
	mesh_pointer->get_local_function_values (p_element, c, local_func);

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
	*value = local_func->Scalar_Product (*solution);
	if (sym_point && symmetrical == 2)
		*value *= -1.0;

	delete local_func;
	delete solution;
	delete[] c;

	return true;
}

template<class Type_Mesh>
bool Task<Type_Mesh>::get_solution_in_point (int k_system, int k_time_layer, double * coordinates, double * value)
{
	// find element that contains the point 
	int p_element = mesh_pointer->point_inside (coordinates);
	if (p_element == -1)
	{
		*value = 0.0;
		return false;
	}

	// get local functions' values of the element
	int func_amount = mesh_pointer->get_amount_non_zero_functions (p_element);
	MathVector * local_func = new MathVector (func_amount);
	mesh_pointer->get_local_function_values (p_element, coordinates, local_func);

	// make vector of solutions in those nodes
	MathVector * solution = new MathVector (func_amount);
	int global_function;
	for (int i = 0; i < func_amount; i++)
	{
		global_function = get_function_global_number (p_element, i);
		if (global_function != -1)
		{
			solution->setElem (i, previous_time_layers_solutions[k_system][k_time_layer].getElem (global_function));
		}
	}

	// multiply them by the solution's values in those nodes
	*value = local_func->Scalar_Product (*solution);

	delete local_func;
	delete solution;

	return true;
}

template<class Type_Mesh>
bool Task<Type_Mesh>::get_non_linear_solution_in_point (int k_system, double * coordinates, double * value)
{
	double * c = new double[mesh_pointer->get_dimentionality ()];
	for (int i = 0; i < mesh_pointer->get_dimentionality (); i++)
		c[i] = coordinates[i];

	bool sym_point = false;
	// find element that contains the point 
	int p_element = mesh_pointer->point_inside (coordinates);
	if (p_element == -1 && symmetrical != 0)
	{
		c[0] = -c[0];
		p_element = mesh_pointer->point_inside (c);
		sym_point = true;
	}
	if (p_element == -1)
	{
		*value = 0.0;
		return false;
	}

	// get local functions' values of the element
	int func_amount = mesh_pointer->get_amount_non_zero_functions (p_element);
	MathVector * local_func = new MathVector (func_amount);
	mesh_pointer->get_local_function_values (p_element, c, local_func);

	// make vector of solutions in those nodes
	MathVector * solution = new MathVector (func_amount);
	int global_function;
	for (int i = 0; i < func_amount; i++)
	{
		global_function = get_function_global_number (p_element, i);
		if (global_function != -1)
		{
			solution->setElem (i, non_linear_layers_solutions[k_system].getElem (global_function));
		}
	}

	// multiply them by the solution's values in those nodes
	*value = local_func->Scalar_Product (*solution);
	if (sym_point && symmetrical == 2)
		*value *= -1.0;

	delete local_func;
	delete solution;
	delete[] c;

	return true;
}

template<class Type_Mesh>
bool Task<Type_Mesh>::get_derivative (int k_system, int k_var, double * coordinates, double *value)
{
	if (SR)
	{
		*value = get_derivative_SR (k_system, k_var, coordinates);
	}
	else
	{
		*value = get_derivative_first (k_system, k_var, coordinates);
	}
	return true;
}

template<class Type_Mesh>
void Task<Type_Mesh>::save_solution_in_points (int k_system, char * file_sour, char * file_dest, bool keep_points)
{
	int dim_node = mesh_pointer->get_dimentionality ();
	double * coordinates = new double[dim_node];
	double value;
	bool exists;
	int n_amount = 0;

	FILE * result_points = fopen (file_dest, "w");
	FILE * points = fopen (file_sour, "r");
	while (!feof (points))
	{
		value = 0.0;
		for (int i = 0; i < dim_node; i++)
		{
			fscanf (points, "%lf", &coordinates[i]);
		}
		exists = get_solution_in_point (k_system, coordinates, &value);
		if (keep_points)
		{
			if (exists)
			{
				for (int i = 0; i < dim_node; i++)
				{
					fprintf (result_points, "%.16lf\t", coordinates[i]);
				}
				fprintf (result_points, "%.16lf\n", value);
				n_amount++;
			}
		}
		else
		{
			fprintf (result_points, "%.16lf\n", value);
		}
	}

	fclose (points);
	fclose (result_points);

	delete[] coordinates;

	//FILE * file = fopen ("Result Files Extra\\points_amount.txt", "w");
	//fprintf (file, "%i", n_amount);
	//fclose (file);

	//printf ("\trequested solution is saved\n");
}

template<class Type_Mesh>
void Task<Type_Mesh>::save_solution_in_points_first_derivative (int k_system, char * file_sour, char * file_dest, bool keep_points, int k_var, double h)
{
	int dim_node = mesh_pointer->get_dimentionality ();
	double * coordinates = new double[dim_node];
	double value;
	int n_amount = 0;

	FILE * result_points = fopen (file_dest, "w");
	FILE * points = fopen (file_sour, "r");
	while (!feof (points))
	{
		value = 0.0;
		for (int i = 0; i < dim_node; i++)
		{
			fscanf (points, "%lf", &coordinates[i]);
		}
		value = get_derivative_first (k_system, k_var, coordinates);
		if (keep_points)
		{
			for (int i = 0; i < dim_node; i++)
			{
				fprintf (result_points, "%.16lf\t", coordinates[i]);
			}
			fprintf (result_points, "%.16lf\n", value);
			n_amount++;
		}
		else
		{
			fprintf (result_points, "%.16lf\n", value);
		}
	}

	fclose (points);
	fclose (result_points);

	delete[] coordinates;

	//printf ("\trequested derivative is saved\n");
}

template<class Type_Mesh>
void Task<Type_Mesh>::save_solution_in_points_second_derivative (int k_system, char * file_sour, char * file_dest, bool keep_points, int k_var1, int k_var2, double h1, double h2)
{
	int dim_node = mesh_pointer->get_dimentionality ();
	double * coordinates = new double[dim_node];
	double value;
	int n_amount = 0;

	FILE * result_points = fopen (file_dest, "w");
	FILE * points = fopen (file_sour, "r");
	while (!feof (points))
	{
		value = 0.0;
		for (int i = 0; i < dim_node; i++)
		{
			fscanf (points, "%lf", &coordinates[i]);
		}
		value = get_derivative_second (k_system, k_var1, k_var2, h1, h2, coordinates);
		if (keep_points)
		{
			for (int i = 0; i < dim_node; i++)
			{
				fprintf (result_points, "%.16lf\t", coordinates[i]);
			}
			fprintf (result_points, "%.16lf\n", value);
			n_amount++;
		}
		else
		{
			fprintf (result_points, "%.16lf\n", value);
		}
	}

	fclose (points);
	fclose (result_points);

	delete[] coordinates;

	//printf ("\trequested second derivative is saved\n");
}

template<class Type_Mesh>
void Task<Type_Mesh>::save_slice (int k_system, char * file_sour, char * file_dest, bool keep_points)
{
	int dim_node = mesh_pointer->get_dimentionality ();
	double * coordinates = new double[dim_node];
	double * coordinatesN = new double[dim_node];
	double value; // value of the solution
	bool exists; // flag of existence of solution's value 
	int n_amount = 0; // amount of points written in file 
	int k_var; // variable to iterate through
	double h; // iteration step
	int iter_amount; // amount of points to iterate through
	double vstart, vend; // starting and ending values of variable 

	FILE * points = fopen (file_sour, "r");
	fscanf (points, "%i %lf", &k_var, &h); // get number of the variable to iterate through and step for it
	if (k_var >= 0 && k_var < dim_node) // if the variable number is valide
	{
		// get other variables' values (for slice)
		for (int i = 0; i < dim_node; i++)
		{
			if (i != k_var)
			{
				fscanf (points, "%lf", &coordinates[i]);
			}
		}
		// get starting value
		mesh_pointer->get_0_boundaries (coordinatesN);
		vstart = coordinatesN[k_var];
		// get ending value
		mesh_pointer->get_N_boundaries (coordinatesN);
		vend = coordinatesN[k_var];

		// for syymetrical tasks
		if ((k_var == 0) && is_symmetrical () != 0)
			vstart = -coordinatesN[k_var];
		
		// get amount of points to iterate through
		iter_amount = (int)((vend - vstart) / h);

		// open file for result input
		FILE * result_points = fopen (file_dest, "w");
		for (int i = 0; i < iter_amount + 1; i++)
		{
			coordinates[k_var] = vstart + i * h;
			//calculate the point
			value = 0.0;
			// get value in the point
			if ( (is_symmetrical() == 0) || (is_symmetrical() != 0 && k_var != 0) )
			{
				exists = get_solution_in_point (k_system, coordinates, &value);
			}
			else
			{
				coordinatesN[k_var] = -coordinates[k_var];
				exists = get_solution_in_point (k_system, coordinatesN, &value);
				if (is_symmetrical () == 2)
				{
					value = -value;
				}
			}
			// if value has to be saved with the point
			if (keep_points)
			{
				// then print in out only of it exists
				if (exists)
				{
					//for (int i = 0; i < dim_node; i++)
					//{
					//	fprintf (result_points, "%.16lf ", coordinates[i]);
					//}
					fprintf (result_points, "%.16lf\t", coordinates[k_var]);

					fprintf (result_points, "%.16lf\n", value);
					n_amount++;
				}
			}
			else // otherwise print it anyway, it'll be zero
			{
				fprintf (result_points, "%.16lf\n", value);
				n_amount++;
			}
		}
		fclose (result_points);
	}

	fclose (points);

	delete[] coordinates;
	delete[] coordinatesN;

	//FILE * file = fopen ("Result Files Extra\\points_amount.txt", "w");
	//fprintf (file, "%i", n_amount);
	//fclose (file);

	//printf ("\trequested solution is saved\n");
}

template<class Type_Mesh>
void Task<Type_Mesh>::save_slice (int k_system, double * boundaries, char * file_sour, char * file_dest, bool keep_points)
{
	int dim_node = get_mesh_dim ();
	double * coordinates = new double[dim_node];
	double * coordinatesN = new double[dim_node];
	double value; // value of the solution
	bool exists; // flag of existence of solution's value 
	int n_amount = 0; // amount of points written in file 
	int k_var; // variable to iterate through
	double h; // iteration step
	int iter_amount; // amount of points to iterate through
	double vstart, vend; // starting and ending values of variable 

	FILE * points = fopen (file_sour, "r");
	fscanf (points, "%i %lf", &k_var, &h); // get number of the variable to iterate through and step for it
	if (k_var >= 0 && k_var < dim_node) // if the variable number is valide
	{
		// get other variables' values (for slice)
		for (int i = 0; i < dim_node; i++)
		{
			if (i != k_var)
			{
				fscanf (points, "%lf", &coordinates[i]);
			}
		}
		// get starting value
		vstart = boundaries[0];
		// get ending value
		vend = boundaries[1];

		// for syymetrical tasks
		if ((k_var == 0) && is_symmetrical () != 0)
			vstart = -boundaries[0];

		// get amount of points to iterate through
		iter_amount = (int)((vend - vstart) / h);

		// open file for result input
		FILE * result_points = fopen (file_dest, "w");
		for (int i = 0; i < iter_amount + 1; i++)
		{
			coordinates[k_var] = vstart + i * h;
			//calculate the point
			value = 0.0;
			// get value in the point
			if ((is_symmetrical () == 0) || (is_symmetrical () != 0 && k_var != 0))
			{
				exists = get_solution_in_point (k_system, coordinates, &value);
			}
			else
			{
				coordinatesN[k_var] = -coordinates[k_var];
				exists = get_solution_in_point (k_system, coordinatesN, &value);
				if (is_symmetrical () == 2)
				{
					value = -value;
				}
			}
			// if value has to be saved with the point
			if (keep_points)
			{
				// then print in out only of it exists
				if (exists)
				{
					//for (int i = 0; i < dim_node; i++)
					//{
					//	fprintf (result_points, "%.16lf ", coordinates[i]);
					//}
					fprintf (result_points, "%.16lf\t", coordinates[k_var]);

					fprintf (result_points, "%.16lf\n", value);
					n_amount++;
				}
			}
			else // otherwise print it anyway, it'll be zero
			{
				fprintf (result_points, "%.16lf\n", value);
				n_amount++;
			}
		}
		fclose (result_points);
	}

	fclose (points);

	delete[] coordinates;
	delete[] coordinatesN;

	//FILE * file = fopen ("Result Files Extra\\points_amount.txt", "w");
	//fprintf (file, "%i", n_amount);
	//fclose (file);

	//printf ("\trequested solution is saved\n");
}

template<class Type_Mesh>
void Task<Type_Mesh>::save_slice (int k_system, double seek_value, int k_var, char * file_dest, bool keep_points)
{
	int dim_node = mesh_pointer->get_dimentionality ();
	double * coordinates = new double[dim_node];
	double * coordinatesN = new double[dim_node];
	double value; // value of the solution
	int n_amount = 0; // amount of points written in file 
	double h = 1e-3; // iteration step
	int iter_amount; // amount of points to iterate through
	double vstart, vend; // starting and ending values of variable 
	double * sol_point = new double[dim_node];
	double * c0 = new double[dim_node];
	double * cN = new double[dim_node];
	double margin;

	if (k_var >= 0 && k_var < dim_node) // if the variable number is valid
	{
		mesh_pointer->get_0_boundaries (coordinatesN);
		for (int i = 0; i < dim_node; i++)
		{
			c0[i] = coordinatesN[i];
		}
		mesh_pointer->get_N_boundaries (coordinatesN);
		for (int i = 0; i < dim_node; i++)
		{
			cN[i] = coordinatesN[i];
		}
		// get starting value
		vstart = c0[k_var];
		// get ending value
		vend = cN[k_var];

		// for syymetrical tasks
		if ((k_var == 0) && is_symmetrical () != 0)
			vstart = -coordinatesN[k_var];

		// get amount of points to iterate through
		iter_amount = (int)((vend - vstart) / h);
		
		//c0[1] = 0.9; // for upper wall
		c0[1] = 0.05;
		cN[1] = 0.4;
		// open file for result input
		FILE * result_points = fopen (file_dest, "w");
		for (int i = 0; i < iter_amount + 1; i++)
		{
			coordinates[k_var] = vstart + i * h;
			//calculate the point
			value = 0.0;
			// get value in the point
			c0[k_var] = cN[k_var] = coordinates[k_var];
			if ((is_symmetrical () == 0) || (is_symmetrical () != 0 && k_var != 0))
			{
				margin = optimization::closest_value_point_GS (this, k_system, c0, cN, seek_value, sol_point);
			}
			else
			{
				c0[k_var] = cN[k_var] = -coordinates[k_var];
				margin = optimization::closest_value_point_GS (this, k_system, c0, cN, seek_value, sol_point);
			}
			// if value has to be saved with the point
			if (keep_points)
			{
				// then print in out only of it exists
				if (fabs(margin) < 9e-4)
				{
					//for (int i = 0; i < dim_node; i++)
					//{
					//	fprintf (result_points, "%.16lf ", coordinates[i]);
					//}
					//fprintf (result_points, "%.e ", margin);
					fprintf (result_points, "%.16lf\t", c0[k_var]);

					fprintf (result_points, "%.16lf\n", sol_point[1]);
					n_amount++;
				}
			}
			else // otherwise print it anyway, it'll be zero
			{
				fprintf (result_points, "%.16lf\n", value);
				n_amount++;
			}
		}
		fclose (result_points);
	}
	
	delete[] coordinates;
	delete[] coordinatesN;
	delete[] sol_point;
	delete[] c0;
	delete[] cN;

	//FILE * file = fopen ("Result Files Extra\\points_amount.txt", "w");
	//fprintf (file, "%i", n_amount);
	//fclose (file);

	//printf ("\trequested solution is saved\n");
}

template<class Type_Mesh>
void Task<Type_Mesh>::save_slice_weighted (int k_system, double * weights, char * file_sour, char * file_dest, bool keep_points)
{
	int dim_node = mesh_pointer->get_dimentionality ();
	double * coordinates = new double[dim_node];
	double * coordinatesN = new double[dim_node];
	double value; // value of the solution
	bool exists; // flag of existence of solution's value 
	int n_amount = 0; // amount of points written in file 
	int k_var; // variable to iterate through
	double h; // iteration step
	int iter_amount; // amount of points to iterate through
	double vstart, vend; // starting and ending values of variable 

	FILE * points = fopen (file_sour, "r");
	fscanf (points, "%i %lf", &k_var, &h); // get number of the variable to iterate through and step for it
	if (k_var >= 0 && k_var < dim_node) // if the variable number is valide
	{
		// get other variables' values (for slice)
		for (int i = 0; i < dim_node; i++)
		{
			if (i != k_var)
			{
				fscanf (points, "%lf", &coordinates[i]);
			}
		}
		// get starting value
		mesh_pointer->get_0_boundaries (coordinatesN);
		vstart = coordinatesN[k_var];
		// get ending value
		mesh_pointer->get_N_boundaries (coordinatesN);
		vend = coordinatesN[k_var];

		// get amount of points to iterate through
		iter_amount = (int)((vend - vstart) / h);

		// open file for result input
		FILE * result_points = fopen (file_dest, "w");


		if (is_symmetrical () != 0)
		{
			for (int i = 0; i < iter_amount; i++)
			{
				coordinates[k_var] = coordinatesN[k_var] - i * h;
				//calculate the point
				// get value in the point
				exists = get_solution_in_point (k_system, coordinates, &value);
				if (is_symmetrical () == 2)
					value = -value;

				// if value has to be saved with the point
				if (exists)
				{
					if (keep_points)
					{
						//for (int i = 0; i < dim_node; i++)
						//{
						//	fprintf (result_points, "%.16lf ", coordinates[i]);
						//}

						fprintf (result_points, "%.16lf\t", -coordinates[k_var] * weights[0]);
						fprintf (result_points, "%.16lf\n", value * weights[1] + weights[2]);
						n_amount++;
					}
					else // otherwise print it anyway, it'll be zero
					{
						fprintf (result_points, "%.16lf\n", value * weights[1] + weights[2]);
						n_amount++;
					}
				}
			}
		}

		for (int i = 0; i < iter_amount + 1; i++)
		{
			coordinates[k_var] = vstart + i * h;
			//calculate the point
			value = 0.0;
			// get value in the point
			exists = get_solution_in_point (k_system, coordinates, &value);
			// if value has to be saved with the point
			if (keep_points)
			{
				// then print in out only of it exists
				if (exists)
				{
					//for (int i = 0; i < dim_node; i++)
					//{
					//	fprintf (result_points, "%.16lf ", coordinates[i]);
					//}
					fprintf (result_points, "%.16lf\t", coordinates[k_var] * weights[0]);

					fprintf (result_points, "%.16lf\n", value * weights[1] + weights[2]);
					n_amount++;
				}
			}
			else // otherwise print it anyway, it'll be zero
			{
				fprintf (result_points, "%.16lf\n", value * weights[1] + weights[2]);
				n_amount++;
			}
		}
		fclose (result_points);
	}

	fclose (points);

	delete[] coordinates;
	delete[] coordinatesN;

	//FILE * file = fopen ("Result Files Extra\\points_amount.txt", "w");
	//fprintf (file, "%i", n_amount);
	//fclose (file);

	//printf ("\trequested solution is saved\n");
}

template<class Type_Mesh>
void Task<Type_Mesh>::save_slice_first_derivative (int k_system, char * file_sour, char * file_dest, bool keep_points)
{
	int dim_node = mesh_pointer->get_dimentionality ();
	double * coordinates = new double[dim_node];
	double * coordinatesN = new double[dim_node];
	double value; // value of the solution
	int n_amount = 0; // amount of points written in file 
	int k_var; // variable to iterate through
	double h; // iteration step
	int iter_amount; // amount of points to iterate through
	double vstart, vend; // starting and ending values of variable 

	FILE * points = fopen (file_sour, "r");
	fscanf (points, "%i %lf", &k_var, &h); // get number of the variable to iterate through and step for it
	if (k_var >= 0 && k_var < dim_node) // if the variable number is valide
	{
		// get other variables' values (for slice)
		for (int i = 0; i < dim_node; i++)
		{
			if (i != k_var)
			{
				fscanf (points, "%lf", &coordinates[i]);
			}
		}
		// get starting value
		mesh_pointer->get_0_boundaries (coordinatesN);
		vstart = coordinatesN[k_var];
		// get ending value
		mesh_pointer->get_N_boundaries (coordinatesN);
		vend = coordinatesN[k_var];
		// get amount of points to iterate through
		iter_amount = (int)((vend - vstart) / h);

		// open file for result input
		FILE * result_points = fopen (file_dest, "w");
		for (int i = 0; i < iter_amount + 1; i++)
		{
			coordinates[k_var] = vstart + i * h;
			//calculate the point
			// get value in the point
			bool found = get_derivative (k_system, 1, coordinates, &value);
			// if value has to be saved with the point
			if (keep_points)
			{
				//for (int i = 0; i < dim_node; i++)
				//{
				//	fprintf (result_points, "%.16lf ", coordinates[i]);
				//}

				fprintf (result_points, "%.16lf\t", coordinates[k_var]);
				fprintf (result_points, "%.16lf\n", value);
				n_amount++;
			}
			else // otherwise print it anyway, it'll be zero
			{
				fprintf (result_points, "%.16lf\n", value);
				n_amount++;
			}
		}
		fclose (result_points);
	}

	fclose (points);

	delete[] coordinates;
	delete[] coordinatesN;

	//printf ("\trequested solution is saved\n");
}

template<class Type_Mesh>
void Task<Type_Mesh>::save_slice_first_derivative (int k_system, int der_var, char * file_sour, char * file_dest, bool keep_points)
{
	int dim_node = mesh_pointer->get_dimentionality ();
	double * coordinates = new double[dim_node];
	double * coordinatesN = new double[dim_node];
	double value; // value of the solution
	bool exists; // flag of existence of solution's value 
	int n_amount = 0; // amount of points written in file 
	int k_var; // variable to iterate through
	double h; // iteration step
	int iter_amount; // amount of points to iterate through
	double vstart, vend; // starting and ending values of variable 

	FILE * points = fopen (file_sour, "r");
	fscanf (points, "%i %lf", &k_var, &h); // get number of the variable to iterate through and step for it
	if (k_var >= 0 && k_var < dim_node) // if the variable number is valide
	{
		// get other variables' values (for slice)
		for (int i = 0; i < dim_node; i++)
		{
			if (i != k_var)
			{
				fscanf (points, "%lf", &coordinates[i]);
			}
		}
		// get starting value
		mesh_pointer->get_0_boundaries (coordinatesN);
		vstart = coordinatesN[k_var];
		// get ending value
		mesh_pointer->get_N_boundaries (coordinatesN);
		vend = coordinatesN[k_var];

		// get amount of points to iterate through
		iter_amount = (int)((vend - vstart) / h);

		// open file for result input
		FILE * result_points = fopen (file_dest, "w");
		// for symmetrical tasks

		if ((k_var == 0) && is_symmetrical () != 0)
		{
			for (int i = 0; i < iter_amount; i++)
			{
				coordinates[k_var] = coordinatesN[k_var] - i * h;
				//calculate the point
				// get value in the point
				exists = get_derivative (k_system, der_var, coordinates, &value);
				if (is_symmetrical () == 2)
					value = -value;

				// if value has to be saved with the point
				if (exists)
				{
					if (keep_points)
					{
						//for (int i = 0; i < dim_node; i++)
						//{
						//	fprintf (result_points, "%.16lf ", coordinates[i]);
						//}

						fprintf (result_points, "%.16lf\t", -coordinates[k_var]);
						fprintf (result_points, "%.16lf\n", value);
						n_amount++;
					}
					else // otherwise print it anyway, it'll be zero
					{
						fprintf (result_points, "%.16lf\n", value);
						n_amount++;
					}
				}
			}
		}
			
			for (int i = 0; i < iter_amount + 1; i++)
		{
			coordinates[k_var] = vstart + i * h;
			//calculate the point
			// get value in the point
			exists = get_derivative (k_system, der_var, coordinates, &value);
			// if value has to be saved with the point
			if (exists)
			{
				if (keep_points)
				{
					//for (int i = 0; i < dim_node; i++)
					//{
					//	fprintf (result_points, "%.16lf ", coordinates[i]);
					//}

					fprintf (result_points, "%.16lf\t", coordinates[k_var]);
					fprintf (result_points, "%.16lf\n", value);
					n_amount++;
				}
				else // otherwise print it anyway, it'll be zero
				{
					fprintf (result_points, "%.16lf\n", value);
					n_amount++;
				}
			}
		}
		fclose (result_points);
	}

	fclose (points);

	delete[] coordinates;
	delete[] coordinatesN;

	//printf ("\trequested solution is saved\n");
}

template<class Type_Mesh>
void Task<Type_Mesh>::save_slice_first_derivative_weighted (int k_system, int der_var, double * weights, char * file_sour, char * file_dest, bool keep_points)
{
	int dim_node = mesh_pointer->get_dimentionality ();
	double * coordinates = new double[dim_node];
	double * coordinatesN = new double[dim_node];
	double value; // value of the solution
	bool exists; // flag of existence of solution's value 
	int n_amount = 0; // amount of points written in file 
	int k_var; // variable to iterate through
	double h; // iteration step
	int iter_amount; // amount of points to iterate through
	double vstart, vend; // starting and ending values of variable 

	FILE * points = fopen (file_sour, "r");
	fscanf (points, "%i %lf", &k_var, &h); // get number of the variable to iterate through and step for it
	if (k_var >= 0 && k_var < dim_node) // if the variable number is valide
	{
		// get other variables' values (for slice)
		for (int i = 0; i < dim_node; i++)
		{
			if (i != k_var)
			{
				fscanf (points, "%lf", &coordinates[i]);
			}
		}
		// get starting value
		mesh_pointer->get_0_boundaries (coordinatesN);
		vstart = coordinatesN[k_var];
		// get ending value
		mesh_pointer->get_N_boundaries (coordinatesN);
		vend = coordinatesN[k_var];

		// get amount of points to iterate through
		iter_amount = (int)((vend - vstart) / h);

		// open file for result input
		FILE * result_points = fopen (file_dest, "w");
		// for symmetrical tasks

		if ((k_var == 0) && is_symmetrical () != 0)
		{
			for (int i = 0; i < iter_amount; i++)
			{
				coordinates[k_var] = coordinatesN[k_var] - i * h;
				//calculate the point
				// get value in the point
				exists = get_derivative (k_system, der_var, coordinates, &value);
				if (is_symmetrical () == 2)
					value = -value;
				// if value has to be saved with the point
				if (exists)
				{
					if (keep_points)
					{
						//for (int i = 0; i < dim_node; i++)
						//{
						//	fprintf (result_points, "%.16lf ", coordinates[i]);
						//}

						fprintf (result_points, "%.16lf\t", -coordinates[k_var] * weights[0]);
						fprintf (result_points, "%.16lf\n", value * weights[1]);
						n_amount++;
					}
					else // otherwise print it anyway, it'll be zero
					{
						fprintf (result_points, "%.16lf\n", value * weights[1]);
						n_amount++;
					}
				}
			}
		}

		for (int i = 0; i < iter_amount + 1; i++)
		{
			coordinates[k_var] = vstart + i * h;
			//calculate the point
			// get value in the point
			exists = get_derivative (k_system, der_var, coordinates, &value);
			// if value has to be saved with the point
			if (exists)
			{
				if (keep_points)
				{
					//for (int i = 0; i < dim_node; i++)
					//{
					//	fprintf (result_points, "%.16lf ", coordinates[i]);
					//}

					fprintf (result_points, "%.16lf\t", coordinates[k_var] * weights[0]);
					fprintf (result_points, "%.16lf\n", value * weights[1]);
					n_amount++;
				}
				else // otherwise print it anyway, it'll be zero
				{
					fprintf (result_points, "%.16lf\n", value * weights[1]);
					n_amount++;
				}
			}
		}
		fclose (result_points);
	}

	fclose (points);

	delete[] coordinates;
	delete[] coordinatesN;

	//printf ("\trequested solution is saved\n");
}

template<class Type_Mesh>
void Task<Type_Mesh>::save_slice_first_derivative_absolute (int k_system, int der_var, double * weights, char * file_sour, char * file_dest, bool keep_points)
{
	int dim_node = mesh_pointer->get_dimentionality ();
	double * coordinates = new double[dim_node];
	double * coordinatesN = new double[dim_node];
	double value; // value of the solution
	bool exists; // flag of existence of solution's value 
	int n_amount = 0; // amount of points written in file 
	int k_var; // variable to iterate through
	double h; // iteration step
	int iter_amount; // amount of points to iterate through
	double vstart, vend; // starting and ending values of variable 

	FILE * points = fopen (file_sour, "r");
	fscanf (points, "%i %lf", &k_var, &h); // get number of the variable to iterate through and step for it
	if (k_var >= 0 && k_var < dim_node) // if the variable number is valide
	{
		// get other variables' values (for slice)
		for (int i = 0; i < dim_node; i++)
		{
			if (i != k_var)
			{
				fscanf (points, "%lf", &coordinates[i]);
			}
		}
		// get starting value
		mesh_pointer->get_0_boundaries (coordinatesN);
		vstart = coordinatesN[k_var];
		// get ending value
		mesh_pointer->get_N_boundaries (coordinatesN);
		vend = coordinatesN[k_var];

		// get amount of points to iterate through
		iter_amount = (int)((vend - vstart) / h);

		// open file for result input
		FILE * result_points = fopen (file_dest, "w");
		// for symmetrical tasks

		if ((k_var == 0) && is_symmetrical () != 0)
		{
			for (int i = 0; i < iter_amount; i++)
			{
				coordinates[k_var] = coordinatesN[k_var] - i * h;
				//calculate the point
				// get value in the point
				exists = get_derivative (k_system, der_var, coordinates, &value);
				if (is_symmetrical () == 2)
					value = -value;
				// if value has to be saved with the point
				if (exists)
				{
					if (keep_points)
					{
						//for (int i = 0; i < dim_node; i++)
						//{
						//	fprintf (result_points, "%.16lf ", coordinates[i]);
						//}

						fprintf (result_points, "%.16lf\t", -coordinates[k_var] * weights[0]);
						fprintf (result_points, "%.16lf\n", value * weights[1]);
						n_amount++;
					}
					else // otherwise print it anyway, it'll be zero
					{
						fprintf (result_points, "%.16lf\n", value * weights[1]);
						n_amount++;
					}
				}
			}
		}

		for (int i = 0; i < iter_amount + 1; i++)
		{
			coordinates[k_var] = vstart + i * h;
			//calculate the point
			// get value in the point
			exists = get_derivative (k_system, der_var, coordinates, &value);
			// if value has to be saved with the point
			if (exists)
			{
				if (keep_points)
				{
					//for (int i = 0; i < dim_node; i++)
					//{
					//	fprintf (result_points, "%.16lf ", coordinates[i]);
					//}

					fprintf (result_points, "%.16lf\t", coordinates[k_var] * weights[0]);
					fprintf (result_points, "%.16lf\n", fabs (value * weights[1]));
					n_amount++;
				}
				else // otherwise print it anyway, it'll be zero
				{
					fprintf (result_points, "%.16lf\n", fabs (value * weights[1]));
					n_amount++;
				}
			}
		}
		fclose (result_points);
	}

	fclose (points);

	delete[] coordinates;
	delete[] coordinatesN;

	//printf ("\trequested solution is saved\n");
}

template<class Type_Mesh>
void Task<Type_Mesh>::save_slice_second_derivative (int k_system, char * file_sour, char * file_dest, bool keep_points, int k_var2, double h1, double h2)
{
	int dim_node = mesh_pointer->get_dimentionality ();
	double * coordinates = new double[dim_node];
	double * coordinatesN = new double[dim_node];
	double value; // value of the solution
	int n_amount = 0; // amount of points written in file 
	int k_var; // variable to iterate through
	double h; // iteration step
	int iter_amount; // amount of points to iterate through
	double vstart, vend; // starting and ending values of variable 

	FILE * points = fopen (file_sour, "r");
	fscanf (points, "%i %lf", &k_var, &h); // get number of the variable to iterate through and step for it
	if (k_var >= 0 && k_var < dim_node) // if the variable number is valide
	{
		// get other variables' values (for slice)
		for (int i = 0; i < dim_node; i++)
		{
			if (i != k_var)
			{
				fscanf (points, "%lf", &coordinates[i]);
			}
		}
		// get starting value
		mesh_pointer->get_0_boundaries (coordinatesN);
		vstart = coordinatesN[k_var];
		// get ending value
		mesh_pointer->get_N_boundaries (coordinatesN);
		vend = coordinatesN[k_var];
		// get amount of points to iterate through
		iter_amount = (int)((vend - vstart) / h);

		// open file for result input
		FILE * result_points = fopen (file_dest, "w");
		for (int i = 0; i < iter_amount + 1; i++)
		{
			coordinates[k_var] = vstart + i * h;
			//calculate the point
			// get value in the point
			value = get_derivative_second (k_system, k_var2, k_var, h1, h2, coordinates);
			// if value has to be saved with the point
			if (keep_points)
			{
				for (int i = 0; i < dim_node; i++)
				{
					fprintf (result_points, "%.16lf\t", coordinates[i]);
				}
				fprintf (result_points, "%.16lf\n", value);
				n_amount++;
			}
			else // otherwise print it anyway, it'll be zero
			{
				fprintf (result_points, "%.16lf\n", value);
				n_amount++;
			}
		}
		fclose (result_points);
	}

	fclose (points);

	delete[] coordinates;
	delete[] coordinatesN;

	//FILE * file = fopen ("Result Files Extra\\points_amount.txt", "w");
	//fprintf (file, "%i", n_amount);
	//fclose (file);

	//printf ("\trequested solution is saved\n");
}

template<class Type_Mesh>
void Task<Type_Mesh>::save_evenly_2D (int k_system, double * h, char * file_dest, bool keep_points)
{
	FILE * file = fopen (file_dest, "w");
	int dim = mesh_pointer->get_dimentionality ();
	double * coordinates0 = new double[dim];
	mesh_pointer->get_0_boundaries (coordinates0);
	double * coordinatesN = new double[dim];
	mesh_pointer->get_N_boundaries (coordinatesN);
	double * coordinates = new double[dim];
	bool exists = false;
	double value;
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
			exists = get_solution_in_point (k_system, coordinates, &value);
			if (keep_points)
			{
				if (exists)
				{
					for (int k = 0; k < dim; k++)
					{
						fprintf (file, "%.16lf\t", coordinates[k]);
					}
					fprintf (file, "%.16lf\n", value);

					n_amount++;
				}
			}
			else
			{
				fprintf (file, "%.16lf\n", value);
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

	//printf ("\trequested solution is saved\n");
}

template<class Type_Mesh>
void Task<Type_Mesh>::save_evenly_2D_first_derivative (int k_system, int k_var, double * h, char * file_dest, bool keep_points)
{
	FILE * file = fopen (file_dest, "w");
	int dim = mesh_pointer->get_dimentionality ();
	double * coordinates0 = new double[dim];
	mesh_pointer->get_0_boundaries (coordinates0);
	double * coordinatesN = new double[dim];
	mesh_pointer->get_N_boundaries (coordinatesN);
	double * coordinates = new double[dim];
	bool exists = false;
	double value;
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
			value = get_derivative_first (k_system, k_var, coordinates);
			if (keep_points)
			{
				for (int k = 0; k < dim; k++)
				{
					fprintf (file, "%.16lf\t", coordinates[k]);
				}
				fprintf (file, "%.16lf\n", value);
			}
			else
			{
				fprintf (file, "%.16lf\n", value);
			}
		}
	}
	fclose (file);

	delete[] coordinates;
	delete[] coordinates0;
	delete[] coordinatesN;

	//printf ("\trequested solution is saved\n");
}

template<class Type_Mesh>
void Task<Type_Mesh>::save_evenly_2D_first_derivative_SR (int k_system, int k_var, double * h, char * file_dest, bool keep_points)
{
	FILE * file = fopen (file_dest, "w");
	int dim = mesh_pointer->get_dimentionality ();
	double * coordinates0 = new double[dim];
	mesh_pointer->get_0_boundaries (coordinates0);
	double * coordinatesN = new double[dim];
	mesh_pointer->get_N_boundaries (coordinatesN);
	double * coordinates = new double[dim];
	bool exists = false;
	double value;
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
			value = get_derivative_SR (k_system, k_var, coordinates);
			if (keep_points)
			{
				for (int k = 0; k < dim; k++)
				{
					fprintf (file, "%.16lf\t", coordinates[k]);
				}
				fprintf (file, "%.16lf\n", value);
			}
			else
			{
				fprintf (file, "%.16lf\n", value);
			}
		}
	}
	fclose (file);

	delete[] coordinates;
	delete[] coordinates0;
	delete[] coordinatesN;

	//printf ("\trequested solution is saved\n");
}

template<class Type_Mesh>
void Task<Type_Mesh>::fprint_edges_by_nodes (char * file_name)
{
	FILE * file = fopen (file_name, "w");
	int e;
	for (int i = 0, i_end = mesh_pointer->get_n_nodes (); i < i_end; i++)
	{
		for (int j = i + 1, j_end = mesh_pointer->get_n_nodes (); j < j_end; j++)
			//for (int j_end = mesh_pointer->get_n_nodes (), j = j_end - 1; j > i; j--)
		{
			e = mesh_pointer->edges->get_edge_number (i, j);
			if (e != -1)
			{
				fprintf (file, "%i:\t%i %i\n", e, i, j);
			}
		}
	}

	fclose (file);
}

template<class Type_Mesh>
void Task<Type_Mesh>::fprint_edges (char * file_name)
{
	FILE * file = fopen (file_name, "w");
	int n1, n2;
	for (int i = 0, i_end = mesh_pointer->edges->get_n_entries (); i < i_end; i++)
	{
		mesh_pointer->edges->get_edge_nodes (i, &n1, &n2);
		fprintf (file, "%i:\t%i %i\n", i, n1, n2);
	}

	fclose (file);
}

template<class Type_Mesh>
void Task<Type_Mesh>::fprint_functions_by_elements (char * file_name)
{
	FILE * file = fopen (file_name, "w");
	int * functions;
	int n_functions;

	for (int i = 0, i_end = mesh_pointer->get_n_elements (); i < i_end; i++)
	{
		fprintf (file, "%i:\t", i);
		n_functions = mesh_pointer->get_amount_non_zero_functions (i);
		functions = new int[n_functions];
		get_element_functions (i, functions);
		for (int j = 0; j < n_functions; j++)
		{
			fprintf (file, "%i ", functions[j]);
		}
		delete[] functions;

		fprintf (file, "\n");
	}
	fclose (file);
}

template<class Type_Mesh>
void Task<Type_Mesh>::fprint_edges_by_elements (char * file_name)
{
	FILE * file = fopen (file_name, "w");
	int * edges_element;
	int n_edges;

	for (int i = 0, i_end = mesh_pointer->get_n_elements (); i < i_end; i++)
	{
		fprintf (file, "%i:\t", i);
		n_edges = mesh_pointer->get_amount_edges (i);
		edges_element = new int[n_edges];
		mesh_pointer->get_element_edges (i, edges_element);
		for (int j = 0; j < n_edges; j++)
		{
			fprintf (file, "%i ", edges_element[j]);
		}
		delete[] edges_element;
		fprintf (file, "\n");
	}
	fclose (file);
}

template<class Type_Mesh>
void Task<Type_Mesh>::fprint_elements_by_edges (char * file_name)
{
	FILE * file = fopen (file_name, "w");

	int * elements = new int[MAX_EDGES_2D];
	int n_elements;
	for (int i = 0, i_end = mesh_pointer->edges->get_n_entries (); i < i_end; i++)
	{
		fprintf (file, "%i:\t", i);
		n_elements = mesh_pointer->get_elements (i, elements);
		for (int j = 0; j < n_elements; j++)
		{
			fprintf (file, "%i ", elements[j]);
		}
		fprintf (file, "\n");
	}
	delete[] elements;

	fclose (file);
}

template<class Type_Mesh>
void Task<Type_Mesh>::get_mesh_boundaries (double * c0, double * cN)
{
	mesh_pointer->get_0_boundaries (c0);
	mesh_pointer->get_N_boundaries (cN);
}

template<class Type_Mesh>
int Task<Type_Mesh>::get_mesh_dim ()
{
	return mesh_pointer->get_dimentionality ();
}

template<class Type_Mesh>
Mesh_Prototype * Task<Type_Mesh>::get_mesh_pointer ()
{
	return mesh_pointer;
}

template<class Type_Mesh>
double Task<Type_Mesh>::get_Q (int k_system, int time_layer, int k_func)
{
	return previous_time_layers_solutions[k_system][time_layer].getElem (k_func);
}

template<class Type_Mesh>
void Task<Type_Mesh>::get_min_max (int k_system, double * min, double * max)
{
}

template<class Type_Mesh>
void Task<Type_Mesh>::read_solution (int k_system, char * file_name)
{
	// read vector from file into pr_lay_sol[k_system][0]
	FILE * file = fopen (file_name, "r");

	double value;
	for (int i = 0; i < N_functions; i++)
	{
		fscanf (file, "%lf", &value);
		previous_time_layers_solutions[k_system][0].setElem (i, value);
	}
	fclose (file);

	file = fopen ("Result Files Extra//test.txt", "w");
	previous_time_layers_solutions[k_system][0].FPrint (file);
	fclose (file);

}

template<class Type_Mesh>
void Task<Type_Mesh>::get_full_solution (int k_system, MathVector * v)
{
	v->Copy (previous_time_layers_solutions[k_system][0]);
}

template<class Type_Mesh>
int Task<Type_Mesh>::is_symmetrical ()
{
	return symmetrical;
}

template<class Type_Mesh>
void Task<Type_Mesh>::get_lambda_values (char * file_name)
{
	lambda_val.clear ();
	FILE * file = fopen (file_name, "r");

	int area;
	double value;
	int n_areas = mesh_pointer->get_amount_of_areas ();
	for (int i = 0; i < n_areas; i++)
	{
		fscanf (file, "%i %lf", &area, &value);
		lambda_val.insert (std::pair<int, double> (area, value));
	}

	fclose (file);
}

template<class Type_Mesh>
void Task<Type_Mesh>::get_gamma_values (char * file_name)
{
	gamma_val.clear ();
	FILE * file = fopen (file_name, "r");

	int area;
	double value;
	int n_areas = mesh_pointer->get_amount_of_areas ();
	for (int i = 0; i < n_areas; i++)
	{
		fscanf (file, "%i %lf", &area, &value);
		gamma_val.insert (std::pair<int, double> (area, value));
	}

	fclose (file);
}

template<class Type_Mesh>
void Task<Type_Mesh>::calc_derivative (int k_system, int k_var, double * c0, double * cN)
{
	// make new system
	compressed_matrix * matrix = new compressed_matrix (eq_matrixes[k_system]);
	SR = true;

	derivative[k_var].Zero ();
	matrix->Clear ();
	matrix->set_x0 (0.0);

	Matrix * M;
	Matrix * D;
	MathVector * B;
	MathVector * Q;
	bool * functions = new bool[N_functions];
	for (int i = 0, i_end = N_functions; i < i_end; i++)
	{
		functions[i] = false;
	}

	int n_functions;
	int iF, jF;
	// go through elements 
	for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
	{
		n_functions = mesh_pointer->get_amount_non_zero_functions (k_element);
		// if elements is in the area
		if (mesh.element_inside_area (k_element, c0, cN))
		{
			M = new Matrix (n_functions, n_functions);
			D = new Matrix (n_functions, n_functions);
			B = new MathVector (n_functions);
			Q = new MathVector (n_functions);
			mesh_pointer->get_M_local_matrix (k_element, M);
			mesh.elements[k_element]->get_D (k_var, D);

			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (k_element, i);
				if (iF != -1)
				{
					Q->setElem (i, previous_time_layers_solutions[k_system][0].getElem (iF));
				}
			}
			D->MultiplyMatrixByVector (*Q, B);

			// add local_M
			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (k_element, i);
				if (iF != -1)
				{
					for (int j = 0; j < n_functions; j++)
					{
						jF = get_function_global_number (k_element, j);
						// add M into respective places (by functions)
						if (jF != -1)
						{
							matrix->add_to_entry (iF, jF, M->Elem (i, j));
						}
					}
					// add local_B
					matrix->add_to_f_entry (iF, B->getElem (i));
					functions[iF] = true;
				}
			}
			delete M;
			delete B;
			delete Q;
			delete D;
		}
	}
	//matrix->fprint ();
	// faux functions
	for (int i = 0, i_end = N_functions; i < i_end; i++)
	{
		if (!functions[i])
		{
			matrix->clear_row (i);
			matrix->set_entry (i, i, 1.0);
			matrix->set_f_entry (i, 0.0);
		}
	}
	// solve
	matrix->solve_CGM_symm (SOLVER_DECOMP_TYPE_LU, 4);
	// save solution into derivative
	matrix->get_solution (&derivative[k_var]);

	delete[] functions;
	delete matrix;
}

template<class Type_Mesh>
void Task<Type_Mesh>::calc_gradient (int k_system)
{
	int dim = mesh_pointer->get_dimentionality ();
	double * c0 = new double[dim];
	double * cN = new double[dim];
	mesh_pointer->get_0_boundaries (c0);
	mesh_pointer->get_N_boundaries (cN);

	// get derivates by each variable in the whole area
	for (int i = 0; i < dim; i++)
	{
		calc_derivative (k_system, i, c0, cN);
	}
	delete[] c0;
	delete[] cN;

	SR = true;
}

template<class Type_Mesh>
void Task<Type_Mesh>::reset_N_functions (int n)
{
	N_functions = n;
}

template<class Type_Mesh>
void Task<Type_Mesh>::remesh (std::vector <int> elements_to_refine)
{
	// so
	// we have a new mesh up there
	// we have to delete old matrices and solutions
	// and derivatives, but not now
	// everything else stays
	
	// conversion matrix
	compressed_matrix * C = new compressed_matrix ();

	// remesh on a new refined mesh
	Type_Mesh new_mesh;
	new_mesh = mesh;
	int n = new_mesh.refine (elements_to_refine);
	N_functions = n;

	// remesh solutions by solving a system
	int dim = new_mesh.get_dimentionality ();
	int n_functions = new_mesh.get_amount_non_zero_functions (0);
	int n_functions_cur;
	Matrix * M = new Matrix (n_functions, n_functions);
	MathVector * B = new MathVector (n_functions);
	int jF, iF;
	// go by systems
	for (int k_system = 0; k_system < n_systems; k_system++)
	{
		// build new portrait
		// resize happens inside technically
		build_portrait (k_system, &new_mesh);
		// copy portrait of that matrix
		(*C) = eq_matrixes[k_system];

		// go by time layers
		for (int k_time_layer = 0; k_time_layer < time_sampling; k_time_layer++)
		{
			// make matrix C
			C->Clear ();
			for (int k_element = 0, k_element_end = new_mesh.get_n_elements (); k_element < k_element_end; k_element++)
			{
				n_functions_cur = new_mesh.get_amount_non_zero_functions (k_element);
				if (n_functions != n_functions_cur)
				{
					delete B;
					delete M;

					n_functions = n_functions_cur;

					M = new Matrix (n_functions, n_functions);
					B = new MathVector (n_functions);
				}

				// get local M from new mesh
				new_mesh.get_M_local_matrix (k_element, M);
				// build local B
				{
					double jac = 0.0; // jacobian
					int n_integr_points = new_mesh.amount_of_integration_points (k_element);
					double * weigths = new double[n_integr_points];
					double ** points = new double *[n_integr_points];
					for (int i = 0; i < n_integr_points; i++)
					{
						points[i] = new double[dim];
					}

					new_mesh.integration_points (k_element, points, weigths, &jac);

					double val;
					double value;
					for (int i = 0; i < n_functions; i++)
					{
						// get vector of integrals f * function_i
						val = 0.0;
						// get theta values * basis_functions in those points
						for (int k = 0; k < n_integr_points; k++)
						{
							// old solution can be picked through old mesh which hasn't been replaces yet and previous_sol 
							get_solution_in_point (k_system, k_time_layer, points[k], &value);
							// sum them multiplying by weigths
							val += new_mesh.get_basis_function_value (k_element, i, points[k]) * value * weigths[k];
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

				// go by element's functions
				for (int i = 0; i < n_functions; i++)
				{
					iF = new_mesh.get_function_global_number (k_element, i);
					if (iF != -1)
					{
						for (int j = 0; j < n_functions; j++)
						{
							jF = new_mesh.get_function_global_number (k_element, j);

							if (jF != -1)
							{
								// add M into respective places (by def_nodes)
								C->add_to_entry (iF, jF, M->Elem (i, j));
							}
						}

						// put B into respective places (by def_nodes)
						C->add_to_f_entry (iF, B->getElem (i));
					}
				}
			}

			C->solve_GMRES (1, 50);
			// save solution
			//previous_time_layers_solutions[k_system][k_time_layer].Print ();
			previous_time_layers_solutions[k_system][k_time_layer].setSize (N_functions);
			C->get_solution (&previous_time_layers_solutions[k_system][k_time_layer]);
			//previous_time_layers_solutions[k_system][k_time_layer].Print ();
		}
	}
	delete M;
	delete B;
	delete C;

	// finish with changing the mesh
	mesh.copy (new_mesh);
	mesh.output ();
}

template<class Type_Mesh>
bool Task<Type_Mesh>::get_isoline_section (int k_system, int var, int k_element, double value, double * c1, double * c2)
{
	// get numbers of global functions from element
	int n_functions = mesh.get_amount_non_zero_functions (k_element);
	// get requested values of those functions (don't forget that number can be -1)
	double * q = new double[n_functions];
	for (int i = 0; i < n_functions; i++)
	{
		int i_function = mesh.get_function_global_number (k_element, i);
		q[i] = 0.0;
		if (i_function != -1)
		{
			// var = 0..n - derivative, < 0 - solution value; gradient norm can not be calculated that way
			if (var < 0)
			{
				q[i] = previous_time_layers_solutions[k_system][0].getElem (i_function);
			}
			else
			{
				if (var < get_mesh_dim ())
				{
					q[i] = derivative[var].getElem (i_function);
				}
			}
		}
	}
	// call mesh to calculate respective points
	bool result = mesh.get_isoline_section (k_element, q, value, c1, c2);
	delete[] q;
	return result;
}

template<class Type_Mesh>
void Task<Type_Mesh>::set_solution (int k_system, int k_time_layer, char * file_name)
{
	FILE * file = fopen (file_name, "r");
	previous_time_layers_solutions[k_system][k_time_layer].getFromFile (file);
	fclose (file);
}

template<class Type_Mesh>
void Task<Type_Mesh>::build_system (int k_system)
{
	Matrix * G;
	Matrix * M;
	MathVector * B;
	int n_functions;
	int n_functions_cur;
	int area;
	double lambda, gamma;
	int iF, jF;
	int dim = mesh_pointer->get_dimentionality ();
	// code for test purposes
	FILE * fileG = fopen ("Result Files Extra//G.txt", "w");
	//FILE * fileM = fopen ("Result Files Extra//M.txt", "w");

	n_functions = mesh_pointer->get_amount_non_zero_functions (0);
	G = new Matrix (n_functions, n_functions);
	M = new Matrix (n_functions, n_functions);
	B = new MathVector (n_functions);

	for (int k_element = 0, k_element_end = mesh_pointer->get_n_elements (); k_element < k_element_end; k_element++) // go through elements
	{
		n_functions_cur = mesh_pointer->get_amount_non_zero_functions (k_element);
		if (n_functions != n_functions_cur)
		{
			delete B;
			delete G;
			delete M;

			n_functions = n_functions_cur;

			G = new Matrix (n_functions, n_functions);
			M = new Matrix (n_functions, n_functions);
			B = new MathVector (n_functions);
		}

		// get G, M matrices 
		mesh_pointer->get_G_local_matrix (k_element, G);
		mesh_pointer->get_M_local_matrix (k_element, M);

		//// code for test purposes, save G, M 
		//if (k_element == 0)
		//{
		//	//G->FPrint (fileG);
		////	M->FPrint (fileM);
		//}

		// get B-vector
		get_local_B (k_system, k_element, B);

		// get average lambda/gamma on the element
		area = mesh.elements[k_element]->get_area ();
		lambda = 0.0;
		gamma = 0.0;
		for (int k_function = 0; k_function < n_functions; k_function++)
		{
			lambda += Lambda (k_system, k_function, area);
			gamma += Gamma (k_system, k_function, area);
		}
		lambda /= (double)n_functions;
		gamma /= (double)n_functions;

		// go by element's functions
		for (int i = 0; i < n_functions; i++)
		{
			iF = get_function_global_number (k_element, i);
			if (iF != -1)
			{
				for (int j = 0; j < n_functions; j++)
				{
					jF = get_function_global_number (k_element, j);

					if (iF != -1)
					{
						// add G into respective places (by def_nodes)
						eq_matrixes[k_system].add_to_entry (iF, jF, lambda * G->Elem (i, j));
						// add M into respective places (by def_nodes)
						eq_matrixes[k_system].add_to_entry (iF, jF, gamma * M->Elem (i, j));
					}
				}

				// put B into respective places (by def_nodes)
				eq_matrixes[k_system].add_to_f_entry (iF, B->getElem (i));
			}
		}
	}

	//eq_matrixes[k_system].fprint ();
	delete B;
	delete G;
	delete M;
	// code for test purposes
	fclose (fileG);
	//fclose (fileM);
}

template<class Type_Mesh>
void Task<Type_Mesh>::apply_boundary_conditions (int k_system)
{
	apply_second_boundary_conditions (k_system);
	apply_first_boundary_conditions (k_system);
}

template<class Type_Mesh>
void Task<Type_Mesh>::apply_first_boundary_conditions (int k_system)
{
	int boundary;
	int iF;
	MathVector * FCondition;
	int * functions;
	int f_amount;
	int k_element;
	int n1, n2;

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

				// get first condition from element
				get_FCondition_edge (k_system, k_element, n1, n2, functions, conditions[k_system][boundary].Function (), FCondition);

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

template<class Type_Mesh>
void Task<Type_Mesh>::apply_second_boundary_conditions (int k_system)
{
	int boundary;
	int iF;
	MathVector * Theta;
	int * functions;
	int f_amount;
	int k_element;
	int n1, n2;

	// second condition by edges 
	for (int i_edge = 0, i_end = mesh_pointer->edges->get_n_entries (); i_edge < i_end; i_edge++)
	{
		// get edges's nodes
		mesh_pointer->edges->get_edge_nodes (i_edge, &n1, &n2);
		// get boundary number for those nodes
		boundary = mesh_pointer->function_boundary_edge (n1, n2);
		if (boundary != -1)
		{
			if (conditions[k_system][boundary].Type () == 2)
			{
				// find element with that edge
				k_element = mesh_pointer->belonging_element (n1, n2);
				// get amount of values to put into f (from element)
				f_amount = mesh_pointer->get_amount_second_condition (k_element);
				Theta = new MathVector (f_amount);
				Theta->Zero ();

				functions = new int[f_amount];
				// get theta from element
				mesh_pointer->get_edge_functions (k_element, n1, n2, functions);
				get_Theta_edge (k_system, k_element, n1, n2, functions, boundary, Theta);
				// add it to f in places by functions
				for (int i = 0; i < f_amount; i++)
				{
					iF = get_function_global_number (k_element, functions[i]);
					if (iF != -1)
					{
						eq_matrixes[k_system].add_to_f_entry (iF, Theta->getElem (i));
					}
				}

				delete Theta;
				delete[] functions;
			}
		}
	}
}

template<class Type_Mesh>
void Task<Type_Mesh>::get_Theta_edge (int k_system, int k_element, int n1, int n2, int * functions, int boundary, MathVector * Theta)
{
	// TEST
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

	double val;
	for (int i = 0; i < n_funcs; i++)
	{
		val = 0.0;
		// get theta values * basis_functions in those points
		for (int j = 0; j < n_integr_points; j++)
		{
			// sum them multiplying by weigths
			val += mesh.elements[k_element]->get_basis_function_value (functions[i], points[j]) * function_Theta (k_system, points[j], area, boundary) * weigths[j];
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

template<class Type_Mesh>
double Task<Type_Mesh>::get_derivative_first (int k_system, int k_var, double * coordinates)
{
	bool found = false;
	double value = 0.0;
	// find element that contains the point 
	int p_element = mesh_pointer->point_inside (coordinates);
	int iF;
	if (p_element != -1)
	{
		int dim = mesh_pointer->get_dimentionality ();

		// get local functions' values of the element
		int n_functions = mesh_pointer->get_amount_non_zero_functions (p_element);
		MathVector * local_func = new MathVector (n_functions);
		mesh_pointer->get_local_function_first_derivative_values (p_element, k_var, coordinates, local_func);

		if (k_var < dim)
		{
			// make vector of solutions for those functions
			MathVector * solution = new MathVector (n_functions);
			for (int i = 0; i < n_functions; i++)
			{
				iF = get_function_global_number (p_element, i);
				solution->setElem (i, previous_time_layers_solutions[k_system][0].getElem (iF));
								
			}

			// multiply them by the solution's values in those nodes
			value = local_func->Scalar_Product (*solution);

			delete local_func;
			delete solution;
		}
		else
		{
			// make vector of solutions in those nodes
			MathVector * solution = new MathVector (n_functions);
			MathVector * gradient = new MathVector (n_functions);
			int global_function;
			double r_var;
			for (int k = 0; k < dim; k++)
			{
				mesh_pointer->get_local_function_first_derivative_values (p_element, k, coordinates, local_func);
				for (int i = 0; i < n_functions; i++)
				{
					global_function = get_function_global_number (p_element, i);
					solution->setElem (i, previous_time_layers_solutions[k_system][0].getElem (global_function));
				}
				// multiply them by the solution's values in those nodes
				r_var = local_func->Scalar_Product (*solution);
				value += pow (r_var, 2.0);
			}
			value = sqrt (value);
			delete solution;
			delete gradient;
		}
	}
	return value;
}

template<class Type_Mesh>
double Task<Type_Mesh>::get_derivative_second (int k_system, int k_var1, int k_var2, double h1, double h2, double * coordinates)
{
	double r = 0.0;
	double val;
	int dim = mesh_pointer->get_dimentionality ();
	double * c = new double[dim];
	double * point = new double[dim];
	// copy coordinates of the requested point 
	for (int i = 0; i < dim; i++)
	{
		point[i] = coordinates[i];
	}

	if (k_var1 < dim)
	{
		mesh_pointer->get_N_boundaries (c);
		if (fabs (c[k_var1] - point[k_var1]) < 1.5 * h1)// if requested point is too close to the right border, then first point should be on the left
		{
			// get derivative in point by another variable
			val = get_derivative_first (k_system, k_var2, point);
			// add it, since we are going inside the area
			r += val;
			// go inside the area
			point[k_var1] -= h1;
			// get derivative there
			val = get_derivative_first (k_system, k_var2, point);
			// substract it
			r -= val;
			// divide by h to get derivative
			r /= h1;
		}
		else
		{
			mesh_pointer->get_0_boundaries (c);
			if (fabs (c[k_var1] - point[k_var1]) < 1.5 * h1)// if requested point is too close to the left border
			{
				// get derivative in point by another variable
				val = get_derivative_first (k_system, k_var2, point);
				// substract it, since we are going inside the area
				r -= val;
				// go inside the area
				point[k_var1] += h1;
				// get derivative there
				val = get_derivative_first (k_system, k_var2, point);
				// add it
				r += val;
				// divide by h to get derivative
				r /= h1;
			}
			else
			{
				point[k_var1] -= h1;
				// get derivative in point by another variable
				val = get_derivative_first (k_system, k_var2, point);
				// substract it
				r -= val;
				// go inside the area
				point[k_var1] += 2.0 * h1;
				// get solution there
				val = get_derivative_first (k_system, k_var2, point);
				// add it
				r += val;
				// divide by h to get derivative
				r /= (2.0 * h1);
			}
		}
	}
	delete[] c;
	delete[] point;
	return r;
}

template<class Type_Mesh>
double Task<Type_Mesh>::get_derivative_SR (int k_system, int k_var, double * coordinates)
{
	double r = 0.0;
	bool found = false;
	// find element that contains the point 
	int p_element = mesh_pointer->point_inside (coordinates);
	if (p_element == -1)
	{
		return false;
	}
	// get local functions' values of the element
	int func_amount = mesh_pointer->get_amount_non_zero_functions (p_element);
	MathVector * local_func = new MathVector (func_amount);
	mesh_pointer->get_local_function_values (p_element, coordinates, local_func);

	int dim = mesh_pointer->get_dimentionality ();

	if (k_var < dim)
	{
		// make vector of solutions in those nodes
		MathVector * solution = new MathVector (func_amount);
		int global_function;
		for (int i = 0; i < func_amount; i++)
		{
			global_function = get_function_global_number (p_element, i);
			solution->setElem (i, derivative[k_var].getElem (global_function));
		}

		// multiply them by the solution's values in those nodes
		r = local_func->Scalar_Product (*solution);

		delete solution;
	}
	else // else make gradient norm
	{
		// make vector of solutions in those nodes
		MathVector * solution = new MathVector (func_amount);
		MathVector * gradient = new MathVector (func_amount);
		int global_function;
		double r_var;
		for (int k = 0; k < dim; k++)
		{
			for (int i = 0; i < func_amount; i++)
			{
				global_function = get_function_global_number (p_element, i);
				solution->setElem (i, derivative[k].getElem (global_function));
			}
			// multiply them by the solution's values in those nodes
			r_var = local_func->Scalar_Product (*solution);
			r += pow (r_var, 2.0);
		}
		r = sqrt (r);
		delete solution;
		delete gradient;
	}
	delete local_func;
	return r;
}

template<class Type_Mesh>
double Task<Type_Mesh>::Lambda (int k_system, int k_function, int area)
{
	return lambda_val[area];
}

template<class Type_Mesh>
double Task<Type_Mesh>::Gamma (int k_system, int k_function, int area)
{
	return gamma_val[area];
}

template<class Type_Mesh>
void Task<Type_Mesh>::get_local_B (int k_system, int k_element, MathVector * B)
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
	for (int i = 0; i < n_funcs; i++)
	{
		// get vector of integrals f * function_i
		val = 0.0;
		// get theta values * basis_functions in those points
		for (int k = 0; k < n_integr_points; k++)
		{
			// sum them multiplying by weigths
			val += mesh.elements[k_element]->get_basis_function_value (i, points[k]) * function_f (k_system, points[k], area) * weigths[k];
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

template<class Type_Mesh>
double Task<Type_Mesh>::function_f (int k_system, double * coordinates, int area)
{
	return coordinates[0] + coordinates[1];
}

template<class Type_Mesh>
double Task<Type_Mesh>::function_FCondition (int k_system, double * coordinates, int area, int boundary)
{
	return 1.0;
}

template<class Type_Mesh>
void Task<Type_Mesh>::get_FCondition_edge (int k_system, int k_element, int n1, int n2, int * functions, int boundary, MathVector * FCondition)
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

	// make matrix 
	Matrix * A = new Matrix (n_funcs, n_funcs);
	MathVector * B = new MathVector (n_funcs);

	double val;	
	for (int i = 0; i < n_funcs; i++)
	{
		for (int j = 0; j < n_funcs; j++)
		{
			val = 0.0;
			// get theta values * basis_functions in those points
			for (int k = 0; k < n_integr_points; k++)
			{
				// sum them multiplying by weigths
				val += mesh.elements[k_element]->get_basis_function_value (functions[i], points[k]) * mesh.elements[k_element]->get_basis_function_value (functions[j], points[k]) * weigths[k];
			}
			// multiply by jac
			val *= jac;
			A->setElem (i, j, val);
		}

		// get vector of integrals f * function_i
		val = 0.0;
		// get theta values * basis_functions in those points
		for (int k = 0; k < n_integr_points; k++)
		{
			// sum them multiplying by weigths
			val += mesh.elements[k_element]->get_basis_function_value (functions[i], points[k]) * function_FCondition (k_system, points[k], area, boundary) * weigths[k];
		}
		// multiply by jac
		val *= jac;
		B->setElem (i, val);
	}
	// solve it
	A->solve (FCondition, *B);
	delete A;
	delete B;
	delete[] weigths;
	for (int i = 0; i < n_integr_points; i++)
	{
		delete[] points[i];
	}
	delete[] points;
}

template<class Type_Mesh>
double Task<Type_Mesh>::function_Theta (int k_system, double * coordinates, int area, int boundary)
{
	return 0.0;
}

template<class Type_Mesh>
int Task<Type_Mesh>::get_function_global_number (int k_element, int k_func)
{
	return mesh.elements[k_element]->get_function_global_number (k_func);
}

template<class Type_Mesh>
void Task<Type_Mesh>::get_element_functions (int k_element, int * functions)
{
	int n_functions = mesh_pointer->get_amount_non_zero_functions (k_element);
	for (int i = 0; i < n_functions; i++)
	{
		functions[i] = get_function_global_number (k_element, i);
	}
}

template<class Type_Mesh>
double Task<Type_Mesh>::get_function_value (int k_system, int k_function, double * coordinates, int area)
{
	return 0.0;
}

template<class Type_Mesh>
void Task<Type_Mesh>::get_function_coefficients (int k_system, int k_element, int function, MathVector * coef)
{
	Matrix * M;
	MathVector * B;
	int dim = mesh_pointer->get_dimentionality ();
	int n_functions = mesh_pointer->get_amount_non_zero_functions (k_element);
	M = new Matrix (n_functions, n_functions);
	B = new MathVector (n_functions);

	mesh_pointer->get_M_local_matrix (k_element, M);

	// get B
	double jac = 0.0; // jacobian
	int n_integr_points = mesh.elements[k_element]->amount_of_integration_points ();
	double * weigths = new double[n_integr_points];
	double ** points = new double *[n_integr_points];
	for (int i = 0; i < n_integr_points; i++)
	{
		points[i] = new double[dim];
	}
	// get integration points for the element
	mesh.elements[k_element]->integration_points (points, weigths, &jac);
	double val;
	int area = mesh_pointer->get_area (k_element);
	// go by element's functions

	double func;
	for (int i = 0; i < n_functions; i++)
	{
		val = 0.0;
		// integrate them * basis_functions on the element
		for (int j = 0; j < n_integr_points; j++)
		{
			// sum them multiplying by weigths
			func = get_function_value (k_system, function, points[j], area);
			val += mesh.elements[k_element]->get_basis_function_value (i, points[j]) * func * weigths[j];
		}
		// multiply by jac
		val *= jac;
		B->setElem (i, val);
	}

	// get coef
	M->solve (coef, *B);

	delete M;
	delete B;

	delete[] weigths;
	for (int i = 0; i < n_integr_points; i++)
	{
		delete[] points[i];
	}
	delete[] points;
}

template<class Type_Mesh>
Task<Type_Mesh>::Task ()
{
	non_linear = false;
	symmetrical = 0;
	N_functions = 0;
	eq_matrixes = NULL;
	time_layers = NULL;
	previous_time_layers_solutions = NULL;
	conditions = NULL;
	non_linear_layers_solutions = NULL;
	lambda_val.insert (std::pair<int, double> (1, 1.0));
	gamma_val.insert (std::pair<int, double> (1, 1.0));

	prepared = false;
	SR = false;

	painter = NULL;
}

template<class Type_Mesh>
Task<Type_Mesh>::Task (const Task & task)
{
	// mesh section
	mesh = task.mesh;

	// time section
	if (n_time_layers != task.n_time_layers)
	{
		if (time_layers != NULL)
			delete[] time_layers;

		n_time_layers = task.n_time_layers;

		if (n_time_layers > 0)
			time_layers = new double[n_time_layers];

		for (int i = 0; i < n_time_layers; i++)
		{
			time_layers[i] = task.time_layers[i];
		}
	}

	// systems' section
	if (eq_matrixes != NULL)
		delete[] eq_matrixes;

	N_functions = task.N_functions;
	n_systems = task.n_systems;
	if (n_systems > 0)
	{
		eq_matrixes = new compressed_matrix[n_systems];
		for (int i = 0; i < n_systems; i++)
		{
			eq_matrixes[i] = task.eq_matrixes[i];
		}
	}

	// conditions section
	if (conditions != NULL)
	{
		for (int i = 0; i < n_systems; i++)
			delete[] conditions[i];
	}
	delete[] conditions;
	if (task.conditions != NULL)
	{
		conditions = new Condition *[n_systems];
		for (int i = 0; i < n_systems; i++)
		{
			conditions[i] = new Condition[mesh_pointer->get_dimentionality () * 2];
			for (int k = 0; k < mesh_pointer->get_dimentionality () * 2; k++)
			{
				conditions[i][k].Copy (task.conditions[i][k]);
			}
		}
	}
}

template<class Type_Mesh>
Task<Type_Mesh>::~Task ()
{
	if (time_layers != NULL)
		delete[] time_layers;
	if (eq_matrixes != NULL)
		delete[] eq_matrixes;
	if (non_linear_layers_solutions != NULL)
		delete[] non_linear_layers_solutions;
	if (previous_time_layers_solutions != NULL)
	{
		for (int k = 0; k < n_systems; k++)
		{
			delete[] previous_time_layers_solutions[k];
		}
		delete[] previous_time_layers_solutions;
	}

	if (conditions != NULL)
	{
		for (int k = 0; k < n_systems; k++)
		{
			delete[] conditions[k];
		}
		delete[] conditions;
	}
	if (derivative != NULL)
		delete[] derivative;
}

template<class Type_Mesh>
void Task<Type_Mesh>::painter_pointer (Painter * p)
{
	painter = p;
}

template<class Type_Mesh>
bool Task<Type_Mesh>::solve_task (int method, int d_type, int depth)
{
	// TEST
	if (!prepared)
	{
		printf ("ERROR: missing preparation\n");
		return false;
	}

	int nonlinear_iter = 1;
	if (non_linear)
		nonlinear_iter = MAX_NONLINEAR_ITER;

	non_linear_layers_solutions = new MathVector[n_systems];
	for (int k_system = 0; k_system < n_systems; k_system++)
	{
		non_linear_layers_solutions[k_system].setSize (eq_matrixes[k_system].Size ()); // default starting solution point is 0
	}

	double * discr = new double[n_systems];
	double * prev_discr = new double[n_systems];
	bool stable = false; // indicator of reaching stable condition
	bool found = false;
	int solver_iterations;
	double sum_discr;
	double T_change;
	int unchanged_discr;
	FILE * file_T = fopen ("Result//Melt_SF//T_change.txt", "w");
	print_solutions ();
	
	for (current_time_layer = 1; current_time_layer < n_time_layers && !stable; current_time_layer++)
	{
		printf ("current time: %.7lf\n", time_layers[current_time_layer]);
		//printf ("%.7lf ", time_layers[current_time_layer]);
		// save previous time layers solution for non-linearity
		for (int k_system = 0; k_system < n_systems; k_system++)
		{
			non_linear_layers_solutions[k_system].Copy ((previous_time_layers_solutions[k_system][0]));
		}
		found = false;
		for (int i = 0; i < nonlinear_iter && !found; i++)
		{
			printf ("\tnonlinear iteration: %i\n", i);
			// iterate through each system of equations
			for (int k_system = 0; k_system < n_systems; k_system++)
			{
				eq_matrixes[k_system].Clear ();
				//eq_matrixes[k_system].set_starting_point (&(non_linear_layers_solutions[k_system])); // set last layer solution as x0 for solver
				eq_matrixes[k_system].set_starting_point (&(previous_time_layers_solutions[k_system][0])); // set last layer solution as x0 for solver

				build_system (k_system);  // reset system for current layer
				apply_boundary_conditions (k_system); // apply boundary conditions
				int solver_param[] = {method, d_type, depth, 0, SOLVER_MKL_NO };
				solver_iterations = eq_matrixes[k_system].solve (solver_param); // solve it	
				discr[k_system] = relax (k_system, &(non_linear_layers_solutions[k_system]));
			}
			if (i > 0)
			{
				sum_discr = 0.0;
				unchanged_discr = 0;
				for (int k_system = 0; k_system < n_systems; k_system++)
				{
					sum_discr += discr[k_system];
					if (fabs (discr[k_system] - prev_discr[k_system]) < DISCR_CHANGE)
						unchanged_discr++;
				}
				if (sum_discr / (double) n_systems < NONLINEAR_DISCR)
					found = true;
				if (unchanged_discr >= n_systems)
					found = true;
				T_change = non_linear_layers_solutions[0].sum_sq_dif (previous_time_layers_solutions[0][0]);
				printf ("%i: %e %e %e\tT: %e %e\n", unchanged_discr, discr[0] - prev_discr[0], discr[1] - prev_discr[1], discr[2] - prev_discr[2], T_change, T_change / non_linear_layers_solutions[0].Norm ());
				fprintf (file_T, "%lf %e %e\n", time_layers[current_time_layer], T_change, T_change / non_linear_layers_solutions[0].Norm ());
			}
			for (int k_system = 0; k_system < n_systems; k_system++)
			{
				prev_discr[k_system] = discr[k_system];
			}
		}
		// move layers
		for (int k_system = 0; k_system < n_systems; k_system++)
			next_time_layer (k_system, &(non_linear_layers_solutions[k_system]));

		print_solutions ();
		print_extra_data ();
	}
	current_time_layer--;
	fclose (file_T);

	delete[] discr;
	delete[] prev_discr;
	delete[] non_linear_layers_solutions;
	non_linear_layers_solutions = NULL;
	return true;
}

template<class Type_Mesh>
bool Task<Type_Mesh>::solve_task (int * solver_param)
{
	// TEST
	if (!prepared)
	{
		printf ("ERROR: missing preparation\n");
		return false;
	}

	int nonlinear_iter = 1;
	if (non_linear)
		nonlinear_iter = MAX_NONLINEAR_ITER;

	non_linear_layers_solutions = new MathVector[n_systems];
	for (int k_system = 0; k_system < n_systems; k_system++)
	{
		non_linear_layers_solutions[k_system].setSize (eq_matrixes[k_system].Size ()); // default starting solution point is 0
	}

	double * discr = new double[n_systems];
	double * prev_discr = new double[n_systems];
	bool stable = false; // indicator of reaching stable condition
	bool found = false;
	int solver_iterations;
	double sum_discr;
	double T_change;
	int unchanged_discr;
	FILE * file_T = fopen ("Result//Melt_SF//T_change.txt", "w");
	for (current_time_layer = 1; current_time_layer < n_time_layers && !stable; current_time_layer++)
	{
		printf ("current time: %.7lf\n", time_layers[current_time_layer]);
		//printf ("%.7lf ", time_layers[current_time_layer]);
		// save previous time layers solution for non-linearity
		for (int k_system = 0; k_system < n_systems; k_system++)
		{
			non_linear_layers_solutions[k_system].Copy ((previous_time_layers_solutions[k_system][0]));
		}
		found = false;
		for (int i = 0; i < nonlinear_iter && !found; i++)
		{
			printf ("\tnonlinear iteration: %i\n", i);
			// iterate through each system of equations
			for (int k_system = 0; k_system < n_systems; k_system++)
			{
				eq_matrixes[k_system].Clear ();
				//eq_matrixes[k_system].set_starting_point (&(non_linear_layers_solutions[k_system])); // set last layer solution as x0 for solver
				eq_matrixes[k_system].set_starting_point (&(previous_time_layers_solutions[k_system][0])); // set last layer solution as x0 for solver

				build_system (k_system);  // reset system for current layer
				apply_boundary_conditions (k_system); // apply boundary conditions
				solver_iterations = eq_matrixes[k_system].solve (solver_param); // solve it	
				discr[k_system] = relax (k_system, &(non_linear_layers_solutions[k_system]));
			}
			if (i > 0)
			{
				sum_discr = 0.0;
				unchanged_discr = 0;
				for (int k_system = 0; k_system < n_systems; k_system++)
				{
					sum_discr += discr[k_system];
					if (fabs (discr[k_system] - prev_discr[k_system]) < DISCR_CHANGE)
						unchanged_discr++;
				}
				if (sum_discr / (double)n_systems < NONLINEAR_DISCR)
					found = true;
				if (unchanged_discr >= n_systems)
					found = true;
				T_change = non_linear_layers_solutions[0].sum_sq_dif (previous_time_layers_solutions[0][0]);
				printf ("%i: %e %e %e\tT: %e %e\n", unchanged_discr, discr[0] - prev_discr[0], discr[1] - prev_discr[1], discr[2] - prev_discr[2], T_change, T_change / non_linear_layers_solutions[0].Norm ());
				fprintf (file_T, "%lf %e %e\n", time_layers[current_time_layer], T_change, T_change / non_linear_layers_solutions[0].Norm ());
			}
			for (int k_system = 0; k_system < n_systems; k_system++)
			{
				prev_discr[k_system] = discr[k_system];
			}
		}
		// move layers
		for (int k_system = 0; k_system < n_systems; k_system++)
			next_time_layer (k_system, &(non_linear_layers_solutions[k_system]));

		print_solutions ();
		print_extra_data ();
	}
	current_time_layer--;
	fclose (file_T);

	delete[] discr;
	delete[] prev_discr;
	delete[] non_linear_layers_solutions;
	non_linear_layers_solutions = NULL;
	return true;
}

template<class Type_Mesh>
bool Task<Type_Mesh>::solve_task (int solver_param[][5])
{
	// TEST
	if (!prepared)
	{
		printf ("ERROR: missing preparation\n");
		return false;
	}

	int nonlinear_iter = 1;
	//if (non_linear)
	//	nonlinear_iter = MAX_NONLINEAR_ITER;
	non_linear_layers_solutions = new MathVector[n_systems];
	for (int k_system = 0; k_system < n_systems; k_system++)
	{
		non_linear_layers_solutions[k_system].setSize (eq_matrixes[k_system].Size ()); // default starting solution point is 0
	}

	double * discr = new double[n_systems];
	double * prev_discr = new double[n_systems];
	bool stable = false; // indicator of reaching stable condition
	bool found = false;
	int solver_iterations;
	double sum_discr;
	double T_change;
	int unchanged_discr;
	FILE * file_T = fopen ("Result//Melt_SF//T_change.txt", "w");
	for (current_time_layer = 1; current_time_layer < n_time_layers && !stable; current_time_layer++)
	{
		printf ("current time: %.7lf\n", time_layers[current_time_layer]);
		//printf ("%.7lf ", time_layers[current_time_layer]);
		// save previous time layers solution for non-linearity
		for (int k_system = 0; k_system < n_systems; k_system++)
		{
			non_linear_layers_solutions[k_system].Copy ((previous_time_layers_solutions[k_system][0]));
		}
		found = false;
		for (int i = 0; i < nonlinear_iter && !found; i++)
		{
			printf ("\tnonlinear iteration: %i\n", i);
			// iterate through each system of equations
			for (int k_system = 0; k_system < n_systems; k_system++)
			{
				eq_matrixes[k_system].Clear ();
				//eq_matrixes[k_system].set_starting_point (&(non_linear_layers_solutions[k_system])); // set last layer solution as x0 for solver
				eq_matrixes[k_system].set_starting_point (&(previous_time_layers_solutions[k_system][0])); // set last layer solution as x0 for solver

				build_system (k_system);  // reset system for current layer
				apply_boundary_conditions (k_system); // apply boundary conditions
				solver_iterations = eq_matrixes[k_system].solve (solver_param[k_system]); // solve it	
				discr[k_system] = relax (k_system, &(non_linear_layers_solutions[k_system]));
			}
			if (i > 0)
			{
				sum_discr = 0.0;
				unchanged_discr = 0;
				for (int k_system = 0; k_system < n_systems; k_system++)
				{
					sum_discr += discr[k_system];
					if (fabs (discr[k_system] - prev_discr[k_system]) < DISCR_CHANGE)
						unchanged_discr++;
				}
				if (sum_discr / (double)n_systems < NONLINEAR_DISCR)
					found = true;
				if (unchanged_discr >= n_systems)
					found = true;
				T_change = non_linear_layers_solutions[0].sum_sq_dif (previous_time_layers_solutions[0][0]);
				printf ("%i: %e %e %e\tT: %e %e\n", unchanged_discr, discr[0] - prev_discr[0], discr[1] - prev_discr[1], discr[2] - prev_discr[2], T_change, T_change / non_linear_layers_solutions[0].Norm ());
				fprintf (file_T, "%lf %e %e\n", time_layers[current_time_layer], T_change, T_change / non_linear_layers_solutions[0].Norm ());
			}
			for (int k_system = 0; k_system < n_systems; k_system++)
			{
				prev_discr[k_system] = discr[k_system];
			}
		}
		// move layers
		for (int k_system = 0; k_system < n_systems; k_system++)
			next_time_layer (k_system, &(non_linear_layers_solutions[k_system]));

		print_solutions ();
		print_extra_data ();
	}
	current_time_layer--;
	fclose (file_T);

	delete[] discr;
	delete[] prev_discr;
	delete[] non_linear_layers_solutions;
	non_linear_layers_solutions = NULL;
	return true;
}

template<class Type_Mesh>
void Task<Type_Mesh>::read_solution (char * file_name, int k_system, int k_time_layer)
{
	double value;
	previous_time_layers_solutions[k_system][k_time_layer].setSize (N_functions);

	FILE * file = fopen (file_name, "r");
	for (int i = 0; i < N_functions; i++)
	{
		fscanf (file, "%lf", &value);
		previous_time_layers_solutions[k_system][k_time_layer].setElem (i, value);
	}
	fclose (file);
}

template<class Type_Mesh>
void Task<Type_Mesh>::set_symmetry (int Symmetrical)
{
	symmetrical = Symmetrical;
}
