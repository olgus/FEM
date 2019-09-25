#include "Task_3D_simple.h"

template class Task_3D <Mesh<Node<Point_3D>, Element>>;
template class Task_3D <Mesh<Node_3D, Element>>;

template class Task_3D <Prismatic_Mesh>;
template class Task_3D <Cubic_Mesh>;

template <class Type_Mesh>
void Task_3D<Type_Mesh>::apply_first_boundary_conditions (int k_system)
{
	bool * replaced = new bool[N_functions];
	for (int i = 0; i < N_functions; i++)
	{
		replaced[i] = false;
	}

	double coordinates[3];
	int * def_nodes;
	int iF;
	double value;
	int boundary;
	int area;
	int n_functions;
	// go by elements
	int n_elements = mesh_pointer->get_n_elements ();
	for (int i_element = 0; i_element < n_elements; i_element++)
	{
		def_nodes = new int[mesh_pointer->get_amount_of_def_nodes (i_element)];
		mesh_pointer->get_def_nodes (i_element, def_nodes);

		// go by functions of the element
		n_functions = mesh_pointer->get_amount_non_zero_functions (i_element);
		for (int i_func = 0; i_func < n_functions; i_func++)	
		{
			boundary = mesh_pointer->function_boundary (i_element, i_func);
			iF = get_function_global_number (i_element, i_func);
			// if the function is on the boundary of the mesh
			if (!replaced[iF] && boundary != -1)
			{
				if (conditions[k_system][boundary].Type () == 1)
				{
					area = mesh_pointer->get_area (i_element);
					mesh_pointer->get_node_coordinates (def_nodes[i_func], coordinates);
					// get value for the def_node
					value = function_FCondition (k_system, coordinates, area, boundary);
					// replace it
					replaced[iF] = true;
					eq_matrixes[k_system].clear_row (iF); // replace row
					eq_matrixes[k_system].set_entry (iF, iF, 1.0);
					eq_matrixes[k_system].set_f_entry (iF, value);
				}
			}
		}
		delete[] def_nodes;
	}

	delete[] replaced;
}

template<class Type_Mesh>
void Task_3D<Type_Mesh>::apply_second_boundary_conditions (int k_system)
{
	// do nothing since edges for 3D elements are not realized
	// funny, i know
}

template<class Type_Mesh>
void Task_3D<Type_Mesh>::get_local_B (int k_system, int k_element, MathVector * B)
{
	// B = M * vector(f(def_nodes))

	int n_functions = mesh_pointer->get_amount_non_zero_functions (k_element);
	Matrix * M = new Matrix (n_functions, n_functions);
	mesh_pointer->get_M_local_matrix (k_element, M);

	MathVector * F = new MathVector (n_functions);
	int * nodes = new int[n_functions];
	mesh_pointer->get_def_nodes (k_element, nodes);
	double coordinates[3];
	int area = mesh_pointer->get_area (k_element);

	for (int i = 0; i < n_functions; i++)
	{
		mesh_pointer->get_node_coordinates (nodes[i], coordinates);
		F->setElem (i, function_f (k_system, coordinates, area));
	}
	M->MultiplyMatrixByVector (*F, B);

	delete F;
	delete M;
	delete[] nodes;
}

template<class Type_Mesh>
Task_3D<Type_Mesh>::Task_3D ()
{
}

template<class Type_Mesh>
Task_3D<Type_Mesh>::Task_3D (const Task_3D & task)
{
}

template<class Type_Mesh>
Task_3D<Type_Mesh>::~Task_3D ()
{
}
