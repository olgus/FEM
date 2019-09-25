#pragma once
#include <vector>

#include "MathVector.h"
#include "Matrix.h"
#include "Support.h"

class Mesh_Prototype
{
public:
	Mesh_Prototype () {}; // constructor
	Mesh_Prototype (const Mesh_Prototype & mp) {}; // constructor-copy
	virtual ~Mesh_Prototype () {}; // destructor

	// edges section
	Edge_Structure * edges;
	virtual int get_elements (int k_edge, int * k_element) { return 0; }; // returns array of elements that have the edge k_edge
	virtual void numerate_edges () {}; // numerates edges
	virtual int get_amount_edges (int k_element) { return 0; }; // returns amount of edges of the element
	virtual void get_element_edges (int k_element, int * element_edges) {}; // returns numbers of edges of the element

	// mesh functions section
	virtual int get_n_nodes () const { return 0; };  // returns amount of nodes
	virtual int get_n_elements () const { return 0; };  // returns amount of elements
	virtual bool build_Mesh (char * file_name, bool save, int * N_functions) { return false; }; // makes a mesh
	virtual bool build_Mesh (char * file_name, char * file_name_nodes, char * file_name_elements, int * N_functions) { return false; }; // makes a mesh
	virtual bool build_Mesh (char * file_name_nodes, char * file_name_elements, int * N_functions) { return false; }; // makes a mesh
	virtual void get_0_boundaries (double * c0) const {}; // returns starting points that form the mesh
	virtual void get_N_boundaries (double * cn) const {}; // returns ending points that form the mesh
	virtual void get_area_0_boundaries (int area, double * c0) const {}; // returns starting points that form the mesh
	virtual void get_area_N_boundaries (int area, double * cn) const {}; // returns ending points that form the mesh
	virtual int get_dimentionality () const { return 0; }; // returns mesh's dimentionality
	virtual void wrap_material (int k_material, double * c0, double * cN) {}; // returns a section that contains all the elements with k_material
	virtual int get_amount_of_areas () const { return 0; };
	// node functions section
	virtual void get_node_coordinates (int k_node, double * c) const {}; // returns coordinates of k-node

	// element functions section
	virtual void prepare_element (int k_element) {}; // calls elements[k_element]->prepare
	virtual int get_amount_of_def_nodes (int k_element) const { return 0; }; // returns the size of defining nodes array of k-element
	virtual int get_node_number (int k_element, int i) { return 0; }; // returns i defining node of the k element
	virtual void get_local_function_values (int k_element, double * coordinates, MathVector * v) {}; // return vector of local functions values
	virtual void get_local_function_values (int k_element, double * coordinates, Matrix * v) {}; // return vector of local functions values
	virtual void get_local_function_first_derivative_values (int k_element, int k_var, double * coordinates, MathVector * v) {}; // return vector of local functions derivative values by k_var
	virtual void get_G_local_matrix (int k_element, Matrix * G_matrix) {}; // calls elements[k_element]->get_G_local_matrix
	virtual void get_DDD_local_matrix (int k_element, MathVector * DDD) {}; // calls elements[k_element]->get_DDD
	virtual void get_M_local_matrix (int k_element, Matrix * G_matrix) {}; // calls elements[k_element]->get_M_local_matrix
	virtual void get_edge_functions (int k_element, int n1, int n2, int * functions) {}; // function for local matrix M
	virtual void get_def_nodes (int k_element, int * d_nodes) const {}; // calls elements[k_element]->get_def_nodes
	virtual void get_base_nodes (int k_element, int * b_nodes) const {}; // calls elements[k_element]->get_base_nodes
	virtual int get_amount_of_base_nodes (int k_element) { return 0; }; // calls elements[k_element]->get_n_base_nodes
	virtual int get_area (int k_element) const { return 0; }; // returns k_element's area
	virtual int get_area (double * coord) { return 0; }; // returns area by point
	virtual double get_geometrical_area (int k_element) const  { return 0; }; // returns k_element's area
	virtual int get_amount_non_zero_functions (int k_element) { return 0; }; // returns G and M matrices' dimentionality for k_element
	virtual int belonging_element (int k_node) { return 0; }; // returns first element that has k_node
	virtual int function_boundary (int k_element, int k_func) { return 0; }; // calls elements[k_element]->get_function boundary
	virtual int function_boundary_edge (int n1, int n2) { return 0; }; // calls elements[k_element]->get_function boundary
	virtual int amount_of_integration_points (int k_element) { return 0; }; // calls elements[k_element]->amount_of_integration_points
	virtual void integration_points (int k_element, double ** points, double * weigths, double * jac) {}; // calls elements[k_element]->integration_points
	virtual double get_basis_function_value (int k_element, int k_func, double * coordinates) { return 0.0; }; // calls elements[k_element]->get_basis_function_value
	virtual void get_mass_center (int k_element, double * c) {};
	virtual int get_order (int k_element) { return 0; };
	virtual int get_function_global_number (int k_element, int k_func) { return 0; };
	virtual void get_area_dividers (std::vector <std::pair<int, int>> * dividers_nodes) {}; // returns boundaries between areas
	virtual int get_function_by_node (int k_element, int n) { return -1; };
	virtual int get_function_by_edge (int k_element, int n1, int n2) { return -1; };

	// functions that connect elements and nodes
	virtual bool point_inside (int k_element, double * coordinates) { return false; }; // calls elements[k_element]->point_inside
	virtual int point_inside (double * coordinates) { return 0; }; // returns element that has point
	virtual bool edge_exists (int k_element, int n1, int n2) { return false; }; // calls elements[k_element]->edge_exists 

	// edges section
	virtual int belonging_element (int n1, int n2) { return 0; }; // finds element that has the node
	virtual int get_amount_second_condition (int k_element) { return 0; }; // calls elements[k_element]->get_amount_second_condition 
	
	// functions that tell important things section
	virtual bool need_to_numerate_edges () { return false; }; // returns true, if the task has to numerate edges 

	// for painter
	virtual void fraction (int k) {};
	virtual void copy (const Mesh_Prototype & mesh) {};
	virtual void copy_nodes (const Mesh_Prototype & mesh) {};

	// folds section
	virtual int get_amount_of_folds () { return 0; };
	virtual void get_fold_coordinates (int k_fold, double * c0, double * cN) {};

	// p-refinement
	virtual int refine (std::vector <int> elements_to_refine) { return 0; };

	// for painter
	virtual bool get_isoline_section (int k_element, double * q, double value, double * c1, double * c2) { return false; };

	// for dynamic mesh
	virtual void move_node (int k_node, double * coordinates) {};
	virtual double nodes_distance (int k1_node, int k2_node) { return 0.0; }; // distance between two nodes
	virtual void refresh_mesh () {}; // recalculates all elements' matrices after moving all nodes
};