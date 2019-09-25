#pragma once
#include <stdio.h>
#include <vector>
#include <queue>
#include <memory>
#include <set>

#include "Mesh_Pointer.h"
#include "Node_Types.h"
#include "Element_Types.h"
#include "Fold.h"

#define ZERO_prototype_mesh 1e-10
#define ZERO_boundary_edge 1e-10

// mesh prototype
template <class Node_Type, class Element_Type>
class Mesh : public Mesh_Prototype
{
private:

protected:
	// mesh section
	int n_nodes; // amount of nodes
	int n_elements; // amount of elements
	int n_material_areas; // amount of areas;

	// area section
	int dim; // dimensionality
	double * coord0, *coordN; // 0 and N coordinates defining the mesh
	int * n_axis; // amount of defining nodes on axis

	// building functions section
	virtual void input_mesh_data (char * file_name); // input mesh data from file
	virtual void input_mesh_data (char * file_name, char * file_density); // input mesh data from file

	virtual bool make_init_Mesh (); // set the initial mesh
	virtual bool make_init_Mesh (char * file_density); // set the initial mesh

	virtual void renumerate (); // renumerates the mesh
	virtual int numerate_functions (); // numerates nodes' functions
	virtual void fold_mesh (); // fold the mesh

	// support functions section
	virtual int find_maximal_node (); // finds node furthest from the origin
public:
	// mesh data
	std::vector <std::unique_ptr<Element_Type>> elements; // mesh elements
	std::vector <Node_Type> nodes; // mesh nodes

	Mesh (); // constructor
	Mesh (const Mesh & mesh); // copy-constructor
	virtual ~Mesh (); // destructor
		
	// mesh functions section
	virtual void get_0_boundaries (double * boundaries) const override; // returns starting points that form the mesh
	virtual void get_N_boundaries (double * boundaries) const override; // returns ending points that form the mesh
	virtual void get_area_0_boundaries (int area, double * c0) const override; // returns starting points that form the mesh
	virtual void get_area_N_boundaries (int area, double * cn) const override; // returns ending points that form the mesh
	virtual int get_amount_of_areas () const override;
	virtual void get_n_axis (int * axis) const; // returns amount of nodes on axis
	virtual int get_n_nodes () const override; // returns amount of nodes
	virtual int get_n_elements () const override; // returns amount of elements
	virtual void get_area_dividers (std::vector <std::pair<int, int>> * dividers_nodes) override; // returns boundaries between areas
	virtual void wrap_material (int k_material, double * c0, double * cN) override; // returns a section that contains all the elements with k_material

	// element functions section
	int get_amount_of_def_nodes (int k_element) const override; // returns the size of defining nodes array of k-element

	// node functions section
	virtual int get_node_number (Node_Type node); // returns the number of the node
	void get_node_coordinates (int k_node, double * c) const override; // returns coordinates of k-node

	// edges section
	virtual int function_boundary_edge (int n1, int n2) override; // calls elements[k_element]->get_function boundary
	virtual int belonging_element (int n1, int n2) override; // finds element that has the node
	virtual void numerate_edges () override; // numerates edges
	virtual int get_elements (int k_edge, int * k_element) override; // returns array of elements that have the edge k_edge
	virtual int get_amount_edges (int k_element) override; // returns amount of edges of the element
	virtual void get_element_edges (int k_element, int * element_edges) override; // returns numbers of edges of the element
	virtual bool need_to_numerate_edges () override; // returns true for 2D tasks

	virtual void get_elements_by_node (int node, int * n_elements);
	// functions that connect nodes and elements section
	virtual int get_node_number (int k_element, int i) override; // returns i defining node of the k element
	virtual int belonging_element (int k_node) override; // returns first element that has the node
	virtual int function_boundary (int k_element, int k_func) override; // calls elements[k_element]->get_function boundary
	virtual int point_inside (double * coordinates) override; // returns element that has point
	
	// folds section
	virtual int get_amount_of_folds () override ;
	virtual void get_fold_coordinates (int k_fold, double * c0, double * cN) override;

	// refinement section
	virtual int refine (std::vector <int> elements_to_refine) override;

	// inside area section
	virtual bool element_inside_area (int k_element, double * c0, double * cN);

	// functions that call element's functions
	virtual void prepare_element (int k_element) override; // calls elements[k_element]->prepare
	virtual bool point_inside (int k_element, double * coordinates) override; // calls elements[k_element]->point_inside
	virtual void get_local_function_values (int k_element, double * coordinates, MathVector * v) override; // calls elements[k_element]->point_inside
	virtual void get_local_function_values (int k_element, double * coordinates, Matrix * v) override; // return vector of local functions values
	virtual void get_local_function_first_derivative_values (int k_element, int k_var, double * coordinates, MathVector * v) override; // return vector of local functions derivative values by k_var
	virtual void get_G_local_matrix (int k_element, Matrix * G_matrix) override; // calls elements[k_element]->get_G_local_matrix
	virtual void get_M_local_matrix (int k_element, Matrix * G_matrix) override; // calls elements[k_element]->get_M_local_matrix
	virtual void get_DDD_local_matrix (int k_element, MathVector * DDD) override; // calls elements[k_element]->get_DDD
	virtual void get_edge_functions (int k_element, int n1, int n2, int * functions) override; // calls elements[k_element]->get_Theta
	virtual void get_def_nodes (int k_element, int * d_nodes) const override; // calls elements[k_element]->get_def_nodes
	virtual void get_base_nodes (int k_element, int * b_nodes) const override; // calls elements[k_element]->get_base_nodes
	virtual int get_amount_of_base_nodes (int k_element) override; // calls elements[k_element]->get_n_base_nodes};
	virtual int get_amount_non_zero_functions (int k_element) override; // returns G and M matrices' dimentionality for k_element
	virtual int get_area (int k_element) const override; // returns k_element's area
	virtual int get_area (double * coord) override; // returns area by point
	virtual double get_geometrical_area (int k_element) const override; // returns k_element's geometrical area
	virtual int get_amount_second_condition (int k_element) override; // calls elements[k_element]->get_amount_second_condition 
	virtual int amount_of_integration_points (int k_element) override; // calls elements[k_element]->amount_of_integration_points
	virtual void integration_points (int k_element, double ** points, double * weigths, double * jac) override; // calls elements[k_element]->integration_points
	virtual double get_basis_function_value (int k_element, int k_func, double * coordinates) override; // calls elements[k_element]->get_basis_function_value
	virtual void get_mass_center (int k_element, double * c);
	virtual int get_order (int k_element) override;
	virtual int get_function_global_number (int k_element, int k_func) override;
	virtual int get_function_by_node (int k_element, int n) override;
	virtual int get_function_by_edge (int k_element, int n1, int n2) override;

	// external functions section
	virtual bool build_Mesh (char * file_name, bool save, int * N_functions) override; // makes a mesh
	virtual bool build_Mesh (char * file_name, char * file_name_nodes, char * file_name_elements, int * N_functions) override; // makes a mesh
	virtual bool build_Mesh (char * file_name_nodes, char * file_name_elements, int * N_functions) override; // makes a mesh
	
	// input-output section
	virtual void output (); // prints mesh data into log file and mesh files
	virtual void output (char * file_nodes, char * file_triangles); // prints mesh data into log file and mesh files

	Mesh<Node_Type, Element_Type> & operator= (const Mesh<Node_Type, Element_Type> & mesh); // assignment operator

	virtual int get_dimentionality () const override; // returns dim
	virtual bool edge_exists (int k_element, int n1, int n2) override; // calls elements[k_element]->edge_exists 	

	// for painter
	virtual bool get_isoline_section (int k_element, double * q, double value, double * c1, double * c2) override;

	// for dynamic mesh
	virtual void move_node (int k_node, double * coordinates) override; // move node and call element functions for affected elements
	virtual double nodes_distance (int k1_node, int k2_node) override; // distance between two nodes
	virtual void refresh_mesh () override; // recalculates all elements' matrices after moving all nodes

	virtual void reset_n_nodes ();
	virtual void reset_n_elements ();

	virtual void copy_nodes (const Mesh_Prototype & mesh);
	virtual void reset_elements ();
};

// for new mesh override:
// input_mesh_data 
// make_init_Mesh 
// output 
// need_to_numerate_edges
// renumerate if needed