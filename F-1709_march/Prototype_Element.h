#pragma once
#include <stdio.h>
#include <algorithm>

#include "Matrix.h"
#include "Mesh_Pointer.h"

#define ZERO_element 1e-10
#define ZERO_boundary 1e-10

// element prototype
class Element
{
private:
protected:
	int * defining_nodes; // numbers of nodes defining the element
	int * base_nodes; // numbers of nodes that define the element, the same as defining_nodes for 1-lvl element
	int * global_functions; // numbers of global functions on the element

	int n_def_nodes; // amount of defining nodes
	int n_base_nodes; // amount of base nodes
	int area; // area number
	int order; // element's order
	int n_functions; // amount of functions on the element 
	int dim; // node_dimentionality
	
	Matrix * G; // G-matrix
	Matrix * M; // M-matrix

public:
	Element (); // constructor
	Element (const Element & element); // copy-constructor
	virtual ~Element (); // destructor
	
	virtual void set_amount_def_nodes (int n); // set n_def_nodes
	virtual void set_element (int * nodes); // set element's nodes and area
	virtual void set_base_nodes (int * nodes); // set element's nodes and area
	virtual void set_nodes (int * nodes); // set element's nodes
	virtual void set_element (int n, int * nodes); // set n_def_nodes and element's nodes
	virtual void set_area (int a); // set node's area
	virtual void set_global_functions (int * funcs); // sets global numbers of functions
	virtual void set_global_functions (); // sets global numbers by def_nodes
	virtual void set_global_function (int k_func_local, int k_func_global); // set k_func_local

	virtual double integrate (int i, int j, int * m); // integrates a function on a triangle, m sets a g(i, j) or m(i, j)
	
	void get_mass_center (const Mesh_Prototype & mesh, double * coordinates); // return coordinates of mass center of the element

	int get_amount_of_def_nodes () const; // get n_def_nodes
	int get_amount_of_base_nodes () const; // get n_base_nodes
	void get_def_nodes (int * n); // get defining_nodes 
	int get_def_node (int k_node); // returns global of defining k_node
	void get_base_nodes (int * n); // get base_nodes 
	int get_area () const; // set node's area
	int get_node (int k_node) const; // return global number of k-defining node
	bool node_belongs (int k_node); // true, if an element has k_node
	void sort_def_nodes (); // sort defining nodes in increasing order
	virtual int get_function_global_number (int k_func); // returns global number of k_func local function
	int get_order ();

	virtual int get_amount_second_condition (); // returns amount of functions that are not 0 on the edge
	virtual double get_geometrical_area (); // returns geometrical area of k-element
	virtual void get_G_local_matrix (Matrix * G_matrix); // function for local matrix G 
	virtual void get_M_local_matrix (Matrix * M_matrix); // function for local matrix M
	virtual void get_edge_functions (const Mesh_Prototype & mesh, int n1, int n2, int * functions); // returns numbers of functions that are not 0 on the edge
	
	virtual bool point_inside (const Mesh_Prototype & mesh, double * coordinates); // check if point defined by coordinates is in the element
	virtual void get_local_function_values (double * coordinates, MathVector * v); // returns local functions' values in point, defined by coordinates
	virtual void get_local_function_values (double * coordinates, Matrix * v); // returns local functions' values in point, defined by coordinates
	virtual void get_local_function_first_derivative_values (double * coordinates, int k_var, MathVector * v); // returns local functions derivatives' values in point, defined by coordinates
	virtual void prepare (const Mesh_Prototype & mesh); // prepares element (sets private variables, if needed), sets G, M
	virtual void inside_prepare (const Mesh_Prototype & mesh); // sets G, M
	virtual void prepare_GM (); // sets private variables, if needed

	virtual double get_function_value (int i, int j, int * m, double * coordinates); // returns value of functions for g, m
	virtual double get_basis_function_value (int k_func, double * coordinates); // returns k_func basis function value 
	virtual void get_basis_function_value (int k_func, double * coordinates, double * f_value); // returns k_func basis function value 
	virtual double get_basis_function_derivative (int k_func, int k_var, double * coordinates); // returns k_func basis derivative by k_var variable  
	virtual int function_boundary (const Mesh_Prototype & mesh, int k_func); // gets number of border for boundary function
	virtual int get_amount_non_zero_functions (); // returns dimentionality for G and M
	virtual bool edge_exists (int n1, int n2); // checks if edge between node n1 and n2 exists

	virtual int get_function_nodes (int k_func, int * nodes); // returns local numbers of nodes where function is not zero

	virtual void get_edges (int ** edges); // returns nodes of edges
	virtual int get_amount_edges (); // returns amount of edges

	virtual void get_D (int var_der, Matrix * D); // returns integrate (bfi * d bfj/ d var_der)
	virtual void get_DDD (MathVector * D); // returns integrate (bfi * bfj * bfk)
	virtual void get_DDF (int var_der1, int var_der2, MathVector * D); //  returns integrate (d bfk/ d var_der1 * d bfj/ d var_der2 * bfi)
	virtual void get_DF (int var_der, Matrix * D); // returns integrate (d bfj/ d var_der * bfi)
	virtual void get_DDF (int k, int var_der1, int var_der2, Matrix * D); // returns integrate (d bfk/ d var_der1 * d bfj/ d var_der2 * bfi)

	// integrate section
	virtual int amount_of_integration_points ();
	virtual void integration_points (double ** points, double * weigths, double * jac);
	
	// for second boundary condition on element's edge 
	virtual int edge_amount_of_integration_points ();
	virtual void edge_integration_points (const Mesh_Prototype & mesh, int n1, int n2, double ** points, double * weigths, double * jac);
	virtual void segment_integration_points (const Mesh_Prototype & mesh, double * coordinates_1, double * coordinates_2, double ** points, double * weigths, double * jac);

	virtual int get_isoline_points (const Mesh_Prototype & mesh, double value, double * q, double * c1, double * c2);

	virtual int get_function_by_node (int n);	
	virtual int get_function_by_edge (int n1, int n2);
	// for dynamic mesh
	virtual void recalc_matrices (); // recalc G and M after moving the node
};

// functions to change

// for new functions:
// get_basis_function_value
// get_basis_function_derivative
// possibly get_function_nodes
// possibly integrate_function_MB
// get_edge_functions
// possibly function_boundary
// possibly get_amount_second_condition

// for new geometry:
// prepare
// point_inside
// get_geometrical_area
// get_function_value
// edge_exists
// get_amount_non_zero_functions
// edge_integrate
// get_edges
// get_amount_edges
// integrate_function_MB
// integrate