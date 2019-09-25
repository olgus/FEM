#pragma once
#include "Prototype_Element.h"

class Triangle : public Element
{
protected:
	Matrix * Alpha; // alpha-coordinates of the triangle. al(i) is a row
	Matrix * Determinant_M; // Determinant_M matrix
	double determinant; // determinant of D matrix
	double length[3]; // edges' lengths

	Matrix * D; // matrix for (bfi * d bfj/ d var_der)
	MathVector * DDx; // matrix for (bfi * d bfj / d x * d bfk / d y)
	MathVector * DDy; // matrix for (bfi * d bfj / d y * d bfk / d x)
public:
	Triangle (); // contructor
	Triangle (const Triangle & triangle); // constructor-copy
	~Triangle (); // destructor
	
	double get_geometrical_area () override; // returns geometrical area
	
	bool point_inside (const Mesh_Prototype & mesh, double * coordinates) override; // check if point defined by coordinates is in the element
	void inside_prepare (const Mesh_Prototype & mesh) override; // sets Alpha, D, G, M matrix
	
	double get_basis_function_value (int k_func, double * coordinates) override; // returns k_func basis function value 
	double get_basis_function_derivative (int k_func, int k_var, double * coordinates) override; // returns k_func basis derivative by k_var variable  
	bool edge_exists (int n1, int n2) override; // checks if edge between node n1 and n2 exists
	virtual int get_amount_second_condition () override; // returns amount of functions that are not 0 on edge
	virtual void get_edge_functions (const Mesh_Prototype & mesh, int n1, int n2, int * functions) override; // returns numbers of functions that are not 0 on the edge

	virtual void get_edges (int ** edges) override; // returns nodes of edges
	virtual int get_amount_edges () override; // returns nodes of edges
	virtual void get_D (int var_der, Matrix * D_copy) override; // returns integrate (bfi * d bfj/ d var_der)
	virtual void get_DDF (int var_der1, int var_der2, MathVector * DDvd) override; //  returns integrate (d bfk/ d var_der1 * d bfj/ d var_der2 * bfi)
	virtual void get_DDF (int k, int var_der1, int var_der2, Matrix * D_copy) override; // returns integrate (d bfk/ d var_der1 * d bfj/ d var_der2 * bfi)

	// integrate section
	double integrate (int i, int j, int * m) override; // integrates a function on a triangle, m sets a g(i, j), m(i, j) 
	int amount_of_integration_points () override;
	void integration_points (double ** points, double * weigths, double * jac) override;
	void edge_integration_points (const Mesh_Prototype & mesh, int n1, int n2, double ** points, double * weigths, double * jac) override;
	void segment_integration_points (const Mesh_Prototype & mesh, double * coordinates_1, double * coordinates_2, double ** points, double * weigths, double * jac) override;
	int edge_amount_of_integration_points () override;

	virtual int get_isoline_points (const Mesh_Prototype & mesh, double value, double * q, double * c1, double * c2) override;

	virtual int get_function_by_node (int n) override;
	virtual int get_function_by_edge (int n1, int n2) override;
};