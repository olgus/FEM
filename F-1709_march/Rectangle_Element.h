#pragma once
#include "Prototype_Element.h"

class Rectangle_Element : public Element
{
protected:
	double x0, xN, y0, yN;
public:
	Rectangle_Element (); // contructor
	Rectangle_Element (const Rectangle_Element & Rectangle_Element); // constructor-copy
	~Rectangle_Element (); // destructor

	virtual double integrate (int i, int j, int * m); // integrates a function on a triangle, m sets a g(i, j) or m(i, j)
	double get_geometrical_area () override; // returns geometrical area

	bool point_inside (const Mesh_Prototype & mesh, double * coordinates) override; // check if point defined by coordinates is in the element

	void inside_prepare (const Mesh_Prototype & mesh) override; // sets G, M
	
	virtual double get_basis_function_value (int k_func, double * coordinates) override; // returns k_func basis function value 
	virtual double get_basis_function_derivative (int k_func, int k_var, double * coordinates) override; // returns k_func basis derivative by k_var variable  
	virtual double get_function_value (int i, int j, int * m, double * coordinates) override; // returns value of functions for g, m

	// also finish these ones
	int get_amount_second_condition () override; // returns dimentionality for G and M
	void get_edge_functions (const Mesh_Prototype & mesh, int n1, int n2, int * functions) override; // function for local matrix M

	// edges
	bool edge_exists (int n1, int n2) override; // checks if edge between node n1 and n2 exists
	int get_amount_edges () override; // returns amount of edges

	// integrate section
	int amount_of_integration_points () override;
	void integration_points (double ** points, double * weigths, double * jac) override;

	// finish work with edges
	int edge_amount_of_integration_points () override;
	void edge_integration_points (const Mesh_Prototype & mesh, int n1, int n2, double ** points, double * weigths, double * jac) override;
};