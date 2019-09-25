#pragma once
#include "Triangle_Element.h"

class Triangle_Hier : public Triangle
{
private:
public:
	Triangle_Hier ();
	Triangle_Hier (const Triangle_Hier & triangle);
	~Triangle_Hier ();

	double get_basis_function_value (int k_func, double * coordinates) override; // returns k_func basis function value 
	double get_basis_function_derivative (int k_func, int k_var, double * coordinates) override; // returns k_func basis derivative by k_var variable  

	virtual int get_function_nodes (int k_func, int * nodes) override; // returns local numbers of nodes where function is not zero
	virtual int get_amount_second_condition () override;  // returns amount of functions that are not 0 on edge
	virtual void get_edge_functions (const Mesh_Prototype & mesh, int n1, int n2, int * functions) override; // returns numbers of functions that are not 0 on the edge
	virtual int function_boundary (const Mesh_Prototype & mesh, int k_func) override; // gets number of border for boundary function

	virtual int get_isoline_points (const Mesh_Prototype & mesh, double value, double * q, double * c1, double * c2) override;

	virtual int get_function_by_node (int n) override;
	virtual int get_function_by_edge (int n1, int n2) override;
};