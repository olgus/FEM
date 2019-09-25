#pragma once
#include "Triangle_Element.h"

class Triangle_Vector : public Triangle
{
private:
	double get_scalar_product (int k1, int k2, double * coordinates);
protected:
public:
	Triangle_Vector ();
	Triangle_Vector (const Triangle_Vector & triangle);
	~Triangle_Vector ();

	// basis functions
	double get_basis_function_value (int k_func, double * coordinates) override; // returns nothing

	void get_basis_function_value (int k_func, double * coordinates, double * f_value) override; // returns k_func basis function value 
	double get_basis_function_derivative (int k_func, int k_var, double * coordinates) override; // returns rot 

	double get_function_value (int i, int j, int * m, double * coordinates) override; // returns value of functions for g, m

	void get_edge_functions (const Mesh_Prototype & mesh, int n1, int n2, int * functions) override; // returns numbers of functions that are not 0 on the edge
};