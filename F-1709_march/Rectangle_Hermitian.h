#pragma once
#include "Rectangle_Element.h"

#define ZERO_Rectangle_Element_hermitian 1e-8

class Rectangle_Element_Hermitian : public Rectangle_Element
{
private:
	Matrix * D; // laplasian * laplasian
public:
	Rectangle_Element_Hermitian ();
	Rectangle_Element_Hermitian (const Rectangle_Element_Hermitian & rect3);
	~Rectangle_Element_Hermitian ();

	virtual void set_global_functions () override; // set element's nodes and area

	bool point_inside (const Mesh_Prototype & mesh, double * coordinates) override; // check if point defined by coordinates is in the element

	void get_D_local_matrix (Matrix * D_matrix); // function for local matrix D

	double get_basis_function_value (int k_func, double * coordinates) override; // returns k_func basis function value 
	double get_basis_function_derivative (int k_func, int k_var, double * coordinates) override; // returns k_func basis derivative by k_var variable 
	double get_basis_function_derivative_second (int k_func, int k_var, double * coordinates); // returns k_func basis derivative by k_var variable 
	double get_function_value (int i, int j, int * m, double * coordinates) override; // returns value of functions for g, m

	int get_amount_non_zero_functions () override; // return n_functions
	void inside_prepare (const Mesh_Prototype & mesh) override; // sets A, G, D
};