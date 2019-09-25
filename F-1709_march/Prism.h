#pragma once
#include "Prototype_Element.h"

#define ZERO_prism 1e-10

class Prism : public Element
{
protected:
	// for triangular faces
	Matrix * Alpha; // alpha-coordinates of the triangle. al(i) is a row
	Matrix * D; // D matrix
	double determinant; // determinant of D matrix

	// height
	double z0, zN;
public:
	Prism ();
	Prism (const Prism & prism);
	~Prism ();

	double get_geometrical_area () override; // returns geometrical area
	bool point_inside (const Mesh_Prototype & mesh, double * coordinates) override; // check if point defined by coordinates is in the element
	void inside_prepare (const Mesh_Prototype & mesh) override; // sets Alpha, D matrix

	double get_basis_function_value (int k_func, double * coordinates) override; // returns k_func basis function value 
	double get_basis_function_derivative (int k_func, int k_var, double * coordinates) override; // returns k_func basis derivative by k_var variable  
	void prepare_GM () override; // sets private variables, if needed
	int get_amount_non_zero_functions () override; // returns dimentionality for G and M
	void get_D (int var_der, Matrix * D) override; // returns integrate (bfi * d bfj/ d var_der)
};