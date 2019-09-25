#pragma once
#include "Prototype_Element.h"

#define ZERO_cube 1e-10

class Cube : public Element
{
protected:
	double c0[3], cN[3]; // element boundaries
public:
	Cube ();
	Cube (const Cube & cube);
	~Cube ();

	double get_geometrical_area () override; // returns geometrical area
	bool point_inside (const Mesh_Prototype & mesh, double * coordinates) override; // check if point defined by coordinates is in the element
	void inside_prepare (const Mesh_Prototype & mesh) override; // sets Alpha, D matrix
	
	double get_basis_function_value (int k_func, double * coordinates) override; // returns k_func basis function value 
	double get_basis_function_derivative (int k_func, int k_var, double * coordinates) override; // returns k_func basis derivative by k_var variable  
	void prepare_GM () override; // sets private variables, if needed
	int get_amount_non_zero_functions () override; // returns dimentionality for G and M
	void get_boundaries (double * C0, double * CN); // return c0 and cN

	// integrate section
	int amount_of_integration_points () override;
	void integration_points (double ** points, double * weigths, double * jac) override;
};