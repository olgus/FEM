#pragma once
#include "Element_1D_L1.h"

class Element_1D_hier : public Element_1D_L1 
{
private:
protected:
public:
	Element_1D_hier ();
	Element_1D_hier (const Element_1D_hier & e);
	~Element_1D_hier ();
	
	virtual double get_basis_function_value (int k_func, double * coordinates) override; // returns k_func basis function value 
	virtual double get_basis_function_derivative (int k_func, int k_var, double * coordinates) override; // returns k_func basis derivative by k_var variable  
};