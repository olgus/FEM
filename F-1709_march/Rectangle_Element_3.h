#pragma once
#include "Rectangle_Element.h"

class Rectangle_Element_3 : public Rectangle_Element
{
private:

public:
	Rectangle_Element_3 ();
	Rectangle_Element_3 (const Rectangle_Element_3 & rect3);
	~Rectangle_Element_3 ();

	virtual double get_basis_function_value (int k_func, double * coordinates) override; // returns k_func basis function value 
	virtual double get_basis_function_derivative (int k_func, int k_var, double * coordinates) override; // returns k_func basis derivative by k_var variable 
};