#pragma once
#include "Triangle_Element.h"

class Triangle_Polar : public Triangle
{
public:
	Triangle_Polar (); // contructor
	Triangle_Polar (const Triangle_Polar & triangle); // constructor-copy
	~Triangle_Polar (); // destructor

	double get_function_value (int i, int j, int * m, double * coordinates) override; // returns value of functions for g, m
};