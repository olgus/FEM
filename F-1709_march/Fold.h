#pragma once
#include <vector>

struct Fold
{
public:
	double c0[2], cN[2]; // corners of the fold
	std::vector <int> elements; // numbers of elements that sit in them
	int k_axis; // axis that was used to create a fold
}; 