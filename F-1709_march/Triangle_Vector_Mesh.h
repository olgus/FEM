#pragma once
#include "Prototype_Mesh.h"

class Triangle_Vector_Mesh : public Mesh <Node_2D, Triangle_Vector>
{
private:
protected:
public:
	Triangle_Vector_Mesh ();
	Triangle_Vector_Mesh (const Triangle_Vector_Mesh & mesh);
	~Triangle_Vector_Mesh ();

	using Mesh::get_basis_function_value;
	void get_basis_function_value (int k_element, int k_func, double * coordinates, double * f_value);
	int numerate_functions () override; // numerates functions
									 
	// sooo
	// parent's build_mesh function should work fine from files ??
};