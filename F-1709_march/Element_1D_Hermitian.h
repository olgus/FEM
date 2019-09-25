#pragma once
#include "Element_1D_L1.h" 

class Element_1D_Hermitian : public Element_1D_L1
{
private:
	Matrix * D; // laplasian * laplasian
protected:
public:
	Element_1D_Hermitian ();
	Element_1D_Hermitian (const Element_1D_Hermitian & e);
	~Element_1D_Hermitian ();

	void get_D_local_matrix (Matrix * D_matrix); // function for local matrix D
	void inside_prepare (const Mesh_Prototype & mesh) override; // sets A, G, D

	double get_basis_function_value (int k_func, double * coordinates) override; // returns k_func basis function value 
	double get_basis_function_derivative (int k_func, int k_var, double * coordinates) override; // returns k_func basis derivative by k_var variable  
}; 