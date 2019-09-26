#pragma once
#include "Prototype_Element.h" 

class Element_1D_L1 : public Element
{
private:
protected:
	double x0, xN;
public:
	Element_1D_L1 ();
	Element_1D_L1 (const Element_1D_L1 & e);
	virtual ~Element_1D_L1 ();

	virtual double get_geometrical_area () override; // returns geometrical area
	virtual bool point_inside (const Mesh_Prototype & mesh, double * coordinates) override; // check if point defined by coordinates is in the element

	virtual void inside_prepare (const Mesh_Prototype & mesh) override;

	virtual double integrate (int i, int j, int * m) override; // integrates a function on a section, m sets a g(i, j), m(i, j) 
	
	// functions section
	virtual double get_basis_function_value (int k_func, double * coordinates) override; // returns k_func basis function value 
	virtual double get_basis_function_derivative (int k_func, int k_var, double * coordinates) override; // returns k_func basis derivative by k_var variable  

	virtual int get_isoline_points (const Mesh_Prototype & mesh, double value, double * q, double * c1, double * c2) override;

	// integrate section
	virtual int amount_of_integration_points () override;
	virtual void integration_points (double ** points, double * weigths, double * jac) override;
};
