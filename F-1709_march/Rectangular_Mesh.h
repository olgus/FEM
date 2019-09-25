#pragma once
#include "Prototype_Mesh.h"
#include "Rectangle_Element.h"

#define ZERO_rectangular_mesh 1e-5

// triangular mesh, 1 lvl approx, triangles
class Rectangular_Mesh : public Mesh <Node_2D, Rectangle_Element>
{
private:
	// uses n_axis as amount of nodes on axises
	int lvl; // nesting lvl 0-infinity
	Point_2D ** tetra_nodes; // nodes that determine the mesh
	int material; // material
	double coef[2]; // coefficients on axises
	int direc[2]; // direction on axises

	void input_mesh_data (char * file_name) override; // input mesh data from file
	bool make_init_Mesh () override; // set the initial mesh
public:
	Rectangular_Mesh (); // constructor
	Rectangular_Mesh (const Rectangular_Mesh & rectangular_mesh); // constructor-copy
	~Rectangular_Mesh (); // destructor

	void make_uniform_mesh (int n1, int n2, double * c0, double * cN, int Lvl);
	void output () override;
};