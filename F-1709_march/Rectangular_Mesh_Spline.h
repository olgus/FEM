#pragma once
#include "Prototype_Mesh.h"
#include "Rectangle_Hermitian.h"

#define ZERO_rectangular_mesh_spline 1e-5

// triangular mesh, 1 lvl approx, triangles
class Rectangular_Mesh_Spline : public Mesh <Node_2D, Rectangle_Element_Hermitian>
{
private:
	// uses n_axis as amount of nodes on axises
	int lvl; // nesting lvl 0-infinity
	Point_2D ** tetra_nodes; // nodes that determine the mesh
	int material; // material
	double coef[2]; // coefficients on axises
	int direc[2]; // direction on axises

	void input_mesh_data (char * file_name) override; // input mesh data from file
	void input_mesh_data (double * c0, double * cN, int * N_axis); // input mesh data
	bool make_init_Mesh () override; // set the initial mesh
	int numerate_functions () override; // numerates nodes' functions
public:
	Rectangular_Mesh_Spline (); // constructor
	Rectangular_Mesh_Spline (const Rectangular_Mesh_Spline & rectangular_mesh); // constructor-copy
	~Rectangular_Mesh_Spline (); // destructor

	using Mesh::build_Mesh;
	virtual bool build_Mesh (double * c0, double * cN, int * N_axis, int * N_functions); // makes a mesh
	void output () override;
};