#pragma once
#include "Prototype_Mesh.h"
#include "Rectangle_Element_3.h"

#define ZERO_rectangular_mesh_3 1e-5

class Rectangular_Mesh_3 : public Mesh <Node_2D, Rectangle_Element_3>
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
	Rectangular_Mesh_3 ();
	Rectangular_Mesh_3 (const Rectangular_Mesh_3 & rm3);
	~Rectangular_Mesh_3 ();
	void output () override;
};