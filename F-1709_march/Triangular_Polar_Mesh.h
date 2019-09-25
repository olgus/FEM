#pragma once
#include "Prototype_Mesh.h"
#include "Triangle_Element_Polar.h"
#include "Node_2D_Polar.h"

#define ZERO_DIFFERENCE_tr_polar 1e-7

struct mesh_side_polar
{
	int sections; // amount of sections on each direction
	double coef; // coef on each side
	int direc; // direction of concentration on each side

	std::vector <Node_2D_Polar> nodes;
};

// triangular mesh, 1 lvl approx, triangles
class Triangular_Polar_Mesh : public Mesh <Node_2D_Polar, Triangle_Polar>
{
private:
	int lvl; // nesting lvl 0-infinity
	mesh_side_polar ** sides; // vertical and horizontal sides
	Point_2D_Polar ** tetra_nodes; // nodes that determine tetragonal areas
	int * materials; // materials data
	int n_areas; // amount of tetragonal areas
	int N, M, max; // amount of horizontal sides, vertical and max out of them

	std::vector <Node_2D_Polar> t0, t1;

	void input_mesh_data (char * file_name) override; // input mesh data from file
	bool make_init_Mesh () override; // set the initial mesh
	void set_triangles (int a); // set initial triangles of the mesh
	void renumerate () override; 
public:
	Triangular_Polar_Mesh (); // constructor
	Triangular_Polar_Mesh (const Triangular_Polar_Mesh & triangular_Mesh); // constructor-copy
	~Triangular_Polar_Mesh (); // destructor

	void output () override;
	bool get_isoline_section (int k_element, double * q, double value, double * c1, double * c2) override;
};