#pragma once
#include "Prototype_Mesh.h"
#include "Triangle_Element_Hier.h"
#include "Node_2D.h"

#define ZERO_DIFFERENCE_mesh_triangle_hier 1e-5

struct mesh_side_hier
{
	int sections; // amount of sections on each direction
	double coef; // coef on each side
	int direc; // direction of concentration on each side

	std::vector <Node_2D> nodes;
};

class Triangular_Mesh_Hier : public Mesh<Node_2D, Triangle_Hier>
{
private:
	int lvl; // nesting lvl 0-infinity
	mesh_side_hier ** sides; // vertical and horizontal sides
	Point_2D ** tetra_nodes; // nodes that determine tetragonal areas
	int * materials; // materials data
	int n_areas; // amount of tetragonal areas
	int N, M, max; // amount of horizontal sides, vertical and max out of them

	std::vector <Node_2D> t0, t1;

	void input_mesh_data (char * file_name) override; // input mesh data from file
	bool make_init_Mesh () override; // set the initial mesh
	void set_triangles (int a); // set initial triangles of the mesh
	void renumerate () override;
	int numerate_functions () override; // numerates nodes' functions

	// fold section
	std::vector <Fold> folds;
	virtual void fold_mesh () override; // fold the mesh
public:
	Triangular_Mesh_Hier ();
	Triangular_Mesh_Hier (const Triangular_Mesh_Hier & tmh);
	~Triangular_Mesh_Hier ();

	// point through folds
	using Mesh::point_inside;
	virtual int point_inside (double * coordinates) override; // returns element that has point
	
	// folds section
	int get_amount_of_folds () override;
	void get_fold_coordinates (int k_fold, double * c0, double * cN) override;

	void output () override;
	bool get_isoline_section (int k_element, double * q, double value, double * c1, double * c2) override;
};