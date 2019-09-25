#pragma once
#include "Prototype_Mesh.h"
#include "Triangle_Element.h"
#include "Triangle_Element_Hier.h"

#define ZERO_DIFFERENCE_Mixed_Mesh 1e-7

struct mesh_side_MM
{
	int sections; // amount of sections on each direction
	double coef; // coef on each side
	int direc; // direction of concentration on each side

	std::vector <Node_2D> nodes;
};

// triangular mesh, 1 lvl approx, triangles
class Mixed_Triangular_Mesh : public Mesh <Node_2D, Element>
{
private:
	int lvl; // nesting lvl 0-infinity
	mesh_side_MM ** sides; // vertical and horizontal sides
	int * order; // order of the elements within the section
	Point_2D ** tetra_nodes; // nodes that determine tetragonal areas
	int * materials; // materials data
	int n_areas; // amount of tetragonal areas
	int N, M, max; // amount of horizontal sides, vertical and max out of them

	std::vector <Node_2D> t0, t1;

	void input_mesh_data (char * file_name) override; // input mesh data from file
	bool make_init_Mesh () override; // set the initial mesh
	void set_triangles (int a, int o); // set initial triangles of the mesh
	void renumerate () override;
	bool need_to_numerate_edges () override; // returns true for 2D tasks
	int numerate_functions () override; // numerates nodes' functions

	// fold section
	std::vector <Fold> folds;
	virtual void fold_mesh () override; // fold the mesh
public:
	Mixed_Triangular_Mesh (); // constructor
	Mixed_Triangular_Mesh (const Mixed_Triangular_Mesh & mesh); // constructor-copy
	~Mixed_Triangular_Mesh (); // destructor

	void build_mesh (double * C0, double * CN, int * N_axis); // input mesh data

	using Mesh::build_Mesh;
	virtual bool build_Mesh (char * file_name_nodes, char * file_name_elements, int * N_functions) override; // makes a mesh
	void wrap_material (int k_material, double * c0, double * cN) override; // returns a section that contains all the elements with k_material

	// point through folds
	using Mesh::point_inside;
	virtual int point_inside (double * coordinates) override; // returns element that has point
															  // folds section
	int get_amount_of_folds () override;
	void get_fold_coordinates (int k_fold, double * c0, double * cN) override;
	void output () override;

	bool get_isoline_section (int k_element, double * q, double value, double * c1, double * c2) override;
};
