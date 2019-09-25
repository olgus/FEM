#pragma once
#include "Prototype_Mesh.h"
#include "Triangle_Element.h"

#define ZERO_DIFFERENCE 1e-7

struct mesh_side 
{
	int sections; // amount of sections on each direction
	double coef; // coef on each side
	int direc; // direction of concentration on each side

	std::vector <Node_2D> nodes;
};

// triangular mesh, 1 lvl approx, triangles
class Triangular_Mesh : public Mesh <Node_2D, Triangle>
{
private:
	int lvl; // nesting lvl 0-infinity
		
	// for frontal methods
	mesh_side ** sides; // vertical and horizontal sides
	Point_2D ** tetra_nodes; // nodes that determine tetragonal areas
	int * materials; // materials data
	int n_areas; // amount of tetragonal areas
	int N, M, max; // amount of horizontal sides, vertical and max out of them
	std::vector <Node_2D> t0, t1;

	void input_mesh_data (char * file_name) override; // input mesh data from file
	
	bool make_init_Mesh () override; // set the initial mesh

	void set_triangles (int a); // set initial triangles of the mesh
	void renumerate () override; 
	bool need_to_numerate_edges () override; // returns true for 2D tasks

	// fold section
	std::vector <Fold> folds; 
	virtual void fold_mesh () override; // fold the mesh
public:
	Triangular_Mesh (); // constructor
	Triangular_Mesh (const Triangular_Mesh & triangular_Mesh); // constructor-copy
	~Triangular_Mesh (); // destructor

	int build_mesh (double * C0, double * CN, int * N_axis); // input mesh data
	int build_mesh (double * C0, double * CN, int * N_axis, char * file_nodes, char * file_triangles); // input mesh data
	void fraction (int k) override; 
	void copy (const Mesh_Prototype & mesh) override;

	// refinement section
	int refine (std::vector <int> elements_to_refine) override;

	// point through folds
	using Mesh::point_inside;
	virtual int point_inside (double * coordinates) override; // returns element that has point

	// folds section
	int get_amount_of_folds () override;
	void get_fold_coordinates (int k_fold, double * c0, double * cN) override;
	Triangular_Mesh & operator= (const Triangular_Mesh & mesh); // assignment operator

	void output () override;
	void output (char * file_nodes, char * file_triangles) override;
	bool get_isoline_section (int k_element, double * q, double value, double * c1, double * c2) override;

	void read_nodes (char * file_nodes);
};
