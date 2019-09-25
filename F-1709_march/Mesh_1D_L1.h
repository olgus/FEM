#pragma once
#include "Prototype_Mesh.h"

#define ZERO_mesh_1D_L1 1e-9

class Mesh_1D_L1 : public Mesh<Node_1D, Element_1D_L1>
{
private:
	int lvl; // nesting lvl 0-infinity
	int * materials; // material

	int * n_sections; // amount of sections within sections
	double * coef; // coefficients of the sections
	double * boundaries; // values of sections' boundaries

	int numerate_functions () override; // numerates nodes' functions
	void input_mesh_data (char * file_name) override; // input mesh data from file
	bool make_init_Mesh () override; // set the initial mesh
	virtual bool need_to_numerate_edges () override;
protected:
public:
	using Mesh::build_Mesh;
	// take that away later
	bool build_Mesh (std::vector <double> coordinates);
	void reset_elements () override;

	Mesh_1D_L1 ();
	Mesh_1D_L1 (const Mesh_1D_L1 & mesh);
	~Mesh_1D_L1 ();
	void output () override;
};
