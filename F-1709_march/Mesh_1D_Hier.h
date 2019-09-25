#pragma once
#include "Prototype_Mesh.h"

#define ZERO_mesh_1D 1e-9

class Mesh_1D_Hier : public Mesh<Node_1D, Element_1D_hier>
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
public:
	using Mesh::build_Mesh;
	bool build_Mesh (std::vector <double> coordinates);

	Mesh_1D_Hier ();
	Mesh_1D_Hier (const Mesh_1D_Hier & mesh);
	~Mesh_1D_Hier ();
	void output () override;
};