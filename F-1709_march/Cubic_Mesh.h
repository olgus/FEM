#pragma once
#include "Prototype_Mesh.h"
#include "Cube_Element.h"
#include <string.h>

class Cubic_Mesh : public Mesh <Node_3D, Cube>
{
private:
	void input_mesh_data (char * file_name) override; // input mesh data from file, file_name is a template
	bool make_init_Mesh () override; // set the initial mesh
	bool need_to_numerate_edges () override; // returns true for 2D tasks
public:
	using Mesh::build_Mesh;
	virtual bool build_Mesh (double * C0, double * CN, int * N_axis); // makes a mesh
	void reset_areas (char * file_areas); // resets area values individually 

	Cubic_Mesh (); // constructor
	Cubic_Mesh (const Cubic_Mesh & cubic_Mesh); // constructor-copy
	~Cubic_Mesh (); // destructor

	void output () override;
	bool check_neighbouring_state (int e1, int e2); // check if e1 and e2 are neighbours
};