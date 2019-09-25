#pragma once
#include "Triangular_Mesh.h"
#include "Mesh_1D_L1.h"
#include "Prism.h"
#include <string.h>

#define ZERO_DIFFERENCE_prism 1e-7

class Prismatic_Mesh : public Mesh <Node_3D, Prism>
{
private:
	Triangular_Mesh * mesh_triangular;
	Mesh_1D_L1 * mesh_1d;
	// separate the file_name into + XY + Z
	void input_mesh_data (char * file_name) override; // input mesh data from file, file_name is a template
	bool make_init_Mesh () override; // set the initial mesh
	bool need_to_numerate_edges () override; // returns true for 2D tasks
public:
	Prismatic_Mesh (); // constructor
	Prismatic_Mesh (const Prismatic_Mesh & prismatic_Mesh); // constructor-copy
	~Prismatic_Mesh (); // destructor
	void output () override;
};