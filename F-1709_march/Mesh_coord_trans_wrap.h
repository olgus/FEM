#pragma once
#include "Prototype_Mesh.h"

/* The idea is that that wrapper should contain:
- fixed computational mesh with nodes and elements
- transformations that hide the computational mesh and give out physical values
	- for boundaries 
	- for nodes
	- for searching the element with a point
	- local function values
	- integration points
	- ... basically for any function that deals with straightforward coordinates
*/