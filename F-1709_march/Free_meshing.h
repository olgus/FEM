#pragma once
#include <vector>
#include <stack>
#include "Spline.h"
#include "Triangle_Element.h"
#include "Region.h"

#define PI_MESH 3.1415926535897932384626433832795
#define SIXTYDEGREESRADIAN PI_MESH / 3.0
#define MAX_MESH_RELAXATION_ITER 50

class Free_mesh
{
private:
	double lvl;
	double alpha, beta;
	Spline * spline;
	std::vector<Node_2D> nodes;
	std::vector<Triangle> elements;
	std::vector<int> boundary;

	void input (char * file_mesh, char * file_density, char * file_order);
	void output (char * file_nodes, char * file_elements);

	void build ();
	void relax ();
	void renumerate ();
	
	std::stack <Region> regions;
	void make_inside_nodes (double * c0, double * cN, int * n1, int * n2);
	double calc_angle (double * c, double * c0, double * c1);

	void star (std::vector<std::vector<int>> * star_nodes);
	double distance (int n1, int n2);
	double density_function_value (double * c);
	int find_maximal_node ();

	// order section
	int N_or;
	double ** or_c0;
	double ** or_cN;
	bool point_in_order_region (int k_region, double * coordinates);
	bool increase_order (int k_element);
public:
	Free_mesh ();
	Free_mesh (const Free_mesh & free_mesh);
	~Free_mesh ();

	void build_Mesh (char * file_mesh, char * file_density, char * file_order, char * file_nodes, char * file_elements);
	void build_Mesh_wo_renumerate (char * file_mesh, char * file_density, char * file_order, char * file_nodes, char * file_elements);
	void set_alpha_beta (double Alpha, double Beta);
};
