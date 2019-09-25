#pragma once
#include "Prototype_Mesh.h"
#include "Triangle_Element.h"
#include "Triangulation.h"
#include <set>
#include <map>
#include <algorithm>

namespace NS_Triangular_Mesh_Delaunay
{
#define PRECIS 1e-10
#define MAX_LISTING_SIZE 16
#define MIN_DISTANCE 1e-10
#define AREA_DEGENERATE_TRIANGLE 1e-10

	struct Triangle_data
	{
		int nodes[3];
		int neighbours[3];

		void get_edge (int k, int * n);
		bool has_edge (int n0, int n1);
		int get_neighbour (int * n); // n - nodes on the edge
		int get_neighbour (int k); // k - number of the neighbour
		void replace_neighbour (int * n, int new_neighbour); // n - nodes on the edge, can also be used as set
		void set_nodes (int * n);

		Triangle_data ();
	};

	// rectangle struct
	struct Rectangle {
		double c0[2], cN[2];
		bool point_inside (double * point);
	};

	// tree node struct
	struct Tree_Node {
		Rectangle rect;
		std::vector<int> triangles;
		Tree_Node * next[4];
	};

	struct BF_triangles {
		std::vector <int> triangles;
	};

	class TMD_pointer
	{
	protected:
	public:
		virtual bool triangle_point_inside (int k_triangle, double * point) { return false; };
		virtual void get_triangle_coordinates (int k_triangle, double (*c)[3][2]) { };
	};

	class Tree_of_triangles
	{
	private:
		Tree_Node * root;
		TMD_pointer * tmd;

		// section-section intersection
		bool SS_intersect (double * p1, double * p2, double * p3, double * p4);
		// triangle-rectangle intersection
		bool TR_intersect (const Rectangle & rect, int k_triangle);
		void search_to_add (Tree_Node * tree_node, int k_triangle, std::vector <Tree_Node *> * tree_nodes);
		void search (Tree_Node * tree_node, int k_triangle, std::vector <Tree_Node *> * tree_nodes);
		int search (Tree_Node * tree_node, double * point); // main search function
		void delete_child_nodes (Tree_Node * tree_node);
	public:
		void initialize (double C0[2], double CN[2], TMD_pointer * TMD); // borders of the mesh

		void add_triangle (int k_triangle); // adds triangle
		void reset_triangle (int k_triangle); // removes triangle from the tree and searches for new place for it
		int tree_search (double * point); // main search function
		void node_division (Tree_Node * tree_node);
		void output ();

		Tree_of_triangles ();
		~Tree_of_triangles ();
	};

	class Triangular_Mesh_Delaunay : public Mesh <Node_2D, Element>, public TMD_pointer
	{
	private:
		// point generator
		Triangulation::Point_Generator * pg;
		// tree for access
		Tree_of_triangles * tot;
		// triangle data
		std::vector <Triangle_data> ib_triangles;
		std::vector <BF_triangles> bf_triangles; // basic figures/ triangles relations

		void renumerate () override;
		bool need_to_numerate_edges () override; // returns true for 2D tasks
		void make_elements ();
		int find_area (int k_triangle);
		void split_triangle (int k_node, int k_triangle, int bf);
		bool point_sphere_relation (int k_node, int k_triangle);
		int find_neighbour (int k_triangle, int n0, int n1);
		int find_neighbour (int tr_start, int tr_end, int k_triangle, int n0, int n1);
		int find_neighbouring_node (int * n0, int * n1, int excluding_node);

		void test ();
		bool bad_angles (int k_triangle_1, int k_triangle_2);
		bool neighbours (int k_triangle_1, int k_triangle_2);
		void tr_config (int k_triangle_1, int k_triangle_2, int * n_edge_nodes, int * unique_nodes);

		bool closeness (double * point);
		int point_in_within_bf (int bf, double * point);
	public:
		Triangular_Mesh_Delaunay ();
		Triangular_Mesh_Delaunay (const Triangular_Mesh_Delaunay & tmd);
		~Triangular_Mesh_Delaunay ();

		// tree section
		void get_triangle_coordinates (int k_triangle, double (*c)[3][2]) override;
		bool triangle_point_inside (int k_triangle, double * point) override;
		bool triangle_degenerate (int * n);

		// build mesh
		int build_mesh (char * file_input, char * file_nodes_output, char * file_elements_output);
		int build_mesh_npg (char * file_input, char * file_nodes_output, char * file_elements_output);
		// tree search
		virtual int point_inside (double * coordinates) override; // returns element that has point
																  // same as for triangular mesh
		bool get_isoline_section (int k_element, double * q, double value, double * c1, double * c2) override;
		void output (char * file_nodes, char * file_triangles) override;
	};
}