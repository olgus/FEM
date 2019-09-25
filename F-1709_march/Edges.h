#pragma once
#include <stdio.h>

#define MAX_EDGES_2D 32

class Edge_Structure
{
private:
	int size; // amount of nodes
	int n_entries; // amount of edges

	int * ig, *jg; // ig and jg arrays 
public:
	Edge_Structure (); // constructor
	Edge_Structure (const Edge_Structure & es); // constructor-copy
	~Edge_Structure (); // destructor
	void Copy (const Edge_Structure & es);

	void set_size (int Size, int N_entries); // set size and n_entries
	bool set_ig_jg (int * Ig, int * Jg); // set ig and kg from arrays
	 
	int get_edge_number (int n1, int n2); // returns the number of the edge between nodes n1, n2
	void get_edge_nodes (int n_edge, int * n1, int * n2); // returns nodes that define edge n_edge
	int get_n_entries (); // returns n_entries
	int amount_of_edges_before (int node); // returns the amount of edges with minimal node less that node
	void get_edges_by_min_node (int node, int * start, int * end); // returns amount of and vector of edges that have the node as minimal
};