#include "Edges.h"

Edge_Structure::Edge_Structure ()
{
	ig = NULL;
	jg = NULL;

	size = n_entries = 0;
}

Edge_Structure::Edge_Structure (const Edge_Structure & es)
{
	size = es.size;
	n_entries = es.n_entries;

	if (ig != NULL)
		delete[] ig;
	if (jg != NULL)
		delete[] jg;

	ig = new int[size + 1];
	for (int i = 0; i < size + 1; i++)
	{
		ig[i] = es.ig[i];
	}

	jg = new int[n_entries];
	for (int i = 0; i < n_entries; i++)
	{
		jg[i] = es.jg[i];
	}

}

Edge_Structure::~Edge_Structure ()
{
	if (ig != NULL)
		delete[] ig;
	if (jg != NULL)
		delete[] jg;
}

void Edge_Structure::Copy (const Edge_Structure & es)
{
	size = es.size;
	n_entries = es.n_entries;

	if (ig != NULL)
		delete[] ig;
	if (jg != NULL)
		delete[] jg;

	ig = new int[size + 1];
	for (int i = 0; i < size + 1; i++)
	{
		ig[i] = es.ig[i];
	}

	jg = new int[n_entries];
	for (int i = 0; i < n_entries; i++)
	{
		jg[i] = es.jg[i];
	}

}

void Edge_Structure::set_size (int Size, int N_entries)
{
	size = Size;
	n_entries = N_entries;

	if (ig != NULL)
		delete[] ig;
	if (jg != NULL)
		delete[] jg;

	ig = new int[size + 1];
	jg = new int[n_entries];
}

bool Edge_Structure::set_ig_jg (int * Ig, int * Jg)
{
	for (int i = 0; i < size + 1; i++)
		ig[i] = Ig[i];
	for (int i = 0; i < n_entries; i++)
		jg[i] = Jg[i];
	if (n_entries != Ig[size])
		return false;
	return true;
}

int Edge_Structure::get_edge_number (int n1, int n2)
{
	int n;
	// check that n1 and n2 are legit
	if (n1 > n_entries || n2 > n_entries)
		return -1;

	// n1 has to be smaller than n2
	if (n1 > n2)
	{
		n = n2;
		n2 = n1;
		n1 = n;
	}

	// check in ig section for n1-node
	for (int i = ig[n1]; i < ig[n1 + 1]; i++)
	{
		// if entry for n2 exists, it's the edge
		if (jg[i] == n2)
			return i;
	}

	return -1;
}

void Edge_Structure::get_edge_nodes (int n_edge, int * n1, int * n2)
{
	*n2 = -1;
	*n1 = -1;
	if (n_edge < n_entries)
	{
		*n2 = jg[n_edge];

		int i;
		for (i = 0; i < n_entries && ig[i + 1] <= n_edge; i++);
		*n1 = i;
	}
}

int Edge_Structure::get_n_entries ()
{
	return n_entries;
}

int Edge_Structure::amount_of_edges_before (int node)
{
	return ig[node];
}

void Edge_Structure::get_edges_by_min_node (int node, int * start, int * end)
{
	*start = ig[node];
	*end = ig[node + 1];
}

