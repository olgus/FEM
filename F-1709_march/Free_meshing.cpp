#include "Free_meshing.h"

void Free_mesh::input (char * file_mesh, char * file_density, char * file_order)
{
	FILE * file = fopen (file_mesh, "r");

	// get lvl
	fscanf (file, "%lf", &lvl);

	// get amount of points that define the mesh
	int n;
	fscanf (file, "%i", &n);

	double c0[2];
	double cN[2];
	c0[0] = c0[1] = 1e+15;
	cN[0] = cN[1] = -1e+15;

	double coordinates[2];
	// push nodes with those points
	Node_2D node;
	for (int i = 0; i < n; i++)
	{
		fscanf (file, "%lf %lf", &coordinates[0], &coordinates[1]);
		for (int k = 0; k < 2; k++)
		{
			if (coordinates[k] < c0[k])
				c0[k] = coordinates[k];
			if (coordinates[k] > cN[k])
				cN[k] = coordinates[k];
		}
		node.set_node (coordinates);
		nodes.push_back (node);
	}

	double sc0[2], scN[2];
	for (int i = 0; i < 2; i++)
	{
		sc0[i] = c0[i] - 0.01 * (cN[i] - c0[i]);
		scN[i] = cN[i] + 0.01 * (cN[i] - c0[i]);
	}
	spline = new Spline ();
	spline->set_alpha_beta (alpha, beta);
	spline->prepare (sc0, scN, lvl, file_density);
	int param[] = {SOLVER_METHOD_CGM_SYMM, SOLVER_DECOMP_TYPE_D, 5, 4, SOLVER_MKL_NO};
	spline->solve_task (param);

	// read regions
	fscanf (file, "%i", &n);
	// push them in queue
	Region region;
	int node_number;

	std::vector<Region> base_regions;
	for (int i = 0; i < n; i++)
	{
		fscanf (file, "%i", &region.N);
		fscanf (file, "%i", &region.material);
		region.def_nodes.clear ();
		for (int k = 0; k < region.N; k++)
		{
			fscanf (file, "%i", &node_number);
			region.def_nodes.push_back (node_number);
			//region.inside_nodes_B.push_back (-1);
			//region.inside_nodes_E.push_back (-1);
		}
		region.def_nodes.push_back (region.def_nodes[0]);
		base_regions.push_back (region);
	}

	// set side nodes for base regions
	int found;
	int start, end;
	int J, M;
	for (size_t i = 0, i_end = base_regions.size (); i < i_end; i++)
	{
		// go through each side of the region
		for (int k = 0; k < base_regions[i].N; k++)
		{
			found = 0;
			start = base_regions[i].def_nodes[k];
			end = base_regions[i].def_nodes[k + 1];
			// if side with the same nodes exists in previous region, just copy numbers of nodes
			for (size_t j = 0; j < i && !found; j++)
			{
				for (int m = 0; m < base_regions[j].N && !found; m++)
				{
					if (start == base_regions[j].def_nodes[m] && end == base_regions[j].def_nodes[m + 1]) 
					{
						found = 1;
						J = (int)j;
						M = m;
					}
					if (start == base_regions[j].def_nodes[m + 1] && end == base_regions[j].def_nodes[m])
					{
						found = 2;
						J = (int)j;
						M = m;
					}
				}
			}
			if (found)
			{
				if (found == 1)
				{
					base_regions[i].inside_nodes_B.push_back (base_regions[J].inside_nodes_B[M]);
					base_regions[i].inside_nodes_E.push_back (base_regions[J].inside_nodes_E[M]);
				}
				else
				{
					base_regions[i].inside_nodes_B.push_back (base_regions[J].inside_nodes_E[M]);
					base_regions[i].inside_nodes_E.push_back (base_regions[J].inside_nodes_B[M]);
				}
			}
			else
			{
				nodes[start].get_coordinates (c0);
				nodes[end].get_coordinates (cN);
				make_inside_nodes (c0, cN, &start, &end);
				base_regions[i].inside_nodes_B.push_back (start);
				base_regions[i].inside_nodes_E.push_back (end);
			}
		}
	}

	// make a vector of boundary nodes
	for (size_t i = 0, i_end = base_regions.size (); i < i_end; i++)
	{
		for (int k = 0; k < base_regions[i].N; k++)
		{
			boundary.push_back (base_regions[i].def_nodes[k]);
			if (base_regions[i].inside_nodes_B[k] != -1)
			{
				// cycle between BE
				int CB = base_regions[i].inside_nodes_B[k];
				int CE = base_regions[i].inside_nodes_E[k];
				int CV = CB > CE ? -1 : 1;
				for (int j = CB; CB <= CE ? j <= CE : j >= CE; j += CV)
				{
					boundary.push_back (j);
				}
			}
		}
	}

	for (size_t i = 0, i_end = base_regions.size (); i < i_end; i++)
	{
		regions.push (base_regions[i]);
	}
	fclose (file);

	file = fopen (file_order, "r");
	fscanf (file, "%i", &N_or);
	or_c0 = new double *[N_or];
	or_cN = new double *[N_or];
	for (int i = 0; i < N_or; i++)
	{
		or_c0[i] = new double [2];
		or_cN[i] = new double [2];
		for (int j = 0; j < 2; j++)
		{
			fscanf (file, "%lf", &or_c0[i][j]);
			fscanf (file, "%lf", &or_cN[i][j]);
		}
	}

	for (int i = 0; i < N_or; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			printf ("%lf ", or_c0[i][j]);
			printf ("%lf ", or_cN[i][j]);
		}
		printf ("\n");
	}
	fclose (file);

}

void Free_mesh::output (char * file_nodes, char * file_elements)
{
	double coordinates[2];
	FILE * file = fopen (file_nodes, "w");
	for (size_t i = 0, i_end = nodes.size (); i < i_end; i++)
	{
		nodes[i].get_coordinates (coordinates);
		fprintf (file, "%.16lf %.16lf\n", coordinates[0], coordinates[1]);
	}
	fclose (file);

	int node[3];
	file = fopen (file_elements, "w");
	for (size_t i = 0, i_end = elements.size (); i < i_end; i++)
	{
		elements[i].get_base_nodes (node);
		for (int k = 0; k < 3; k++)
		{
			fprintf (file, "%i ", node[k]);
		}
		fprintf (file, "%i ", elements[i].get_area ());
		int order = 1;
		if (increase_order ((int)i))
			order++;
		fprintf (file, "%i\n", order);
	}
	fclose (file);
}

void Free_mesh::build ()
{
	// go through stack of regions
	Region region; 
	Triangle triangle;
	int tr_nodes[3];
	while (!regions.empty ())
	{
		region = regions.top ();
		regions.pop ();
		printf ("%i", (int)regions.size ());
		int k = 0;
		for (int i = 0; i < region.N; i++)
		{
			if (region.inside_nodes_E[i] == -1)
				k++;
		}
		// if amount of sections == 3 and no side nodes are present
		if (region.N == 3 && k == region.N)
		{
			// then it's a triangle, push it into elements
			for (int i = 0; i < region.N; i++)
			{
				tr_nodes[i] = region.def_nodes[i];
			}
			std::sort (tr_nodes, tr_nodes + 3);
			triangle.set_area (region.material);
			triangle.set_element (tr_nodes);
			triangle.set_base_nodes (tr_nodes);
			elements.push_back (triangle);
		}
		else
		{
			// otherwise
			// build vector of node numbers
			std::vector <int> region_nodes;
			for (int i = 0; i < region.N; i++)
			{
				region_nodes.push_back (region.def_nodes[i]);
				if (region.inside_nodes_B[i] != -1)
				{
					// cycle between BE
					int CB = region.inside_nodes_B[i];
					int CE = region.inside_nodes_E[i];
					int CV = CB > CE ? -1 : 1;
					for (int j = CB; CB <= CE ? j <= CE : j >= CE; j += CV)
					{
						region_nodes.push_back (j);
					}
				}
			}
			region_nodes.push_back (region_nodes[0]);
			region_nodes.push_back (region_nodes[1]);

			// look for the line that separates the region into two with 4 angles closest to pi / 3
			// cycle through nodes
			double c[2], c0[2], c1[2];
			double discr;
			double min_discr = 1e+15;
			int I, J;
			double a[4];
			for (size_t i = 1, i_end = region_nodes.size () - 1; i < i_end; i++)
			{
				// cycle through nodes after i-node
				// start with node after current side
				for (size_t j = i + 1; j < i_end; j++)
				{
					// calculate all 4 angles
					discr = 0.0;
					// 1 angle between i and previous to it
					nodes[region_nodes[i]].get_coordinates (c);
					nodes[region_nodes[i - 1]].get_coordinates (c0);
					nodes[region_nodes[j]].get_coordinates (c1);
					a[0] = calc_angle (c, c0, c1);
					// 2 angle between i and next to it
					nodes[region_nodes[i + 1]].get_coordinates (c0);
					a[1] = calc_angle (c, c0, c1);
					// 3 angle between j and previous to it
					nodes[region_nodes[j]].get_coordinates (c);
					nodes[region_nodes[j - 1]].get_coordinates (c0);
					nodes[region_nodes[i]].get_coordinates (c1);
					a[2] = calc_angle (c, c0, c1);
					// 4 angle between j and next to it
					nodes[region_nodes[j + 1]].get_coordinates (c0);
					a[3] = calc_angle (c, c0, c1);
					bool line = true;
					for (int k = 0; k < 4 && line; k++)
					{
						if (fabs (a[k] - PI_MESH) > 1e-5 && fabs (a[k]) > 1e-5)
							discr += fabs (a[k] - SIXTYDEGREESRADIAN);
						else
							line = false;
					}
					if (line)
					{
						// if all of them are closer to 60 then previous ones, save them and the line
						if (discr < min_discr)
						{
							min_discr = discr;
							I = (int)i;
							J = (int)j;
						}
					}
				}
			}

			// line found
			// build nodes on the line
			int n1, n2;
			nodes[region_nodes[I]].get_coordinates (c0);
			nodes[region_nodes[J]].get_coordinates (c1);
			make_inside_nodes (c0, c1, &n1, &n2);

			// build base nodes
			std::vector <int> base_nodes;
			std::vector <int> new_B;
			std::vector <int> new_E;
			int R0, R1;
			// go through edges
			int p;
			for (int i = 0; i < region.N; i++)
			{
				base_nodes.push_back (region.def_nodes[i]);
				p = -1;
				// check if I or J is a base node
				// if it, save them and push B and E later
				if (region_nodes[I] == region.def_nodes[i])
				{
					R0 = (int)base_nodes.size () - 1;
				}
				if (region_nodes[J] == region.def_nodes[i])
				{
					R1 = (int)base_nodes.size () - 1;
				}
				// if I or J on the edge
				if ((region.inside_nodes_B[i] <= region_nodes[I] && region_nodes[I] <= region.inside_nodes_E[i])
					|| (region.inside_nodes_E[i] <= region_nodes[I] && region_nodes[I] <= region.inside_nodes_B[i]))
					p = I;
				if ((region.inside_nodes_B[i] <= region_nodes[J] && region_nodes[J] <= region.inside_nodes_E[i])
					|| (region.inside_nodes_E[i] <= region_nodes[J] && region_nodes[J] <= region.inside_nodes_B[i]))
					p = J;

				if (p != -1)
				{
					// push new base node
					base_nodes.push_back (region_nodes[p]);
					// check if it is inside
					if ((region.inside_nodes_B[i] < region_nodes[p] && region_nodes[p] < region.inside_nodes_E[i])
						|| (region.inside_nodes_E[i] < region_nodes[p] && region_nodes[p] < region.inside_nodes_B[i]))
					{
						// push B
						new_B.push_back (region.inside_nodes_B[i]);
						// push new E
						new_E.push_back (region_nodes[p - 1]);
						// push new B
						new_B.push_back (region_nodes[p + 1]);
						// push E
						new_E.push_back (region.inside_nodes_E[i]);
					}
					else
					{
						// if there is only node
						if (region.inside_nodes_B[i] == region.inside_nodes_E[i])
						{
							new_B.push_back (-1);
							new_E.push_back (-1);
							new_B.push_back (-1);
							new_E.push_back (-1);
						}
						else
						{
							// or on one of the boundaries
							if (region.inside_nodes_B[i] == region_nodes[p])
							{
								// if it's left boundary
								new_B.push_back (-1);
								new_E.push_back (-1);
								// push new B
								new_B.push_back (region_nodes[p + 1]);
								// push E
								new_E.push_back (region.inside_nodes_E[i]);
							}
							else
							{
								// push B
								new_B.push_back (region.inside_nodes_B[i]);
								// push new E
								new_E.push_back (region_nodes[p - 1]);
								// if it's right boundary
								new_B.push_back (-1);
								new_E.push_back (-1);
							}
						}
					}
					// save I or J positions
					if (p == I)
						R0 = (int)base_nodes.size () - 1;
					if (p == J)
						R1 = (int)base_nodes.size () - 1;
				}
				else
				{
					// otherwise push B and E
					new_B.push_back (region.inside_nodes_B[i]);
					new_E.push_back (region.inside_nodes_E[i]);
				}
			}

			// build two new regions
			Region region1, region2;
			region1.material = region.material;
			region2.material = region.material;

			if (R0 > R1)
			{
				int R = R1;
				R1 = R0;
				R0 = R;
				R = n1;
				n1 = n2;
				n2 = R;
			}
			// push sides before R0 into region1
			for (int i = 0; i < R0; i++)
			{
				region1.def_nodes.push_back (base_nodes[i]);
				region1.inside_nodes_B.push_back (new_B[i]);
				region1.inside_nodes_E.push_back (new_E[i]);
			}
			// push side between R0 and R1 into region1
			region1.def_nodes.push_back (base_nodes[R0]);
			region1.inside_nodes_B.push_back (n1);
			region1.inside_nodes_E.push_back (n2);
			// push side between R1 and R0 into region2
			region2.def_nodes.push_back (base_nodes[R1]);
			region2.inside_nodes_B.push_back (n2);
			region2.inside_nodes_E.push_back (n1);
			// push sides between R0 and R1 into region2
			for (int i = R0; i < R1; i++)
			{
				region2.def_nodes.push_back (base_nodes[i]);
				region2.inside_nodes_B.push_back (new_B[i]);
				region2.inside_nodes_E.push_back (new_E[i]);
			}
			// push sides after R1 into region1
			for (int i = R1, i_end = (int)base_nodes.size (); i < i_end; i++)
			{
				region1.def_nodes.push_back (base_nodes[i]);
				region1.inside_nodes_B.push_back (new_B[i]);
				region1.inside_nodes_E.push_back (new_E[i]);
			}
			region1.def_nodes.push_back (region1.def_nodes[0]);
			region2.def_nodes.push_back (region2.def_nodes[0]);
			region1.N = (int)region1.def_nodes.size () - 1;
			region2.N = (int)region2.def_nodes.size () - 1;
			regions.push (region1);
			regions.push (region2);
		}
	}
}

void Free_mesh::relax ()
{
	std::vector<std::vector<int>> star_nodes;
	star (&star_nodes);
	for (int k = 0; k < MAX_MESH_RELAXATION_ITER; k++)
	{
		for (size_t i = 0, i_end = nodes.size (); i < i_end; i++)
		{
			// if the nodes is not on the boundary
			if (std::find (boundary.begin (), boundary.end (), i) == boundary.end ())
			{
				double h;
				double c[2];
				nodes[i].get_coordinates (c);
				h = density_function_value (c);
				// go by star nodes
				double hn, hm, l;
				double cn[2];
				double sx, sy, s;
				sx = sy = s = 0;
				for (size_t j = 0, j_end = star_nodes[i].size (); j < j_end; j++)
				{
					nodes[star_nodes[i][j]].get_coordinates (cn);
					hn = density_function_value (cn);
					hm = (h + hn) / 2.0;
					l = distance ((int)i, star_nodes[i][j]);
					
					sx += cn[0] * l / hm;
					sy += cn[1] * l / hm;
					s += l / hm;
				}

				c[0] = sx / s;
				c[1] = sy / s;
				nodes[i].set_coordinates (c);
			}
		}
	}
}

void Free_mesh::renumerate ()
{
	std::queue <int> Q;
	std::vector < std::vector <int> > node_links ((int)nodes.size(), std::vector <int> (0));
	int d_nodes;
	int * n;

	for (std::size_t i = 0, vsize = elements.size (); i < vsize; i++)
	{
		d_nodes = elements[i].get_amount_of_def_nodes ();
		n = new int[d_nodes];
		for (int j = 0; j < d_nodes; j++)
		{
			for (int k = 0; k < d_nodes; k++)
			{
				if (k != j)
				{
					elements[i].get_def_nodes (n);
					if (std::find (node_links[n[j]].begin (), node_links[n[j]].end (), n[k]) == node_links[n[j]].end ())
						node_links[n[j]].push_back (n[k]);
				}
			}
		}
		delete[] n;
	}

	for (std::size_t i = 0, vsize = node_links.size (); i < vsize; i++)
	{
		std::sort (node_links[i].begin (), node_links[i].end ());
	}

	std::vector <Node_2D> NF;
	int node_max = find_maximal_node ();
	NF.push_back (nodes[node_max]);

	for (std::size_t i = 0, vsize = node_links[node_max].size (); i < vsize; i++)
	{
		Q.push (node_links[node_max][i]);
	}

	int node_current;
	while (!Q.empty ())
	{
		node_current = Q.front ();
		Q.pop ();
		if (std::find (NF.begin (), NF.end (), nodes[node_current]) == NF.end ()) // if node hasn't been added to NF
		{
			NF.push_back (nodes[node_current]);
			for (unsigned int i = 0; i < node_links[node_current].size (); i++)
			{
				Q.push (node_links[node_current][i]);
			}
		}
	}
	std::reverse (NF.begin (), NF.end ());

	bool found;
	int * node_defining;
	int n_defining;
	for (unsigned int i = 0; i < elements.size (); i++)
	{
		n_defining = elements[i].get_amount_of_def_nodes ();
		node_defining = new int[n_defining];
		elements[i].get_def_nodes (node_defining);

		for (int j = 0; j < n_defining; j++)
		{
			// change j-node from nodes to new node from NF
			found = false;
			for (std::size_t k = 0, vsize = NF.size (); k < vsize && !found; k++)
			{
				if (NF[k] == nodes[node_defining[j]])
				{
					found = true;
					node_defining[j] = static_cast<int> (k);
				}
			}
			elements[i].set_nodes (node_defining);
		}
		delete[] node_defining;

		n_defining = elements[i].get_amount_of_base_nodes ();
		node_defining = new int[n_defining];
		elements[i].get_base_nodes (node_defining);
		for (int j = 0; j < n_defining; j++)
		{
			// change j-node from base nodes to new node from NF
			found = false;
			for (std::size_t k = 0, vsize = NF.size (); k < vsize && !found; k++)
			{
				if (NF[k] == nodes[node_defining[j]])
				{
					found = true;
					node_defining[j] = static_cast<int> (k);
				}
			}
			elements[i].set_base_nodes (node_defining);
		}

		elements[i].sort_def_nodes ();
		delete[] node_defining;
	}

	nodes.clear ();
	nodes.insert (nodes.begin (), NF.begin (), NF.end ());
}

void Free_mesh::make_inside_nodes (double * c0, double * cN, int * n1, int * n2)
{
	std::vector <Node_2D> local_nodes;
	Node_2D node;
	double coordinates[2];
	double h, h1;

	*n1 = (int)nodes.size ();
	h = density_function_value (c0);
	h1 = density_function_value (cN);
	
	// start with the smallest step
	bool switched = false;
	if (h1 < h)
	{
		for (int i = 0; i < 2; i++)
		{
			h = c0[i];
			c0[i] = cN[i];
			cN[i] = h;
		}
		switched = true;
	}

	for (int i = 0; i < 2; i++)
	{
		coordinates[i] = c0[i];
	}

	// length of the section 
	double L = sqrt (pow (cN[0] - c0[0], 2.0) + pow (cN[1] - c0[1], 2.0));

	h = density_function_value (coordinates);

	// cycle till last bit
	while (L > h)
	{
		h = density_function_value (coordinates);
		// add a node with h distance from current
		for (int i = 0; i < 2; i++)
		{
			coordinates[i] = (cN[i] - coordinates[i]) * h / L + coordinates[i];
		}
		node.set_coordinates (coordinates);
		local_nodes.push_back (node);
		L -= h;
	}
	// if there weren't any new nodes, do nothing at all
	if ((int)local_nodes.size () > 0)
	{
		// check how L and h relate
		if (L / h > 0.8 /*&& L / h < 1.0*/) 
		{
			// put it in the middle 
			if ((int)local_nodes.size () == 1)
			{
				coordinates[0] = c0[0];
				coordinates[1] = c0[1];
			}
			else
			{
				local_nodes[(int)local_nodes.size () - 2].get_coordinates (coordinates);

			}

			//L = sqrt (pow (cN[0] - coordinates[0], 2.0) + pow (cN[1] - coordinates[1], 2.0));
			//coordinates[0] = (cN[0] - coordinates[0]) / 2.0 + coordinates[0];
			//coordinates[1] = (cN[1] - coordinates[1]) / 2.0 + coordinates[1];
			coordinates[0] = (coordinates[0] + cN[0]) / 2.0;
			coordinates[1] = (coordinates[1] + cN[1]) / 2.0;

			local_nodes[(int)local_nodes.size () - 1].set_coordinates (coordinates);
		}
		else
		{
			//if (L / h < 0.5)
				local_nodes.pop_back ();
		}
	}

	// save these nodes
	nodes.insert (nodes.end (), local_nodes.begin (), local_nodes.end ());

	if ((int)nodes.size () == *n1)
	{
		*n1 = -1;
		*n2 = -1;
	}
	else
	{
		*n2 = (int)nodes.size () - 1;
	}
	if (switched)
	{
		int n = *n1;
		*n1 = *n2;
		*n2 = n;
	}
}

double Free_mesh::calc_angle (double * c, double * c0, double * c1)
{
	double a[2], b[2];
	for (int i = 0; i < 2; i++)
	{
		a[i] = c0[i] - c[i];
		b[i] = c1[i] - c[i];
	}
	double aL = sqrt (pow (a[0], 2.0) + pow (a[1], 2.0));
	double bL = sqrt (pow (b[0], 2.0) + pow (b[1], 2.0));
	double sp = a[0] * b[0] + a[1] * b[1];
	double ca = sp / (aL * bL);
	return acos (ca);
}

void Free_mesh::star (std::vector<std::vector<int>> * star_nodes)
{
	for (size_t i = 0, i_end = nodes.size (); i < i_end; i++)
	{
		std::vector <int> v;
		star_nodes->push_back (v);
	}

	// go by elements
	int base_nodes[3];
	for (size_t i = 0, i_end = elements.size (); i < i_end; i++)
	{
		// for each base_node
		elements[i].get_base_nodes (base_nodes);
		for (int k = 0; k < 3; k++)
		{
			for (int m = 0; m < 3; m++)
			{
				// add other two nodes intro star_nodes, if they weren't there before
				if (k != m)
				{
					if (std::find ((*star_nodes)[base_nodes[k]].begin (), (*star_nodes)[base_nodes[k]].end (), base_nodes[m]) == (*star_nodes)[base_nodes[k]].end())
					{
						(*star_nodes)[base_nodes[k]].push_back (base_nodes[m]);
					}
				}
			}
		}
	}
}

double Free_mesh::distance (int n1, int n2)
{
	double c1[2], c2[2];
	nodes[n1].get_coordinates (c1);
	nodes[n2].get_coordinates (c2);

	return sqrt (pow (c1[0] - c2[0], 2.0) + pow (c1[1] - c2[1], 2.0));
}

double Free_mesh::density_function_value (double * c)
{
	double h;
	spline->get_solution_in_point (0, c, &h);
	return h;
}

int Free_mesh::find_maximal_node ()
{
	double l, r;
	l = -1e+15;
	int res = 0;

	for (unsigned int i = 0; i < nodes.size (); i++)
	{
		r = nodes[i].distance_origin ();
		if (r > l)
		{
			res = (int)i;
			l = r;
		}
	}
	return res;
}

bool Free_mesh::point_in_order_region (int k_region, double * coordinates)
{
	if ((or_c0[k_region][0] < coordinates[0] && coordinates[0] < or_cN[k_region][0]) &&
		(or_c0[k_region][1] < coordinates[1] && coordinates[1] < or_cN[k_region][1]))
		return true;
	return false;
}

bool Free_mesh::increase_order (int k_element)
{
	double coordinates[2];
	double cn[2];
	for (int k = 0; k < 2; k++)
	{
		coordinates[k] = 0.0;
	}
	int base_nodes[3];

	// get mass center
	elements[k_element].get_base_nodes (base_nodes);
	for (int j = 0; j < 3; j++)
	{
		nodes[base_nodes[j]].get_coordinates (cn);
		for (int k = 0; k < 2; k++)
		{
			coordinates[k] += cn[k];
		}
	}
	for (int k = 0; k < 2; k++)
	{
		coordinates[k] /= 3.0;
	}

	// go through regions 
	for (int i = 0; i < N_or; i++)
	{	
		// check if point is in the region
		if (point_in_order_region (i, coordinates))
			return true;
	}
	return false;
}

Free_mesh::Free_mesh ()
{
	or_c0 = NULL;
	or_cN = NULL;
	alpha = 1e-5;
	beta = 1e-10;
}

Free_mesh::Free_mesh (const Free_mesh & free_mesh)
{
}

Free_mesh::~Free_mesh ()
{
	if (or_c0 != NULL)
	{
		for (int i = 0; i < N_or; i++)
		{
			delete[] or_c0[i];
		}
		delete[] or_c0;
	}
	if (or_cN != NULL)
	{
		for (int i = 0; i < N_or; i++)
		{
			delete[] or_cN[i];
		}
		delete[] or_cN;
	}
}

void Free_mesh::build_Mesh (char * file_mesh, char * file_density, char * file_order, char * file_nodes, char * file_elements)
{
	input (file_mesh, file_density, file_order);
	printf ("\tbuilding the mesh\n");
	build ();
	printf ("\n\trelaxing the mesh\n");
	relax ();
	printf ("\trenumerating the mesh\n");
	renumerate ();
	output (file_nodes, file_elements);
	printf ("\tfree-meshing section finished building\n");
}

void Free_mesh::build_Mesh_wo_renumerate (char * file_mesh, char * file_density, char * file_order, char * file_nodes, char * file_elements)
{
	input (file_mesh, file_density, file_order);
	printf ("\tbuilding the mesh\n");
	build ();
	printf ("\n\trelaxing the mesh\n");
	relax ();
	printf ("\trenumerating the mesh\n");
	output (file_nodes, file_elements);
	printf ("\tfree-meshing section finished building\n");
}

void Free_mesh::set_alpha_beta (double Alpha, double Beta)
{
	alpha = Alpha;
	beta = Beta;
}
