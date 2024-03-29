#include "Mixed_Triangular_Mesh.h"

void Mixed_Triangular_Mesh::input_mesh_data (char * file_name)
{
	FILE * file = fopen (file_name, "r");

	fscanf (file, "%i %i", &n_axis[0], &n_axis[1]);
	fscanf (file, "%i", &lvl);

	N = (n_axis[0]) * (n_axis[1] + 1);
	M = (n_axis[0] + 1) * (n_axis[1]);

	M > N ? max = M : max = N;

	sides = new mesh_side_MM *[2]; // allocate the memory for mesh tetragonal areas
	for (int i = 0; i < 2; i++)
		sides[i] = new mesh_side_MM[max];

	tetra_nodes = new Point_2D *[n_axis[0] + 1]; // allocate the memory for tetragonal areas' nodes
	for (int i = 0; i < n_axis[0] + 1; i++)
		tetra_nodes[i] = new Point_2D[n_axis[1] + 1];

	// get nodes
	double c0, c1;

	for (int j = 0; j < n_axis[1] + 1; j++)
	{
		for (int i = 0; i < n_axis[0] + 1; i++)
		{
			fscanf (file, "%lf %lf", &c0, &c1);
			tetra_nodes[i][j].set_dim (2);
			tetra_nodes[i][j].set_point (c0, c1);
		}
	}
	// set boundary values
	coord0 = new double[2];
	coordN = new double[2];
	coord0[0] = tetra_nodes[0][0].X ();
	coord0[1] = tetra_nodes[0][0].Y ();
	coordN[0] = tetra_nodes[n_axis[0]][n_axis[1]].X ();
	coordN[1] = tetra_nodes[n_axis[0]][n_axis[1]].Y ();

	n_areas = n_axis[0] * n_axis[1];
	order = new int[n_areas];
	// get order
	for (int i = 0; i < n_areas; i++)
	{
		fscanf (file, "%i", &order[i]);
	}
	
	materials = new int[n_areas];
	// get materials
	for (int i = 0; i < n_areas; i++)
	{
		fscanf (file, "%i", &materials[i]);
	}

	// get data about vertical sides
	// amounts of nodes
	for (int i = 0; i < M; i++)
	{
		fscanf (file, "%i", &(sides[0][i].sections));
		sides[0][i].sections *= (int)round (pow (2, lvl));
	}
	// coef
	for (int i = 0; i < M; i++)
	{
		fscanf (file, "%lf", &(sides[0][i].coef));
	}
	// direction
	for (int i = 0; i < M; i++)
	{
		fscanf (file, "%i", &(sides[0][i].direc));
	}

	// get data about horizontal sides
	// amounts of nodes
	for (int i = 0; i < N; i++)
	{
		fscanf (file, "%i", &(sides[1][i].sections));
		sides[1][i].sections *= (int)round (pow (2, lvl));
	}
	// coef
	for (int i = 0; i < N; i++)
	{
		fscanf (file, "%lf", &(sides[1][i].coef));
	}
	// direction
	for (int i = 0; i < N; i++)
	{
		fscanf (file, "%i", &(sides[1][i].direc));
	}

	// change sides' coefs if direction is down the side
	for (int i = 0; i < N; i++)
	{
		if (sides[1][i].direc == -1)
			sides[1][i].coef = 1.0 / sides[1][i].coef;
	}
	for (int i = 0; i < M; i++)
	{
		if (sides[0][i].direc == -1)
			sides[0][i].coef = 1.0 / sides[0][i].coef;
	}

	fclose (file);
}

void Mixed_Triangular_Mesh::build_mesh (double * C0, double * CN, int * N_axis)
{
	n_axis[0] = 1;
	n_axis[1] = 1;

	M > N ? max = M : max = N;

	sides = new mesh_side_MM *[2]; // allocate the memory for mesh tetragonal areas
	for (int i = 0; i < 2; i++)
		sides[i] = new mesh_side_MM[2];

	tetra_nodes = new Point_2D *[n_axis[0] + 1]; // allocate the memory for tetragonal areas' nodes
	for (int i = 0; i < n_axis[0] + 1; i++)
		tetra_nodes[i] = new Point_2D[n_axis[1] + 1];

	for (int i = 0; i < 2; i++)
	{
		coord0[i] = C0[i];
		coordN[i] = CN[i];
	}

	// set nodes of the area
	tetra_nodes[0][0].set_dim (2);
	tetra_nodes[0][0].set_point (coord0[0], coord0[1]);
	tetra_nodes[0][1].set_dim (2);
	tetra_nodes[0][1].set_point (coordN[0], coord0[1]);
	tetra_nodes[1][0].set_dim (2);
	tetra_nodes[1][0].set_point (coord0[0], coordN[1]);
	tetra_nodes[1][1].set_dim (2);
	tetra_nodes[1][1].set_point (coordN[0], coordN[1]);

	// set sides' data
	for (int i = 0; i < 2; i++)
	{
		sides[0][i].sections = N_axis[1];
		sides[0][i].coef = 1.0;
		sides[0][i].direc = 1;

		sides[1][i].sections = N_axis[0];
		sides[1][i].coef = 1.0;
		sides[1][i].direc = 1;
	}

	// materials
	materials = new int[1];
	materials[0] = 1;

	// make that mesh
	make_init_Mesh ();
}

bool Mixed_Triangular_Mesh::build_Mesh (char * file_name_nodes, char * file_name_elements, int * N_functions)
{
	for (int i = 0; i < dim; i++)
	{
		coord0[i] = 1e+15;
		coordN[i] = -1e+15;
	}

	// get nodes
	FILE * file = fopen (file_name_nodes, "r");
	double * coordinates = new double[dim];
	Node_2D node;
	while (!feof (file))
	{
		for (int i = 0; i < dim; i++)
		{
			fscanf (file, "%lf", &coordinates[i]);
			if (coordinates[i] > coordN[i])
			{
				coordN[i] = coordinates[i];
			}
			if (coordinates[i] < coord0[i])
			{
				coord0[i] = coordinates[i];
			}
		}
		node.set_coordinates (coordinates);
		nodes.push_back (node);
	}
	nodes.pop_back ();
	n_nodes = (int)nodes.size ();
	fclose (file);
	printf ("\tgot %i nodes\n", n_nodes);

	// get elements
	file = fopen (file_name_elements, "r");

	int order;
	while (!feof (file))
	{
		int n_def_nodes = 3;
		int area;
		int * element_nodes = new int[n_def_nodes];
		for (int i = 0; i < n_def_nodes; i++)
		{
			fscanf (file, "%i", &element_nodes[i]);
		}
		fscanf (file, "%i", &area);
		// scanfs order which is not needed typically, except mixed_tr_mesh
		fscanf (file, "%i", &order);

		if (order == 1)
		{
			std::unique_ptr <Triangle> element = std::make_unique<Triangle> ();
			element->set_base_nodes (element_nodes);
			element->set_element (element_nodes);
			element->set_area (area);
			elements.push_back (std::move (element));
		}
		if (order == 2)
		{
			std::unique_ptr <Triangle> element = std::make_unique<Triangle> ();
			element->set_base_nodes (element_nodes);
			element->set_element (element_nodes);
			element->set_area (area);
			elements.push_back (std::move (element));
		}
		delete[] element_nodes;

	}
	elements.pop_back ();
	n_elements = (int)elements.size ();
	fclose (file);
	printf ("\tgot %i elements\n", n_elements);

	if (need_to_numerate_edges ())
		numerate_edges ();
	*N_functions = numerate_functions ();
	fold_mesh ();

	delete[] coordinates;
	return true;
}

void Mixed_Triangular_Mesh::wrap_material (int k_material, double * c0, double * cN)
{	
	int n_base_nodes = elements[0]->get_amount_of_base_nodes ();
	int * base_nodes = new int[n_base_nodes];
	int material;

	for (int i = 0; i < dim; i++)
	{
		c0[i] = 1e+15;
		cN[i] = -1e+15;
	}

	double * coordinates = new double[dim];
	for (int i = 0; i < n_elements; i++)
	{
		material = elements[i]->get_area ();

		if (material == k_material)
		{
			elements[i]->get_base_nodes (base_nodes);
			for (int j = 0; j < n_base_nodes; j++)
			{
				nodes[base_nodes[j]].get_coordinates (coordinates);
				for (int k = 0; k < dim; k++)
				{
					if (coordinates[k] < c0[k])
					{
						c0[k] = coordinates[k];
					}
					if (coordinates[k] > cN[k])
					{
						cN[k] = coordinates[k];
					}
				}
			}
		}
	}
	delete[] coordinates;
	delete[] base_nodes;
}

int Mixed_Triangular_Mesh::point_inside (double * coordinates)
{
	int p_element = -1;
	// look through folds
	for (int i = 0, i_end = (int)folds.size (); i < i_end; i++)
	{
		// if the point is inside a fold (&&)
		if ((folds[i].c0[0] - 1e-10 < coordinates[0] && coordinates[0] < folds[i].cN[0] + 1e-10) &&
			(folds[i].c0[1] - 1e-10 < coordinates[1] && coordinates[1] < folds[i].cN[1] + 1e-10))
		{
			// call elements->points_inside for only elements in the fold 
			for (int k = 0, k_end = (int)folds[i].elements.size (); k < k_end && p_element == -1; k++)
			{
				if (elements[folds[i].elements[k]]->point_inside (*this, coordinates))
				{
					p_element = folds[i].elements[k];
				}
			}
		}
	}
	return p_element;
}

int Mixed_Triangular_Mesh::get_amount_of_folds ()
{
	return (int)folds.size ();
}

void Mixed_Triangular_Mesh::get_fold_coordinates (int k_fold, double * c0, double * cN)
{
	if (k_fold < folds.size ())
	{
		for (int i = 0; i < dim; i++)
		{
			c0[i] = folds[k_fold].c0[i];
			cN[i] = folds[k_fold].cN[i];
		}
	}
}

bool Mixed_Triangular_Mesh::make_init_Mesh ()
{
	double q; // current coef
	int n; // amount of nodes on the side
	double GPS; // geometry progression sum
	double L; // full length of the side
	double l; // current length
	double l0; // start length
	Node_2D node; // new node
	double c[2]; // x, y coordinates of the new node
	double unit[2]; // unit vector on the side
	int sl, sr, sd, su; // numbers of sides for the area
	Point_2D * p1, *p2; // points for the side building
	int area; // area numer through i and j
	int or; // order number

			  // do i need that?
			  //p1 = new Point_2D ();
			  //p2 = new Point_2D ();

	for (int j = 0; j < n_axis[1]; j++) // go through areas
	{
		for (int i = 0; i < n_axis[0]; i++)
		{
			area = materials[j * (n_axis[0]) + i];
			or = order[j * (n_axis[0]) + i];

			// figure out which sides make the area
			sl = j * (n_axis[0] + 1) + i;
			sr = j * (n_axis[0] + 1) + (i + 1);
			sd = j + i * (n_axis[1] + 1);
			su = (j + 1) + i * (n_axis[1] + 1);

			if (sides[1][su].sections == sides[1][sd].sections) // up and down become the basis
			{
				// do up and down side fully
				// direction should technically always be up the axis

				// down
				// get length of the down side
				L = sqrt (pow (tetra_nodes[i][j].X () - tetra_nodes[i + 1][j].X (), 2.0) +
					pow (tetra_nodes[i][j].Y () - tetra_nodes[i + 1][j].Y (), 2.0));

				// get coef for the side
				q = sides[1][sd].coef;
				// get geometric sum
				if (fabs (q - 1.0) < ZERO_DIFFERENCE_Mixed_Mesh) // or if coef == 1, just the amount of sections
					GPS = sides[1][sd].sections;
				else
					GPS = (1.0 - pow (q, sides[1][sd].sections)) / (1.0 - q);

				c[0] = tetra_nodes[i][j].X ();
				c[1] = tetra_nodes[i][j].Y ();
				node.set_coordinates (c);
				if (std::find (nodes.begin (), nodes.end (), node) == nodes.end ())
				{
					nodes.push_back (node);
					n_nodes++;
				}
				if (j == 0) // fill down side if the row is first
				{
					sides[1][sd].nodes.push_back (node);
				}

				// unit vector on the side
				unit[0] = (tetra_nodes[i + 1][j].X () - tetra_nodes[i][j].X ()) / L;
				unit[1] = (tetra_nodes[i + 1][j].Y () - tetra_nodes[i][j].Y ()) / L;

				// get length of the first section, again, always the left one for the down side
				l0 = L / GPS;
				for (int k = 1; k < sides[1][sd].sections + 1; k++) // push all the nodes 
				{
					if (fabs (q - 1.0) < ZERO_DIFFERENCE_Mixed_Mesh)
						l = l0 * k;
					else
						l = l0 * (1.0 - pow (q, k)) / (1.0 - q);

					c[0] = unit[0] * l + tetra_nodes[i][j].X ();
					c[1] = unit[1] * l + tetra_nodes[i][j].Y ();

					node.set_coordinates (c);
					if (std::find (nodes.begin (), nodes.end (), node) == nodes.end ())
					{
						nodes.push_back (node);
						n_nodes++;
					}
					if (j == 0)
						sides[1][sd].nodes.push_back (node);
				}
				// down side should *technically be done
				// do the up one then
				// get length of the up side
				L = sqrt (pow (tetra_nodes[i][j + 1].X () - tetra_nodes[i + 1][j + 1].X (), 2.0) +
					pow (tetra_nodes[i][j + 1].Y () - tetra_nodes[i + 1][j + 1].Y (), 2.0));

				// get coef for the side
				q = sides[1][su].coef;
				// get geometric sum
				if (fabs (q - 1.0) < ZERO_DIFFERENCE_Mixed_Mesh) // or if coef == 1, just the amount of sections
					GPS = sides[1][su].sections;
				else
					GPS = (1.0 - pow (q, sides[1][su].sections)) / (1.0 - q);

				// starting point, always left one for the up side, so has to be pushed only if i = 0
				c[0] = tetra_nodes[i][j + 1].X ();
				c[1] = tetra_nodes[i][j + 1].Y ();
				node.set_coordinates (c);
				sides[1][su].nodes.push_back (node);
				if (std::find (nodes.begin (), nodes.end (), node) == nodes.end ())
				{
					nodes.push_back (node);
					n_nodes++;
				}

				// unit vector on the side
				unit[0] = (tetra_nodes[i + 1][j + 1].X () - tetra_nodes[i][j + 1].X ()) / L;
				unit[1] = (tetra_nodes[i + 1][j + 1].Y () - tetra_nodes[i][j + 1].Y ()) / L;

				// get length of the first section, again, always the left one for the down side
				l0 = L / GPS;
				for (int k = 1; k < sides[1][su].sections + 1; k++) // push all the nodes 
				{
					if (fabs (q - 1.0) < ZERO_DIFFERENCE_Mixed_Mesh) // or if coef == 1, just the amount of sections
						l = l0 * k;
					else
						l = l0 * (1.0 - pow (q, k)) / (1.0 - q);

					c[0] = unit[0] * l + tetra_nodes[i][j + 1].X ();
					c[1] = unit[1] * l + tetra_nodes[i][j + 1].Y ();

					node.set_coordinates (c);
					if (std::find (nodes.begin (), nodes.end (), node) == nodes.end ())
					{
						nodes.push_back (node);
						n_nodes++;
					}
					sides[1][su].nodes.push_back (node);
				}
				// go through nodes on the up and down side
				for (int t = 0; t < sides[1][sd].sections + 1; t++)
				{
					// clear triangle sections
					t0.clear ();
					// figure out starting point and finishing point
					// always go up the side!
					p1 = sides[1][sd].nodes[t].get_Point ();
					p2 = sides[1][su].nodes[t].get_Point ();
					// count q and n for the side through qs and ns of left and right
					n = (int)((1.0 - (double)(t) / (double)(sides[1][sd].sections)) * sides[0][sl].sections +
						((double)t / (double)(sides[1][sd].sections)) * sides[0][sr].sections);
					q = ((1.0 - (double)(t) / (double)(sides[1][sd].sections)) * sides[0][sl].coef +
						((double)t / (double)(sides[1][sd].sections)) * sides[0][sr].coef);

					// get length of the side to slice
					L = sqrt (pow (p1->X () - p2->X (), 2.0) +
						pow (p1->Y () - p2->Y (), 2.0));

					// get geometric sum
					if (fabs (q - 1.0) < ZERO_DIFFERENCE_Mixed_Mesh) // or if coef == 1, just the amount of sections
						GPS = n;
					else
						GPS = (1.0 - pow (q, n)) / (1.0 - q);

					// unit vector on the side
					unit[0] = (p2->X () - p1->X ()) / L;
					unit[1] = (p2->Y () - p1->Y ()) / L;

					// get length of the first section, again, always the down one for the left side
					l0 = L / GPS;

					node.set_Point (*p1);
					t0.push_back (node);

					if (t == 0) // if current side is left, push down left node there
					{
						sides[0][sl].nodes.push_back (node);
					}
					if (t == sides[1][su].sections) // if current side is right, push down right node there
					{
						sides[0][sr].nodes.push_back (node);
					}
					for (int k = 1; k < n; k++) // push nodes between basises
					{
						if (fabs (q - 1.0) < ZERO_DIFFERENCE_Mixed_Mesh) // or if coef == 1, just the amount of sections
							l = l0 * k;
						else
							l = l0 * (1.0 - pow (q, k)) / (1.0 - q);
						c[0] = unit[0] * l + p1->X ();
						c[1] = unit[1] * l + p1->Y ();

						node.set_coordinates (c);
						if (std::find (nodes.begin (), nodes.end (), node) == nodes.end ())
						{
							nodes.push_back (node);
							n_nodes++;
						}

						t0.push_back (node);

						if (t == 0) // if current side is left, push node there
						{
							sides[0][sl].nodes.push_back (node);
						}
						if (t == sides[1][su].sections) // if current side is right, push node there
						{
							sides[0][sr].nodes.push_back (node);
						}
					}
					node.set_Point (*p2);
					t0.push_back (node);

					if (t == 0) // if current side is left, push up left node 
					{
						sides[0][sl].nodes.push_back (node);
					}
					if (t == sides[1][su].sections) // if current side is right, push up right node there
					{
						sides[0][sr].nodes.push_back (node);
					}

					if (t > 0)
					{
						set_triangles (area, or);
					}
					t1.clear ();
					t1.insert (t1.begin (), t0.begin (), t0.end ());
				}
			}
			else
			{
				if (sides[0][sr].sections == sides[0][sl].sections) // right and left become the basis
				{
					// do left and right side fully
					// direction should technically always be up the axis

					// left
					// get length of the left side
					L = sqrt (pow (tetra_nodes[i][j].X () - tetra_nodes[i][j + 1].X (), 2.0) +
						pow (tetra_nodes[i][j].Y () - tetra_nodes[i][j + 1].Y (), 2.0));

					// get coef for the side
					q = sides[0][sl].coef;
					// get geometric sum
					if (fabs (q - 1.0) < ZERO_DIFFERENCE_Mixed_Mesh) // or if coef == 1, just the amount of sections
						GPS = sides[0][sl].sections;
					else
						GPS = (1.0 - pow (q, sides[0][sl].sections)) / (1.0 - q);

					c[0] = tetra_nodes[i][j].X ();
					c[1] = tetra_nodes[i][j].Y ();
					node.set_coordinates (c);
					if (i == 0)
						sides[0][sl].nodes.push_back (node);
					if (std::find (nodes.begin (), nodes.end (), node) == nodes.end ())
					{
						nodes.push_back (node);
						n_nodes++;
					}

					// unit vector on the side
					unit[0] = (tetra_nodes[i][j + 1].X () - tetra_nodes[i][j].X ()) / L;
					unit[1] = (tetra_nodes[i][j + 1].Y () - tetra_nodes[i][j].Y ()) / L;

					// get length of the first section, again, always the down one for the left side
					l0 = L / GPS;
					for (int k = 1; k < sides[0][sl].sections + 1; k++) // push all the nodes 
					{
						if (fabs (q - 1.0) < ZERO_DIFFERENCE_Mixed_Mesh)
							l = l0 * k;
						else
							l = l0 * (1.0 - pow (q, k)) / (1.0 - q);

						c[0] = unit[0] * l + tetra_nodes[i][j].X ();
						c[1] = unit[1] * l + tetra_nodes[i][j].Y ();

						node.set_coordinates (c);
						if (std::find (nodes.begin (), nodes.end (), node) == nodes.end ())
						{
							nodes.push_back (node);
							n_nodes++;
						}
						if (i == 0)
							sides[0][sl].nodes.push_back (node);
					}
					// left side should *technically be done
					// do the right one then
					// get length of the right side
					L = sqrt (pow (tetra_nodes[i + 1][j + 1].X () - tetra_nodes[i + 1][j].X (), 2.0) +
						pow (tetra_nodes[i + 1][j + 1].Y () - tetra_nodes[i + 1][j].Y (), 2.0));

					// get coef for the side
					q = sides[0][sr].coef;
					// get geometric sum
					if (fabs (q - 1.0) < ZERO_DIFFERENCE_Mixed_Mesh) // or if coef == 1, just the amount of sections
						GPS = sides[0][sr].sections;
					else
						GPS = (1.0 - pow (q, sides[0][sr].sections)) / (1.0 - q);

					// starting point, always down one for the right side, so has to be pushed always
					c[0] = tetra_nodes[i + 1][j].X ();
					c[1] = tetra_nodes[i + 1][j].Y ();
					node.set_coordinates (c);
					if (std::find (nodes.begin (), nodes.end (), node) == nodes.end ())
					{
						nodes.push_back (node);
						n_nodes++;
					}
					sides[0][sr].nodes.push_back (node);

					// unit vector on the side
					unit[0] = (tetra_nodes[i + 1][j + 1].X () - tetra_nodes[i + 1][j].X ()) / L;
					unit[1] = (tetra_nodes[i + 1][j + 1].Y () - tetra_nodes[i + 1][j].Y ()) / L;

					// get length of the first section, again, always the left one for the down side
					l0 = L / GPS;
					for (int k = 1; k < sides[0][sr].sections + 1; k++) // push all the nodes 
					{
						if (fabs (q - 1.0) < ZERO_DIFFERENCE_Mixed_Mesh) // or if coef == 1, just the amount of sections
							l = l0 * k;
						else
							l = l0 * (1.0 - pow (q, k)) / (1.0 - q);

						c[0] = unit[0] * l + tetra_nodes[i + 1][j].X ();
						c[1] = unit[1] * l + tetra_nodes[i + 1][j].Y ();

						node.set_coordinates (c);
						if (std::find (nodes.begin (), nodes.end (), node) == nodes.end ())
						{
							nodes.push_back (node);
							n_nodes++;
						}
						sides[0][sr].nodes.push_back (node);
					}
					// go through nodes on the up and down side
					for (int t = 0; t < sides[0][sl].sections + 1; t++)
					{
						t0.clear ();
						// figure out starting point and finishing point
						// always go up the side!
						p1 = sides[0][sl].nodes[t].get_Point ();
						p2 = sides[0][sr].nodes[t].get_Point ();
						// count q and n for the side through qs and ns of down and left
						n = (int)((1.0 - (double)(t) / (double)(sides[0][sl].sections)) * sides[1][sd].sections +
							((double)t / (double)(sides[0][sl].sections)) * sides[1][su].sections);
						q = ((1.0 - (double)(t) / (double)(sides[0][sl].sections)) * sides[1][sd].coef +
							((double)t / (double)(sides[0][sl].sections)) * sides[1][su].coef);

						// get length of the side to slice
						L = sqrt (pow (p1->X () - p2->X (), 2.0) +
							pow (p1->Y () - p2->Y (), 2.0));

						// get geometric sum
						if (fabs (q - 1.0) < ZERO_DIFFERENCE_Mixed_Mesh) // or if coef == 1, just the amount of sections
							GPS = n;
						else
							GPS = (1.0 - pow (q, n)) / (1.0 - q);

						// unit vector on the side
						unit[0] = (p2->X () - p1->X ()) / L;
						unit[1] = (p2->Y () - p1->Y ()) / L;

						node.set_Point (*p1);
						t0.push_back (node);

						// get length of the first section, again, always the down one for the left side
						l0 = L / GPS;
						if (t == 0) // if current side is down, push left down node there
						{
							sides[1][sd].nodes.push_back (node);
						}
						if (t == sides[0][sl].sections) // if current side is up, push right up node there
						{
							sides[1][su].nodes.push_back (node);
						}
						for (int k = 1; k < n; k++) // push nodes between basises
						{
							if (fabs (q - 1.0) < ZERO_DIFFERENCE_Mixed_Mesh) // or if coef == 1, just the amount of sections
								l = l0 * k;
							else
								l = l0 * (1.0 - pow (q, k)) / (1.0 - q);
							c[0] = unit[0] * l + p1->X ();
							c[1] = unit[1] * l + p1->Y ();

							node.set_coordinates (c);
							if (std::find (nodes.begin (), nodes.end (), node) == nodes.end ())
							{
								nodes.push_back (node);
								n_nodes++;
							}

							t0.push_back (node);

							if (t == 0) // if current side is down, push node there
							{
								sides[1][sd].nodes.push_back (node);
							}
							if (t == sides[0][sl].sections) // if current side is up, push node there
							{
								sides[1][su].nodes.push_back (node);
							}
						}

						node.set_Point (*p2);
						t0.push_back (node);

						if (t == 0) // if current side is down, push left down node 
						{
							sides[1][sd].nodes.push_back (node);
						}
						if (t == sides[0][sl].sections) // if current side is up, push right up node there
						{
							sides[1][su].nodes.push_back (node);
						}

						if (t > 0)
						{
							set_triangles (area, or);
						}
						t1.clear ();
						t1.insert (t1.begin (), t0.begin (), t0.end ());
					}
				}
				else // what a shame
				{
					return false;
				}
			}
		}
	}
	//delete p1;
	//delete p2;
	return true;
}

void Mixed_Triangular_Mesh::set_triangles (int a, int o)
{
	int n[3]; // nodes that form the triangle
	std::size_t n0, n1; // amount of nodes on each side
	unsigned int cn0, cn1; // current node on each side
	double d0, d1; // distances

	n0 = t0.size ();
	n1 = t1.size ();
	// go by sides
	for (cn0 = 0, cn1 = 0; cn0 < n0 - 1 && cn1 < n1 - 1; )
	{
		// get distance between current node on one side and next node on the other
		d0 = t1[cn1].distance (t0[cn0 + 1]);
		d1 = t0[cn0].distance (t1[cn1 + 1]);

		if (d0 < d1) // make triangle with cn0, cn1 and next t0 node
		{
			n[0] = get_node_number (t0[cn0]);
			n[1] = get_node_number (t1[cn1]);
			n[2] = get_node_number (t0[cn0 + 1]);
			cn0++;
		}
		else // make triangle with cn0, cn1 and next t1 node
		{
			n[0] = get_node_number (t0[cn0]);
			n[1] = get_node_number (t1[cn1]);
			n[2] = get_node_number (t1[cn1 + 1]);
			cn1++;
		}
		if (o == 1) // 1-st order triangle
		{
			std::unique_ptr <Triangle> triangle = std::make_unique<Triangle> ();
			triangle->set_element (n);
			triangle->set_area (a);
			triangle->set_base_nodes (n);
			elements.push_back (std::move(triangle));
			n_elements++;
		}		
		if (o == 2) // hierar
		{
			std::unique_ptr <Triangle_Hier> triangle_H = std::make_unique<Triangle_Hier> ();
			triangle_H->set_element (n);
			triangle_H->set_area (a);
			triangle_H->set_base_nodes (n);
			elements.push_back (std::move (triangle_H));
			n_elements++;
		}


	}

	// if there are still nodes left on one of the sides, make triangles with them like a fan
	while (cn0 < n0 - 1)
	{
		n[0] = get_node_number (t0[cn0]);
		n[1] = get_node_number (t1[cn1]);
		n[2] = get_node_number (t0[cn0 + 1]);
		cn0++;
		if (o == 1) // 1-st order triangle
		{
			std::unique_ptr <Triangle> triangle = std::make_unique<Triangle> ();
			triangle->set_element (n);
			triangle->set_area (a);
			triangle->set_base_nodes (n);
			elements.push_back (std::move (triangle));
			n_elements++;
		}
		if (o == 2) // hierar
		{
			std::unique_ptr <Triangle_Hier> triangle_H = std::make_unique<Triangle_Hier> ();
			triangle_H->set_element (n);
			triangle_H->set_area (a);
			triangle_H->set_base_nodes (n);
			elements.push_back (std::move (triangle_H));
			n_elements++;
		}
	}
	while (cn1 < n1 - 1)
	{
		n[0] = get_node_number (t0[cn0]);
		n[1] = get_node_number (t1[cn1]);
		n[2] = get_node_number (t1[cn1 + 1]);
		cn1++;
		if (o == 1) // 1-st order triangle
		{
			std::unique_ptr <Triangle> triangle = std::make_unique<Triangle> ();
			triangle->set_element (n);
			triangle->set_area (a);
			triangle->set_base_nodes (n);
			elements.push_back (std::move (triangle));
			n_elements++;
		}
		if (o == 2) // hierar
		{
			std::unique_ptr <Triangle_Hier> triangle_H = std::make_unique<Triangle_Hier> ();
			triangle_H->set_element (n);
			triangle_H->set_area (a);
			triangle_H->set_base_nodes (n);
			elements.push_back (std::move (triangle_H));
			n_elements++;
		}
	}
}

void Mixed_Triangular_Mesh::output ()
{
	FILE * log = fopen ("Result Files Extra//log_mesh.txt", "w");

	// print mesh data
	fprintf (log, "%i %i\n", n_nodes, n_elements);
	printf ("Nodes:\t%i\tElements:\t%i\n", n_nodes, n_elements);

	fclose (log);

	log = fopen ("MeshData\\Nodes.txt", "w");

	double coordinates[2];
	for (int i = 0; i < n_nodes; i++)
	{
		nodes[i].get_coordinates (coordinates);
		fprintf (log, "%.16lf %.16lf\n", coordinates[0], coordinates[1]);
	}

	fclose (log);

	log = fopen ("MeshData\\Triangles.txt", "w");

	int nodes[3];
	for (int i = 0; i < n_elements; i++)
	{
		elements[i]->get_def_nodes (nodes);
		fprintf (log, "%i %i %i %i\n", nodes[0], nodes[1], nodes[2], elements[i]->get_area ());
	}

	fclose (log);
}

bool Mixed_Triangular_Mesh::get_isoline_section (int k_element, double * q, double value, double * c1, double * c2)
{
	int amount = elements[k_element]->get_isoline_points (*this, value, q, c1, c2);
	if (amount == 2)
		return true;
	return false;
}

void Mixed_Triangular_Mesh::renumerate ()
{
	std::queue <int> Q;
	std::vector < std::vector <int> > node_links (n_nodes, std::vector <int> (0));
	int d_nodes;
	int * n;

	nodes.size ();
	for (std::size_t i = 0, vsize = elements.size (); i < vsize; i++)
	{
		d_nodes = elements[i]->get_amount_of_def_nodes ();
		n = new int[d_nodes];
		for (int j = 0; j < d_nodes; j++)
		{
			for (int k = 0; k < d_nodes; k++)
			{
				if (k != j)
				{
					elements[i]->get_def_nodes (n);
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
		n_defining = elements[i]->get_amount_of_def_nodes ();
		node_defining = new int[n_defining];
		elements[i]->get_def_nodes (node_defining);

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
			elements[i]->set_nodes (node_defining);
		}
		delete[] node_defining;

		n_defining = elements[i]->get_amount_of_base_nodes ();
		node_defining = new int[n_defining];
		elements[i]->get_base_nodes (node_defining);
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
			elements[i]->set_base_nodes (node_defining);
		}

		elements[i]->sort_def_nodes ();
		delete[] node_defining;
	}

	nodes.clear ();
	nodes.insert (nodes.begin (), NF.begin (), NF.end ());
}

bool Mixed_Triangular_Mesh::need_to_numerate_edges ()
{
	return true;
}

int Mixed_Triangular_Mesh::numerate_functions ()
{
	int counter = 0;
	std::vector <int> edges_w_func;
	int e[16];
	int n_elem;
	bool add_edge;

	// store numbers of edges that have functions on them
	for (size_t i = 0, i_end = edges->get_n_entries (); i < i_end; i++)
	{
		add_edge = true;
		// figure out which edges have functions
		n_elem = get_elements (int (i), e);
		for (int k = 0; k < n_elem && add_edge; k++)
		{
			// if all elements that have the edge are 2-order
			if (elements[e[k]]->get_order () == 1)
				add_edge = false;
		}
		// save its number
		if (add_edge)
			edges_w_func.push_back ((int)i);
	}
	
	int * noded_funcs = new int[n_nodes]; // functions' numbers by nodes
	int * edged_funcs = new int[(int)edges_w_func.size()]; // functions' numbers by edges_w_func

	// go by nodes
	int es, ee;
	int edge_number;
	for (int i = 0; i < n_nodes; i++)
	{
		// node's function is counter
		noded_funcs[i] = counter;
		counter++;
		// get all the edges that have that node as minimal
		edges->get_edges_by_min_node (i, &es, &ee);
		for (int e = es; e < ee; e++)
		{
			auto edge =	std::find (edges_w_func.begin (), edges_w_func.end (), e);
			if (edge != edges_w_func.end())
			{
				// if the edge's number in edges_w_func, counter it
				edge_number = (int)std::distance (edges_w_func.begin (), edge);
				edged_funcs[edge_number] = counter;
				counter++;
			}
		}
	}
	
	// go by elements
	int def_nodes[3];
	int functions[6];
	int n1, n2;
	int gl_edge_number;
	for (int i = 0; i < n_elements; i++)
	{
		elements[i]->get_def_nodes (def_nodes);
		// copy numbers from functions' numbers by nodes
		for (int k = 0; k < 3; k++)
		{
			functions[k] = noded_funcs[def_nodes[k]];
		}

		// if order is 2, also copy numbers from functions' numbers by edges_w_func 
		if (elements[i]->get_order () == 2)
		{
			for (int k = 0; k < 3; k++)
			{
				functions[k + 3] = -1;

				n1 = def_nodes[k / 2];
				n2 = def_nodes[(k + 3) / 2];
				gl_edge_number = edges->get_edge_number (n1, n2);

				auto edge = std::find (edges_w_func.begin (), edges_w_func.end (), gl_edge_number);
				if (edge != edges_w_func.end ())
				{
					edge_number = (int)std::distance (edges_w_func.begin (), edge);
					functions[k + 3] = edged_funcs[edge_number];
				}
			}
		}
		// set element's funcs
		elements[i]->set_global_functions (functions);
	}

	delete[] noded_funcs;
	delete[] edged_funcs;
	// return amount of global functions
	return counter;
}

void Mixed_Triangular_Mesh::fold_mesh ()
{
	// define maximum amount of elements in a fold
	int max_elements = (int)sqrt (n_elements) + 1;
	// start with the whole mesh as a fold
	Fold fold;

	for (int i = 0; i < 2; i++)
	{
		fold.c0[i] = coord0[i];
		fold.cN[i] = coordN[i];
	}
	fold.k_axis = 0;
	folds.push_back (fold);

	double coordinates[2];

	bool added;
	// go through elements
	for (int k_element = 0; k_element < n_elements; k_element++)
	{
		elements[k_element]->get_mass_center (*this, coordinates);

		added = false;
		for (int k_fold = 0, k_fold_end = (int)folds.size (); k_fold < k_fold_end && !added; k_fold++)
		{
			// if a mass center of element lies in the fold, add it to the list
			if ((folds[k_fold].c0[0] - 1e-10 < coordinates[0] && coordinates[0] < folds[k_fold].cN[0] + 1e-10) &&
				(folds[k_fold].c0[1] - 1e-10 < coordinates[1] && coordinates[1] < folds[k_fold].cN[1] + 1e-10))
			{
				folds[k_fold].elements.push_back (k_element);
				added = true;
			}
		}
		// if amount of elements exceeds maximum
		for (int k_fold = 0; k_fold < (int)folds.size (); k_fold++)
		{
			while (folds[k_fold].elements.size () > max_elements)
			{
				// change fold's axis
				folds[k_fold].k_axis = (folds[k_fold].k_axis + 1) % 2;
				// split by the axis
				for (int i = 0; i < 2; i++)
				{
					fold.c0[i] = folds[k_fold].c0[i];
					fold.cN[i] = folds[k_fold].cN[i];
				}
				double cs = (folds[k_fold].cN[folds[k_fold].k_axis] + folds[k_fold].c0[folds[k_fold].k_axis]) / 2.0;
				// one is a changed old one
				folds[k_fold].cN[folds[k_fold].k_axis] = cs;
				// another adds to the list
				fold.c0[folds[k_fold].k_axis] = cs;

				std::vector <int> fold_elements;
				std::vector <int> k_fold_elements;
				// split elements
				for (int i = 0, i_end = (int)folds[k_fold].elements.size (); i < i_end; i++)
				{
					elements[folds[k_fold].elements[i]]->get_mass_center (*this, coordinates);
					if ((folds[k_fold].c0[folds[k_fold].k_axis] - 1e-10 < coordinates[folds[k_fold].k_axis] && coordinates[folds[k_fold].k_axis] < folds[k_fold].cN[folds[k_fold].k_axis] + 1e-10))
					{
						k_fold_elements.push_back (folds[k_fold].elements[i]);
					}
					else
					{
						fold_elements.push_back (folds[k_fold].elements[i]);
					}
				}
				// reset old fold
				folds[k_fold].elements.clear ();
				folds[k_fold].elements.insert (folds[k_fold].elements.begin (), k_fold_elements.begin (), k_fold_elements.end ());
				//set new fold
				fold.elements.clear ();
				fold.elements.insert (fold.elements.begin (), fold_elements.begin (), fold_elements.end ());
				fold.k_axis = folds[k_fold].k_axis + 1;
				// add new fold
				folds.push_back (fold);
			}
		}
	}

	int base_nodes[3];
	// move fold's boundaries to fit triangles
	for (int k_fold = 0, k_fold_end = (int)folds.size (); k_fold < k_fold_end; k_fold++)
	{
		double c0[2], cN[2];
		// current boundaries
		for (int i = 0; i < 2; i++)
		{
			c0[i] = folds[k_fold].c0[i];
			cN[i] = folds[k_fold].cN[i];
		}
		// go through elements
		for (int k_element = 0, end_elements = (int)folds[k_fold].elements.size (); k_element < end_elements; k_element++)
		{
			elements[folds[k_fold].elements[k_element]]->get_base_nodes (base_nodes);
			for (int i = 0; i < 3; i++)
			{
				get_node_coordinates (base_nodes[i], coordinates);
				for (int j = 0; j < 2; j++)
				{
					if (coordinates[j] < c0[j])
						c0[j] = coordinates[j];
					if (coordinates[j] > cN[j])
						cN[j] = coordinates[j];
				}
			}
		}
		// reset boundaries
		for (int i = 0; i < 2; i++)
		{
			folds[k_fold].c0[i] = c0[i];
			folds[k_fold].cN[i] = cN[i];
		}
	}

}

Mixed_Triangular_Mesh::Mixed_Triangular_Mesh ()
{
	dim = 2;

	coord0 = new double[dim];
	coordN = new double[dim];
	n_axis = new int[dim];

	materials = NULL;
	order = NULL;
}

Mixed_Triangular_Mesh::Mixed_Triangular_Mesh (const Mixed_Triangular_Mesh & mixed_Triangular_Mesh)
{
	// if memory has been allocated, free it
	if (coord0 != NULL)
		delete[] coord0;
	if (coordN != NULL)
		delete[] coordN;
	if (n_axis != NULL)
		delete[] n_axis;
	if (materials != NULL)
		delete[] materials;
	// reset dimentionality
	dim = mixed_Triangular_Mesh.dim;

	// allocate the memory
	// and copy Mixed_Triangular_Mesh's data
	coord0 = new double[dim];
	for (int i = 0; i < dim; i++)
		coord0[i] = mixed_Triangular_Mesh.coord0[i];

	coordN = new double[dim];
	for (int i = 0; i < dim; i++)
		coordN[i] = mixed_Triangular_Mesh.coordN[i];

	n_axis = new int[dim];
	for (int i = 0; i < dim; i++)
		n_axis[i] = mixed_Triangular_Mesh.n_axis[i];

	materials = new int[n_axis[0] * n_axis[1]];
	for (int i = 0; i < n_axis[0] * n_axis[1]; i++)
		materials[i] = mixed_Triangular_Mesh.materials[i];

	// clear nodes and copy Mixed_Triangular_Mesh's nodes
	nodes.clear ();
	n_nodes = mixed_Triangular_Mesh.n_nodes;
	nodes.insert (nodes.begin (), mixed_Triangular_Mesh.nodes.begin (), mixed_Triangular_Mesh.nodes.end ());

	// clear elements and copy Mixed_Triangular_Mesh's elements
	elements.clear ();
	n_elements = mixed_Triangular_Mesh.n_elements;
	for (const auto& e : mixed_Triangular_Mesh.elements)
		elements.push_back (std::make_unique<Element> (*e));
}

Mixed_Triangular_Mesh::~Mixed_Triangular_Mesh ()
{
	if (materials != NULL)
		delete[] materials;
	if (order != NULL)
		delete[] order;
}



