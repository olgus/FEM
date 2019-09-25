#include "Triangular_Mesh_Delaunay.h"

void NS_Triangular_Mesh_Delaunay::Triangle_data::get_edge (int k, int * n)
{
	if (k < 2)
		n[0] = nodes[0];
	else
		n[0] = nodes[1];
	if (k > 0)
		n[1] = nodes[2];
	else
		n[1] = nodes[1];
}

bool NS_Triangular_Mesh_Delaunay::Triangle_data::has_edge (int n0, int n1)
{
	bool match[2] = { false, false };
	for (int i = 0; i < 3; i++)
	{
		if (nodes[i] == n0)
			match[0] = true;
		if (nodes[i] == n1)
			match[1] = true;
	}
	return match[0] && match[1];
}

int NS_Triangular_Mesh_Delaunay::Triangle_data::get_neighbour (int * n)
{
	if (n[0] == nodes[0])
	{
		if (n[1] == nodes[1])
			return neighbours[0];
		else
			return neighbours[1];
	}
	else
		return neighbours[2];
}

int NS_Triangular_Mesh_Delaunay::Triangle_data::get_neighbour (int k)
{
	return neighbours[k];
}

void NS_Triangular_Mesh_Delaunay::Triangle_data::replace_neighbour (int * n, int new_neighbour)
{
	if (n[0] == nodes[0])
	{
		if (n[1] == nodes[1])
		{
			neighbours[0] = new_neighbour;
		}
		else
		{
			neighbours[1] = new_neighbour;
		}
	}
	else
	{
		neighbours[2] = new_neighbour;
	}
}

void NS_Triangular_Mesh_Delaunay::Triangle_data::set_nodes (int * n)
{
	for (int i = 0; i < 3; i++)
		nodes[i] = n[i];
	std::sort (nodes, nodes + 3);
}

NS_Triangular_Mesh_Delaunay::Triangle_data::Triangle_data ()
{
	for (int i = 0; i < 3; i++)
	{
		nodes[i] = -1;
		neighbours[i] = -1;
	}
}

bool NS_Triangular_Mesh_Delaunay::Tree_of_triangles::SS_intersect (double * p1, double * p2, double * p3, double * p4)
{
	double k = (p4[1] - p3[1]) * (p2[0] - p1[0]) - (p4[0] - p3[0]) * (p2[1] - p1[1]);
	double k1 = (p4[0] - p3[0]) * (p1[1] - p3[1]) - (p4[1] - p3[1]) * (p1[0] - p3[0]);
	double k2 = (p2[0] - p1[0]) * (p1[1] - p3[1]) - (p2[1] - p1[1]) * (p1[0] - p3[0]);
	if (fabs (k) < PRECIS) // parallel
	{
		if (fabs (k1) < PRECIS) // coincide
		{
			if (
				(((p3[0] - PRECIS < p1[0]) && (p1[0] < p4[0] + PRECIS)) && ((p3[1] - PRECIS < p1[1]) && (p1[1] < p4[1] + PRECIS))) ||
				(((p3[0] - PRECIS < p2[0]) && (p2[0] < p4[0] + PRECIS)) && ((p3[1] - PRECIS < p2[1]) && (p2[1] < p4[1] + PRECIS))))
				return true;
		}
	}
	else
	{
		double c1 = k1 / k;
		double c2 = k2 / k;
		if (((0 - PRECIS < c1) && (c1 < 1 + PRECIS)) &&
			((0 - PRECIS < c2) && (c2 < 1 + PRECIS)))
			return true;
	}
	return false;
}

bool NS_Triangular_Mesh_Delaunay::Tree_of_triangles::TR_intersect (const Rectangle & rect, int k_triangle)
{
	double c[3][2];
	tmd->get_triangle_coordinates (k_triangle, &c);
	for (int i = 0; i < 3; i++)
	{
		// if at least one triangle node is in the rectangle, intersection exists
		if (((rect.c0[0] - PRECIS < c[i][0]) && (c[i][0] < rect.cN[0] + PRECIS)) &&
			((rect.c0[1] - PRECIS < c[i][1]) && (c[i][1] < rect.cN[1] + PRECIS)))
			return true;
	}
	// otherwise check sides intersections
	double p[4][2];
	p[0][0] = rect.c0[0];
	p[0][1] = rect.c0[1];
	p[1][0] = rect.c0[0];
	p[1][1] = rect.cN[1];
	p[2][0] = rect.cN[0];
	p[2][1] = rect.cN[1];
	p[3][0] = rect.cN[0];
	p[3][1] = rect.c0[1];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (SS_intersect (c[i], c[(i + 1) % 3], p[j], p[(j + 1) % 4]))
				return true;
		}
	}
	return false;
}

void NS_Triangular_Mesh_Delaunay::Tree_of_triangles::search_to_add (Tree_Node * tree_node, int k_triangle, std::vector<Tree_Node*>* tree_nodes)
{
	bool CR_intersection[] = { false, false, false, false };
	int CR_intersection_counter;
	// check how many intersections
	CR_intersection_counter = 0;
	for (int k = 0; k < 4; k++)
	{
		if (tree_node->next[k] != NULL)
		{
			CR_intersection[k] = TR_intersect (tree_node->next[k]->rect, k_triangle);
			if (CR_intersection[k])
				CR_intersection_counter++;
		}
	}
	if (CR_intersection_counter == 0)
		tree_nodes->push_back (tree_node);
	// if more than 2, node found
	if (CR_intersection_counter >= 3)
	{
		tree_nodes->push_back (tree_node);
	}
	else
	{
		// otherwise go into those nodes
		for (int k = 0; k < 4; k++)
		{
			if (CR_intersection[k])
				search_to_add (tree_node->next[k], k_triangle, tree_nodes);
		}
	}
}

void NS_Triangular_Mesh_Delaunay::Tree_of_triangles::search (Tree_Node * tree_node, int k_triangle, std::vector<Tree_Node*>* tree_nodes)
{
	bool CR_intersection[] = { false, false, false, false };
	int CR_intersection_counter;
	// check if we have found the triangle
	if (std::find (tree_node->triangles.begin (), tree_node->triangles.end (), k_triangle) != tree_node->triangles.end ())
		tree_nodes->push_back (tree_node);

	// check how many intersections
	CR_intersection_counter = 0;
	for (int k = 0; k < 4; k++)
	{
		if (tree_node->next[k] != NULL)
		{
			CR_intersection[k] = TR_intersect (tree_node->next[k]->rect, k_triangle);
			if (CR_intersection[k])
				CR_intersection_counter++;
		}
	}
	// in that case triangle will not be written in child nodes
	if (CR_intersection_counter >= 3)
		return;

	// otherwise go into those nodes
	for (int k = 0; k < 4; k++)
	{
		if (CR_intersection[k])
			search (tree_node->next[k], k_triangle, tree_nodes);
	}
}

int NS_Triangular_Mesh_Delaunay::Tree_of_triangles::search (Tree_Node * tree_node, double * point)
{
	bool belong;
	int r = -1;
	// if there are any triangle, check them
	for (int k = 0; k < (int)tree_node->triangles.size (); k++)
	{
		belong = tmd->triangle_point_inside (tree_node->triangles[k], point);
		if (belong)
			return tree_node->triangles[k];
	}

	// go through leaves 
	// find rectangles that have the point (2 or 4 in case of point lying on the boundary)
	for (int k = 0; k < 4 && r == -1; k++)
	{
		if (tree_node->next[k] != NULL)
		{
			belong = tree_node->next[k]->rect.point_inside (point);
			if (belong)
				r = search (tree_node->next[k], point);
		}
	}
	// if no leaves and point doesn't belong to triangle return -1
	return r;
}

void NS_Triangular_Mesh_Delaunay::Tree_of_triangles::delete_child_nodes (Tree_Node * tree_node)
{
	for (int k = 0; k < 4; k++)
	{
		if (tree_node->next[k] != NULL)
		{
			delete_child_nodes (tree_node->next[k]);
			delete tree_node->next[k];
		}
	}
}

void NS_Triangular_Mesh_Delaunay::Tree_of_triangles::initialize (double C0[2], double CN[2], TMD_pointer * TMD)
{
	root = new Tree_Node;
	tmd = TMD;

	for (int i = 0; i < 4; i++)
		root->next[i] = NULL;
	root->rect = Rectangle{ C0[0], C0[1], CN[0], CN[1] };
}

void NS_Triangular_Mesh_Delaunay::Tree_of_triangles::add_triangle (int k_triangle)
{
	std::vector <Tree_Node *> tree_nodes;
	search_to_add (root, k_triangle, &tree_nodes);
	for (int k = 0, k_end = (int)tree_nodes.size (); k < k_end; k++)
	{
		tree_nodes[k]->triangles.push_back (k_triangle);
		// if amount > MAX_LISTING_SIZE, divide the node
		if ((int)tree_nodes[k]->triangles.size () > MAX_LISTING_SIZE)
			node_division (tree_nodes[k]);
	}
}

void NS_Triangular_Mesh_Delaunay::Tree_of_triangles::reset_triangle (int k_triangle)
{
	std::vector <Tree_Node *> tree_nodes;
	// find nodes that have the triangle
	search (root, k_triangle, &tree_nodes);
	// remove the triangle
	for (int k = 0, k_end = (int)tree_nodes.size (); k < k_end; k++)
	{
		tree_nodes[k]->triangles.erase (std::remove (tree_nodes[k]->triangles.begin (), tree_nodes[k]->triangles.end (), k_triangle), tree_nodes[k]->triangles.end ());
	}
	// find the new place for it
	add_triangle (k_triangle);
}

int NS_Triangular_Mesh_Delaunay::Tree_of_triangles::tree_search (double * point)
{
	bool belong;
	// check that the point is in the root rectangle
	if (root != NULL)
	{
		belong = root->rect.point_inside (point);
		if (belong)
			return search (root, point);
	}
	return -1;
}

bool NS_Triangular_Mesh_Delaunay::Rectangle::point_inside (double * point)
{
	if ((c0[0] - PRECIS < point[0]) && (point[0] < cN[0] + PRECIS) &&
		(c0[1] - PRECIS < point[1]) && (point[1] < cN[1] + PRECIS))
		return true;
	return false;
}

void NS_Triangular_Mesh_Delaunay::Tree_of_triangles::node_division (Tree_Node * tree_node)
{
	// if node was already divided then leave it
	if (tree_node->next[0] != NULL)
		return;
	// divide node geometrically making 4 new rectangles
	for (int i = 0; i < 4; i++)
	{
		tree_node->next[i] = new Tree_Node;
		for (int k = 0; k < 4; k++)
		{
			tree_node->next[i]->next[k] = NULL;
		}
	}
	double div[] = { (tree_node->rect.c0[0] + tree_node->rect.cN[0]) / 2.0 , (tree_node->rect.c0[1] + tree_node->rect.cN[1]) / 2.0 };
	tree_node->next[0]->rect.c0[0] = tree_node->rect.c0[0];
	tree_node->next[0]->rect.c0[1] = div[1];
	tree_node->next[0]->rect.cN[0] = div[0];
	tree_node->next[0]->rect.cN[1] = tree_node->rect.cN[1];

	tree_node->next[1]->rect.c0[0] = div[0];
	tree_node->next[1]->rect.c0[1] = div[1];
	tree_node->next[1]->rect.cN[0] = tree_node->rect.cN[0];
	tree_node->next[1]->rect.cN[1] = tree_node->rect.cN[1];

	tree_node->next[2]->rect.c0[0] = div[0];
	tree_node->next[2]->rect.c0[1] = tree_node->rect.c0[1];
	tree_node->next[2]->rect.cN[0] = tree_node->rect.cN[0];
	tree_node->next[2]->rect.cN[1] = div[1];

	tree_node->next[3]->rect.c0[0] = tree_node->rect.c0[0];
	tree_node->next[3]->rect.c0[1] = tree_node->rect.c0[1];
	tree_node->next[3]->rect.cN[0] = div[0];
	tree_node->next[3]->rect.cN[1] = div[1];

	bool CR_intersection[] = { false, false, false, false };
	int counter;
	// go through triangle in the original node
	for (int i = 0; i < (int)tree_node->triangles.size (); i++)
	{
		// calc relationships between new rectangles and the Triangle
		counter = 0;
		for (int k = 0; k < 4; k++)
		{
			CR_intersection[k] = TR_intersect (tree_node->next[k]->rect, tree_node->triangles[i]);
			if (CR_intersection[k])
				counter++;
		}
		// if the Triangle belongs to the 3+ rectangles, keep it the original node
		// otherwise put into respective new nodes
		// zero out the Triangle number in the original node
		if (counter <= 2)
		{
			for (int k = 0; k < 4; k++)
			{
				if (CR_intersection[k])
				{
					tree_node->next[k]->triangles.push_back (tree_node->triangles[i]);
				}
			}
			tree_node->triangles[i] = -1;
		}
	}

	// delete all "-1"s
	tree_node->triangles.erase (std::remove (tree_node->triangles.begin (), tree_node->triangles.end (), -1), tree_node->triangles.end ());

	for (int i = 0; i < 4; i++)
	{
		if ((int)tree_node->next[i]->triangles.size () > MAX_LISTING_SIZE)
			node_division (tree_node->next[i]);
	}
}

void NS_Triangular_Mesh_Delaunay::Tree_of_triangles::output ()
{
	FILE * file_output = fopen ("Source Files//Triangulation//tree.txt", "w");
	Tree_Node * current;
	std::queue <Tree_Node *> queue;
	queue.push (root);

	while (!queue.empty ())
	{
		current = queue.front ();
		queue.pop ();

		if (current != NULL)
		{
			fprintf (file_output, "(%.7lf, %.7lf) ; (%.7lf, %.7lf) : ", current->rect.c0[0], current->rect.c0[1], current->rect.cN[0], current->rect.cN[1]);
			for (int i = 0; i < (int)current->triangles.size (); i++)
				fprintf (file_output, "%i ", current->triangles[i]);
			fprintf (file_output, "\n");

			for (int i = 0; i < 4; i++)
				queue.push (current->next[i]);
		}
	}

	fclose (file_output);
}

NS_Triangular_Mesh_Delaunay::Tree_of_triangles::Tree_of_triangles ()
{
	root = NULL;
	tmd = NULL;
}

NS_Triangular_Mesh_Delaunay::Tree_of_triangles::~Tree_of_triangles ()
{
	if (root != NULL)
	{
		delete_child_nodes (root);
		delete root;
	}
}

void NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::renumerate ()
{
}

bool NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::need_to_numerate_edges ()
{
	return true;
}

void NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::make_elements ()
{
	coord0[0] = coord0[1] = 1e+10;
	coordN[0] = coordN[1] = -1e+10;
	// basic figures nodes
	for (int i = 0, i_end = pg->n_basic_figures; i < i_end; i++)
	{
		// push new nodes
		int n_points = pg->basic_figures[i]->get_amount_of_base_points ();
		int * figure_nodes = new int[n_points];
		for (int k = 0; k < n_points; k++)
		{
			Point p = pg->basic_figures[i]->get_point (k);
			// check if it exists
			Node_2D node;
			double coord[2] = { p.x, p.y };
			node.set_node (coord);
			std::vector <Node_2D>::iterator it = std::find (nodes.begin (), nodes.end (), node);
			if (it != nodes.end ())
			{
				figure_nodes[k] = (int)std::distance (nodes.begin (), it);
			}
			else
			{
				figure_nodes[k] = (int)nodes.size ();
				nodes.push_back (node);
				if (p.x < coord0[0])
					coord0[0] = p.x;
				if (p.y < coord0[1])
					coord0[1] = p.y;
				if (p.x > coordN[0])
					coordN[0] = p.x;
				if (p.y > coordN[1])
					coordN[1] = p.y;
			}
		}
		// add the relation
		BF_triangles bf_tr;
		bf_triangles.push_back (bf_tr);

		// make triangles on them
		{
			int local_nodes[2][3];
			int tr_nodes[3];
			int amount;
			int tr_start = (int)ib_triangles.size ();
			pg->basic_figures[i]->triangulate (&amount, local_nodes);
			for (int k = 0; k < amount; k++)
			{
				for (int j = 0; j < 3; j++)
				{
					tr_nodes[j] = figure_nodes[local_nodes[k][j]];
				}
				Triangle_data triangle;
				triangle.set_nodes (tr_nodes);
				ib_triangles.push_back (triangle);
			}
			// set neighbours between those triangles
			{
				int n[2];
				for (int k_triangle = tr_start, k_triangle_end = (int)ib_triangles.size (); k_triangle < k_triangle_end; k_triangle++)
				{
					for (int k = 0; k < 3; k++)
					{
						ib_triangles[k_triangle].get_edge (k, n);
						ib_triangles[k_triangle].replace_neighbour (n, find_neighbour (tr_start, (int)ib_triangles.size (), k_triangle, n[0], n[1]));
					}
					bf_triangles[i].triangles.push_back (k_triangle);
				}
			}
		}
		delete[] figure_nodes;
	}
	TMD_pointer * tmd = this;
	tot = new Tree_of_triangles;
	tot->initialize (coord0, coordN, tmd);
	for (int i = 0, i_end = (int)ib_triangles.size (); i < i_end; i++)
	{
		tot->add_triangle (i);
	}

	// set neighbours for those triangles
	//{
	//	int n[2];
	//	for (int i = 0, i_end = (int)ib_triangles.size (); i < i_end; i++)
	//	{
	//		for (int k = 0; k < 3; k++)
	//		{
	//			ib_triangles[i].get_edge (k, n);
	//			ib_triangles[i].replace_neighbour (n, find_neighbour (i, n[0], n[1]));
	//		}
	//	}
	//}
	// go through generated points
	{
		int k_triangle;
		double c[2];
		double dist, l_dist;
		for (int i = 0, i_end = (int)pg->gen_points.size (); i < i_end; i++)
		{
			double point[2] = { pg->gen_points[i].x, pg->gen_points[i].y };
			// check if isn't too close to other nodes
			if (!closeness (point))
			{
				// if it isn't, add it
				Node_2D node;
				node.set_coordinates (point);
				nodes.push_back (node);

				// if the point is in the basic figure
				for (int k = 0; k < pg->n_basic_figures; k++)
				{
					Basic_Geometry::Point p = { point[0], point[1] };
					if (pg->basic_figures[k]->point_inside (p))
					{
						// find triangle within that basic figure
						int tr_split = point_in_within_bf (k, point);
						// split it
						split_triangle ((int)nodes.size () - 1, tr_split, k);
					}
				}
			}

			// old version
			//k_triangle = tot->tree_search (point);
			//// don't add new node if it is too close to any of the triangle's nodes
			//dist = 1e+10;
			//if (k_triangle != -1)
			//{
			//	for (int k = 0; k < 3; k++)
			//	{
			//		nodes[ib_triangles[k_triangle].nodes[k]].get_coordinates (c);
			//		Point p = { c[0], c[1] };
			//		l_dist = pg->gen_points[i].distance (p);
			//		if (l_dist < dist)
			//			dist = l_dist;
			//	}
			//	if (dist > MIN_DISTANCE)
			//	{
			//		// add node
			//		Node_2D node;
			//		node.set_coordinates (point);
			//		nodes.push_back (node);
			//		// split triagnle
			//		split_triangle ((int)nodes.size () - 1, k_triangle);
			//	}
			//}
		}
	}
	int appr_n_elements = 0;
	{
		for (int i = 0, i_end = (int)bf_triangles.size (); i < i_end; i++)
		{
			appr_n_elements += (int)bf_triangles[i].triangles.size ();
		}
	}
	// test purposes
	for (int i = 0; i < 0/*(int)(appr_n_elements / 250 + 1)*/; i++)
	{
		printf ("%i ", i);
		test ();
	}
	{
		// make elements 
		for (int i = 0, i_end = (int)bf_triangles.size (); i < i_end; i++)
		{
			for (int k = 0, k_end = (int)bf_triangles[i].triangles.size (); k < k_end; k++)
			{
				std::unique_ptr <Triangle> triangle = std::make_unique<Triangle> ();
				triangle->set_base_nodes (ib_triangles[bf_triangles[i].triangles[k]].nodes);
				triangle->set_nodes (ib_triangles[bf_triangles[i].triangles[k]].nodes);
				triangle->set_area (pg->basic_figures[i]->get_material ());
				elements.push_back (std::move (triangle));
			}
		}
		//// make elements 
		//for (int i = 0, i_end = (int)ib_triangles.size (); i < i_end; i++)
		//{
		//	std::unique_ptr <Triangle> triangle = std::make_unique<Triangle> ();
		//	triangle->set_base_nodes (ib_triangles[i].nodes);
		//	triangle->set_nodes (ib_triangles[i].nodes);
		//	triangle->set_area (find_area (i));
		//	elements.push_back (std::move (triangle));
		//}
	}
	n_elements = (int)elements.size ();
	n_nodes = (int)nodes.size ();

	tot->output ();
}

int NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::find_area (int k_triangle)
{
	// get_center
	double c[3][2];
	get_triangle_coordinates (k_triangle, &c);
	double center[2] = { 0.0, 0.0 };
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 2; j++)
			center[j] += c[i][j];
	for (int j = 0; j < 2; j++)
		center[j] /= 3.0;
	// find belonging basic figure
	for (int i = 0, i_end = pg->n_basic_figures; i < i_end; i++)
	{
		Point p = { center[0], center[1] };
		if (pg->basic_figures[i]->point_inside (p))
			return pg->basic_figures[i]->get_material ();
	}
	return 0;
}

void NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::split_triangle (int k_node, int k_triangle, int bf)
{
	if (k_triangle == -1)
		return;
	// build neighbourhood
	// fuck set
	std::vector <int> neighbourhood;
	neighbourhood.push_back (k_triangle);
	//for (int i = 0; i < 3; i++)
	//{
	//	// first level
	//	if (point_sphere_relation (k_node, ib_triangles[k_triangle].neighbours[i]))
	//	{
	//		neighbourhood.push_back (ib_triangles[k_triangle].neighbours[i]);
	//		for (int j = 0; j < 3; j++)
	//		{
	//			// second level
	//			if (point_sphere_relation (k_node, ib_triangles[ib_triangles[k_triangle].neighbours[i]].neighbours[j]))
	//			{
	//				neighbourhood.push_back (ib_triangles[ib_triangles[k_triangle].neighbours[i]].neighbours[j]);
	//			}
	//		}
	//	}
	//}

	struct neibourhood_relation
	{
		int came_from;
		int triangle;
	};

	std::queue<neibourhood_relation> neighbour_queue;
	neibourhood_relation cur_triangle;
	for (int i = 0; i < 3; i++)
	{
		neighbour_queue.push ({ k_triangle, ib_triangles[k_triangle].neighbours[i] });
		while (!neighbour_queue.empty ())
		{
			cur_triangle = neighbour_queue.front ();
			neighbour_queue.pop ();
			// if point in in the sphere
			if (point_sphere_relation (k_node, cur_triangle.triangle))
			{
				// keep triangle to split
				neighbourhood.push_back (cur_triangle.triangle);
				// add neighbours except came_from
				for (int k = 0; k < 3; k++)
				{
					if (ib_triangles[cur_triangle.triangle].neighbours[k] != cur_triangle.came_from)
						neighbour_queue.push ({ cur_triangle.triangle, ib_triangles[cur_triangle.triangle].neighbours[k] });
				}
			}
		}
	}

	std::sort (neighbourhood.begin (), neighbourhood.end ());
	neighbourhood.erase (unique (neighbourhood.begin (), neighbourhood.end ()), neighbourhood.end ());
	// edges struct
	struct Outside_Edge
	{
		int n[2];
		int neighbour;
		Outside_Edge ()
		{
			n[0] = n[1] = -1;
			neighbour = -1;
		};
		Outside_Edge (int N[2], int Neighbour)
		{
			n[0] = N[0];
			n[1] = N[1];
			neighbour = Neighbour;
		};
		void get_nodes (int * N)
		{
			N[0] = n[0];
			N[1] = n[1];
		}
		bool connected (const Outside_Edge & edge)
		{
			if ((n[0] == edge.n[1]) || ((n[1] == edge.n[0])) || (n[1] == edge.n[1]) || ((n[0] == edge.n[0])))
				return true;
			return false;
		};
		Outside_Edge & operator=(const Outside_Edge & edge)
		{
			n[0] = edge.n[0];
			n[1] = edge.n[1];
			neighbour = edge.neighbour;
			return *this;
		};
	};
	std::vector<Outside_Edge> outside_edges;
	// get outside edges
	{
		int n[2];
		for (int i = 0, i_end = (int)neighbourhood.size (); i < i_end; i++)
		{
			for (int k = 0; k < 3; k++)
			{
				if (std::find (neighbourhood.begin (), neighbourhood.end (), ib_triangles[neighbourhood[i]].neighbours[k]) == neighbourhood.end ())
				{
					ib_triangles[neighbourhood[i]].get_edge (k, n);
					outside_edges.push_back (Outside_Edge (n, ib_triangles[neighbourhood[i]].neighbours[k]));
				}
			}
		}
	}
	// sort them in a circle
	{
		int replace_pos;
		for (int i = 0, i_end = (int)outside_edges.size () - 1; i < i_end; i++)
		{
			replace_pos = i;
			for (int j = i + 1, j_end = (int)outside_edges.size (); j < j_end && replace_pos == i; j++)
			{
				if (outside_edges[i].connected (outside_edges[j]))
					replace_pos = j;
			}
			Outside_Edge e;
			e = outside_edges[i + 1];
			outside_edges[i + 1] = outside_edges[replace_pos];
			outside_edges[replace_pos] = e;

		}
	}
	// define numbers of new triangles
	std::vector<int> new_triangles;
	{
		int tr_nodes[3];
		tr_nodes[2] = k_node;
		int counter = 0;
		for (int i = 0, i_end = (int)outside_edges.size (); i < i_end; i++)
		{
			outside_edges[i].get_nodes (tr_nodes);
			// if degenerate triangle
			if (triangle_degenerate (tr_nodes))
			{
				// no triangle
				new_triangles.push_back (-1);
			}
			else
			{
				// otherwise reset old ones
				if (counter < (int)neighbourhood.size ())
				{
					ib_triangles[neighbourhood[counter]].set_nodes (tr_nodes);
					new_triangles.push_back (neighbourhood[counter]);
					counter++;
				}
				else
				{
					// or add a new one
					new_triangles.push_back ((int)ib_triangles.size ());
					ib_triangles.push_back (Triangle_data ());
					ib_triangles[new_triangles[i]].set_nodes (tr_nodes);
				}

				// reset triangle in a tree
				tot->reset_triangle (new_triangles[i]);
			}
		}
	}
	// fill them 
	{
		int nt;
		int tr_nodes[3];
		int next_tr_nodes[3];
		int connection_edge[2];
		int neighbouring_node;
		tr_nodes[2] = k_node;
		next_tr_nodes[2] = k_node;
		for (int i = 0, i_end = (int)new_triangles.size (); i < i_end; i++)
		{
			// next triangle in a fan
			i == (int)new_triangles.size () - 1 ? nt = 0 : nt = i + 1;
			// get nodes of the next triangle
			outside_edges[nt].get_nodes (next_tr_nodes);
			outside_edges[i].get_nodes (tr_nodes);
			// find second node that connects them
			neighbouring_node = find_neighbouring_node (tr_nodes, next_tr_nodes, k_node);
			// k_node technically should always bigger number than neighbouring_node
			connection_edge[0] = neighbouring_node;
			connection_edge[1] = k_node;

			if (new_triangles[i] != -1)
			{
				// keep the neighbour next to the outside edge
				ib_triangles[new_triangles[i]].replace_neighbour (tr_nodes, outside_edges[i].neighbour);
				ib_triangles[new_triangles[i]].replace_neighbour (connection_edge, new_triangles[nt]);
			}
			// set it as a neighbour anyway
			if (new_triangles[nt] != -1)
			{
				ib_triangles[new_triangles[nt]].replace_neighbour (connection_edge, new_triangles[i]);
			}
			// replace neighbour for current edge
			if (outside_edges[i].neighbour != -1)
			{
				ib_triangles[outside_edges[i].neighbour].replace_neighbour (tr_nodes, new_triangles[i]);
			}
		}
	}

	// update bf_triangles
	for (int i = 0, i_end = (int)new_triangles.size (); i < i_end; i++)
	{
		if (new_triangles[i] != -1)
		{
			if (std::find (bf_triangles[bf].triangles.begin (), bf_triangles[bf].triangles.end (), new_triangles[i]) == bf_triangles[bf].triangles.end ())
				bf_triangles[bf].triangles.push_back (new_triangles[i]);
		}
	}
}

bool NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::point_sphere_relation (int k_node, int k_triangle)
{
	if (k_triangle != -1)
	{
		double c[3][2];
		get_triangle_coordinates (k_triangle, &c);
		double p[2];
		nodes[k_node].get_coordinates (p);

		//double A10 = pow (c[1][0] - p[0], 2.0) + pow (c[1][1] - p[1], 2.0);
		//double A20 = pow (c[2][0] - p[0], 2.0) + pow (c[2][1] - p[1], 2.0);
		//double r = pow (c[0][0] - p[0], 2.0) + pow (c[0][1] - p[1], 2.0);
		//r *= ((c[1][0] - p[0]) * (c[2][1] - p[1])) -
		//	((c[1][1] - p[1]) * (c[2][0] - p[0]));
		//r -= (c[0][0] - p[0]) * 
		//	(A10 * (c[2][1] - p[1]) - A20 * (c[1][1] - p[1]));
		//r += (c[0][1] - p[1]) * 
		//	(A10 * (c[2][0] - p[0]) - A20 * (c[1][0] - p[0]));
		double C0 = pow (c[0][0], 2.0) + pow (c[0][1], 2.0);
		double C1 = pow (c[1][0], 2.0) + pow (c[1][1], 2.0);
		double C2 = pow (c[2][0], 2.0) + pow (c[2][1], 2.0);
		double D0 = c[1][0] * c[2][1] - c[1][1] * c[2][0];
		double D1 = C1 * c[2][1] - C2 * c[1][1];
		double D2 = C1 * c[2][0] - C2 * c[1][0];
		double A = c[0][0] * (c[1][1] - c[2][1]) - c[0][1] * (c[1][0] - c[2][0]) + D0;
		double B = C0 * (c[1][1] - c[2][1]) - c[0][1] * (C1 - C2) + D1;
		double C = C0 * (c[1][0] - c[2][0]) - c[0][0] * (C1 - C2) + D2;
		double D = C0 * D0 - c[0][0] * D1 + c[0][1] * D2;
		double r = A * (pow (p[0], 2.0) + pow (p[1], 2.0)) - p[0] * B + p[1] * C - D;
		r *= A > 0 ? 1.0 : -1.0;
		return r < PRECIS;
	}
	return false;
}

int NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::find_neighbour (int k_triangle, int n0, int n1)
{
	for (int i = 0, i_end = (int)ib_triangles.size (); i < i_end; i++)
	{
		if (i != k_triangle)
			if (ib_triangles[i].has_edge (n0, n1))
				return i;
	}
	return -1;
}

int NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::find_neighbour (int tr_start, int tr_end, int k_triangle, int n0, int n1)
{
	for (int i = tr_start; i < tr_end; i++)
	{
		if (i != k_triangle)
			if (ib_triangles[i].has_edge (n0, n1))
				return i;
	}
	return -1;
}

int NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::find_neighbouring_node (int * n0, int * n1, int excluding_node)
{
	for (int i = 0; i < 3; i++)
	{
		if (n0[i] != excluding_node)
			for (int j = 0; j < 3; j++)
			{
				if (n0[i] == n1[j])
					return n0[i];
			}
	}
	return -1;
}

void NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::test ()
{
	struct BA_pairs {
		int t1, t2;
		bool repeat;
	};

	std::vector <BA_pairs> bap;
	FILE * file = fopen ("Source Files//Triangulation//angles.txt", "w");
	fprintf (file, "first run\n");
	for (int i = 0, i_end = (int)ib_triangles.size (); i < i_end; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (bad_angles (i, ib_triangles[i].neighbours[j]))
			{
				fprintf (file, "%i %i\n", i, ib_triangles[i].neighbours[j]);
				bap.push_back ({ i, ib_triangles[i].neighbours[j], false });
			}
		}
	}

	// fix
	{
		// remove entries if they contain repeats
		for (int i = 0, i_end = (int)bap.size (); i < i_end; i++)
		{
			if (!bap[i].repeat)
			{
				for (int j = i + 1, j_end = (int)bap.size (); j < j_end; j++)
				{
					if ((bap[i].t1 == bap[j].t1) || (bap[i].t1 == bap[j].t2))
					{
						bap[j].repeat = true;
					}
					if ((bap[i].t2 == bap[j].t1) || (bap[i].t2 == bap[j].t2))
					{
						bap[j].repeat = true;
					}
				}
			}
		}
		// go through legit pairs
		for (int i = 0, i_end = (int)bap.size (); i < i_end; i++)
		{
			if (!bap[i].repeat)
			{
				// get nodes
				int n_edge_nodes[2];
				int unique_nodes[2];
				tr_config (bap[i].t1, bap[i].t2, n_edge_nodes, unique_nodes);
				std::sort (unique_nodes, unique_nodes + 2);

				// save edges with neighbours
				struct Outside_Edge
				{
					int n[2];
					int neighbour;
					Outside_Edge (int * N, int Neighbour)
					{
						n[0] = N[0];
						n[1] = N[1];
						neighbour = Neighbour;
					};
				};
				std::vector <Outside_Edge> oes;
				{
					int n[2];
					for (int k = 0; k < 3; k++)
					{
						ib_triangles[bap[i].t1].get_edge (k, n);
						if (!((n[0] == n_edge_nodes[0]) && (n[1] == n_edge_nodes[1]))) // not an edge to remove
						{
							oes.push_back (Outside_Edge (n, ib_triangles[bap[i].t1].get_neighbour (k)));
						}
						ib_triangles[bap[i].t2].get_edge (k, n);
						if (!((n[0] == n_edge_nodes[0]) && (n[1] == n_edge_nodes[1]))) // not an edge to remove
						{
							oes.push_back (Outside_Edge (n, ib_triangles[bap[i].t2].get_neighbour (k)));
						}
					}
				}
				// build new nodes of triangles
				{
					int n[3] = { n_edge_nodes[0], unique_nodes[0], unique_nodes[1] };
					ib_triangles[bap[i].t1].set_nodes (n);
				}
				{
					int n[3] = { n_edge_nodes[1], unique_nodes[0], unique_nodes[1] };
					ib_triangles[bap[i].t2].set_nodes (n);
				}
				// reset neighbours
				for (int k = 0, k_end = (int)oes.size (); k < k_end; k++)
				{
					// egde detected
					if (ib_triangles[bap[i].t1].has_edge (oes[k].n[0], oes[k].n[1]))
					{
						ib_triangles[bap[i].t1].replace_neighbour (oes[k].n, oes[k].neighbour);
						if (oes[k].neighbour != -1)
							ib_triangles[oes[k].neighbour].replace_neighbour (oes[k].n, bap[i].t1);

					}
					// egde detected
					if (ib_triangles[bap[i].t2].has_edge (oes[k].n[0], oes[k].n[1]))
					{
						ib_triangles[bap[i].t2].replace_neighbour (oes[k].n, oes[k].neighbour);
						if (oes[k].neighbour != -1)
							ib_triangles[oes[k].neighbour].replace_neighbour (oes[k].n, bap[i].t2);

					}
				}
				ib_triangles[bap[i].t1].replace_neighbour (unique_nodes, bap[i].t2);
				ib_triangles[bap[i].t2].replace_neighbour (unique_nodes, bap[i].t1);

			}
		}
	}

	fprintf (file, "\nsecond run\n");
	for (int i = 0, i_end = (int)ib_triangles.size (); i < i_end; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (bad_angles (i, ib_triangles[i].neighbours[j]))
			{
				fprintf (file, "%i %i\n", i, ib_triangles[i].neighbours[j]);
			}
		}
	}
	fclose (file);
}

bool NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::bad_angles (int k_triangle_1, int k_triangle_2)
{
	//if (neighbours (k_triangle_1, k_triangle_2))
	if (k_triangle_2 != -1)
	{
		double c[3][2];
		double p[2];
		get_triangle_coordinates (k_triangle_1, &c);
		int n = -1;
		for (int i = 0; i < 3 && n == -1; i++)
			if (std::find (ib_triangles[k_triangle_2].nodes, ib_triangles[k_triangle_2].nodes + 3, ib_triangles[k_triangle_1].nodes[i])
				== std::end (ib_triangles[k_triangle_2].nodes)) // node not found
				n = i;
		if (n != 1 && n != -1)
		{
			p[0] = c[1][0];
			p[1] = c[1][1];

			c[1][0] = c[n][0];
			c[1][1] = c[n][1];

			c[n][0] = p[0];
			c[n][1] = p[1];
		}
		n = -1;
		for (int i = 0; i < 3 && n == -1; i++)
			if (std::find (ib_triangles[k_triangle_1].nodes, ib_triangles[k_triangle_1].nodes + 3, ib_triangles[k_triangle_2].nodes[i])
				== std::end (ib_triangles[k_triangle_1].nodes)) // node not found
				n = ib_triangles[k_triangle_2].nodes[i];
		get_node_coordinates (n, p);

		double l01 = sqrt (pow (p[0] - c[0][0], 2.0) + pow (p[1] - c[0][1], 2.0));
		double l03 = sqrt (pow (p[0] - c[2][0], 2.0) + pow (p[1] - c[2][1], 2.0));
		double l23 = sqrt (pow (c[1][0] - c[2][0], 2.0) + pow (c[1][1] - c[2][1], 2.0));
		double l21 = sqrt (pow (c[1][0] - c[0][0], 2.0) + pow (c[1][1] - c[0][1], 2.0));

		double ca = (p[0] - c[0][0]) * (p[0] - c[2][0]) + (p[1] - c[0][1]) * (p[1] - c[2][1]);
		ca /= l01 * l03;
		double cb = (c[1][0] - c[0][0]) * (c[1][0] - c[2][0]) + (c[1][1] - c[0][1]) * (c[1][1] - c[2][1]);
		cb /= l21 * l23;
		if ((ca < -PRECIS) && (cb < -PRECIS))
			return true;
		if ((ca > -PRECIS) && (cb > -PRECIS))
			return false;
		double sa = (p[0] - c[0][0]) * (p[1] - c[2][1]) - (p[0] - c[2][0]) * (p[1] - c[0][1]);
		sa /= l01 * l03;
		double sb = (c[1][0] - c[2][0]) * (c[1][1] - c[0][1]) - (c[1][0] - c[0][0]) * (c[1][1] - c[2][1]);
		sb /= l21 * l23;
		sa = sqrt (1.0 - ca * ca);
		sb = sqrt (1.0 - cb * cb);
		if ((sa * cb + ca * sb) > PRECIS)
			return false;
		else
			return true;
	}
	return false;
}

bool NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::neighbours (int k_triangle_1, int k_triangle_2)
{
	int counter = 0;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			if (ib_triangles[k_triangle_1].nodes[i] == ib_triangles[k_triangle_2].nodes[j])
				counter++;
	return counter > 1;
}

void NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::tr_config (int k_triangle_1, int k_triangle_2, int * n_edge_nodes, int * unique_nodes)
{
	for (int i = 0; i < 3; i++)
	{
		if (std::find (ib_triangles[k_triangle_1].nodes, ib_triangles[k_triangle_1].nodes + 3, ib_triangles[k_triangle_2].nodes[i]) == std::end (ib_triangles[k_triangle_1].nodes))
			unique_nodes[0] = ib_triangles[k_triangle_2].nodes[i];
	}
	for (int i = 0; i < 3; i++)
	{
		if (std::find (ib_triangles[k_triangle_2].nodes, ib_triangles[k_triangle_2].nodes + 3, ib_triangles[k_triangle_1].nodes[i]) == std::end (ib_triangles[k_triangle_2].nodes))
			unique_nodes[1] = ib_triangles[k_triangle_1].nodes[i];
	}

	int counter = 0;
	for (int i = 0; i < 3; i++)
	{
		if (ib_triangles[k_triangle_2].nodes[i] != unique_nodes[0])
		{
			n_edge_nodes[counter] = ib_triangles[k_triangle_2].nodes[i];
			counter++;
		}
	}
}

bool NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::closeness (double * point)
{
	for (int i = 0, i_end = (int)nodes.size (); i < i_end; i++)
	{
		if (nodes[i].distance (point) < PRECIS)
			return true;
	}
	return false;
}

int NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::point_in_within_bf (int bf, double * point)
{
	for (int i = 0, i_end = (int)bf_triangles[bf].triangles.size (); i < i_end; i++)
	{
		if (triangle_point_inside (bf_triangles[bf].triangles[i], point))
			return bf_triangles[bf].triangles[i];
	}
	return -1;
}

NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay ()
{
	dim = 2;

	coord0 = new double[dim];
	coordN = new double[dim];
	pg = NULL;
	tot = NULL;
}

NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay (const Triangular_Mesh_Delaunay & tmd)
{
}

NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::~Triangular_Mesh_Delaunay ()
{
	if (tot != NULL)
		delete tot;
	if (pg != NULL)
		delete pg;
}

void NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::get_triangle_coordinates (int k_triangle, double (*c)[3][2])
{
	for (int k = 0;k < 3; k++)
		nodes[ib_triangles[k_triangle].nodes[k]].get_coordinates ((*c)[k]);
}

bool NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::triangle_point_inside (int k_triangle, double * point)
{
	double c[3][2];
	// get element's coordinates
	get_triangle_coordinates (k_triangle, &c);

	// get its geometrical area
	double area = (c[0][0] - c[2][0]) * (c[1][1] - c[2][1]) -
		(c[1][0] - c[2][0]) * (c[0][1] - c[2][1]);
	area = fabs (area) / 2.0;

	// calculate geometrical areas of triangles made with point (coordinates)
	double square = 0.0;
	for (int i = 0; i < 3; i++)
	{
		square += fabs ((c[(i + 1) % 3][0] - c[i][0]) * (point[1] - c[i][1]) -
			(point[0] - c[i][0]) * (c[(i + 1) % 3][1] - c[i][1]));
	}
	square /= 2.0;

	// if half of their sum is equal to element's geometrical area
	if (abs (square - area) < PRECIS)
	{
		// point is in the element
		return true;
	}
	return false;
}

bool NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::triangle_degenerate (int * n)
{
	double c[3][2];
	for (int k = 0; k < 3; k++)
		get_node_coordinates (n[k], c[k]);
	double area = (c[0][0] - c[2][0]) * (c[1][1] - c[2][1]) -
		(c[1][0] - c[2][0]) * (c[0][1] - c[2][1]);
	area = fabs (area) / 2.0;

	return area < AREA_DEGENERATE_TRIANGLE;
}

int NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::build_mesh (char * file_input, char * file_nodes_output, char * file_elements_output)
{
	printf ("\tgenerating points\n");
	pg = new Triangulation::Point_Generator ();
	pg->generate_points (file_input);
	printf ("\tmaking elements\n");
	make_elements ();
	numerate_functions ();
	printf ("\trenumerating mesh\n");
	renumerate ();
	printf ("\toutput\n");
	output (file_nodes_output, file_elements_output);
	return 1;
}

int NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::build_mesh_npg (char * file_input, char * file_nodes_output, char * file_elements_output)
{
	printf ("\tgenerating points\n");
	pg = new Triangulation::Point_Generator ();
	pg->read_points (file_input);
	printf ("\tmaking elements\n");
	make_elements ();
	numerate_functions ();
	printf ("\trenumerating mesh\n");
	renumerate ();
	printf ("\toutput\n");
	output (file_nodes_output, file_elements_output);
	return 1;
}

int NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::point_inside (double * coordinates)
{
	if (tot != NULL)
		return tot->tree_search (coordinates);
	else return -1;
}

bool NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::get_isoline_section (int k_element, double * q, double value, double * c1, double * c2)
{
	int amount = elements[k_element]->get_isoline_points (*this, value, q, c1, c2);
	if (amount == 2)
		return true;
	return false;
}

void NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay::output (char * file_nodes, char * file_triangles)
{
	FILE * log = fopen (file_nodes, "w");

	double coordinates[2];
	for (int i = 0; i < n_nodes; i++)
	{
		nodes[i].get_coordinates (coordinates);
		fprintf (log, "%.16lf %.16lf\n", coordinates[0], coordinates[1]);
	}

	fclose (log);

	log = fopen (file_triangles, "w");

	//int nodes[3];
	for (int i = 0; i < n_elements; i++)
	{
		int n_functions = elements[i]->get_amount_of_base_nodes ();
		for (int k = 0; k < n_functions; k++)
		{
			fprintf (log, "%i ", elements[i]->get_node (k));
		}
		fprintf (log, "%i ", elements[i]->get_area ());
		fprintf (log, "1");
		fprintf (log, "\n");
	}

	fclose (log);
}

