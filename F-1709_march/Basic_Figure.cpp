#include "Basic_Figure.h"

Basic_Geometry::Basic_Figure::Basic_Figure ()
{
	base_points_amount = 0;
	base_points = NULL;
}

Basic_Geometry::Basic_Figure::Basic_Figure (const Basic_Figure & bf)
{
	if (base_points != NULL)
		delete base_points;
	base_points_amount = bf.base_points_amount;
	base_points = new Point[base_points_amount];
	for (int i = 0; i < base_points_amount; i++)
		base_points[i] = bf.base_points[i];
}

Basic_Geometry::Basic_Figure::~Basic_Figure ()
{
	if (base_points != NULL)
		delete[] base_points;
}

void Basic_Geometry::Basic_Figure::set_point (int k_point, const Point & point)
{
	if (k_point < base_points_amount)
	{
		base_points[k_point] = point;
	}
}

void Basic_Geometry::Basic_Figure::get_point (int k_point, Point * figure_point)
{
	if (k_point < base_points_amount)
	{
		*figure_point = base_points[k_point];
	}
}

int Basic_Geometry::Basic_Figure::get_amount_of_base_points ()
{
	return base_points_amount;
}

void Basic_Geometry::Basic_Figure::set_material (int Material)
{
	material = Material;
}

void Basic_Geometry::Basic_Figure::get_point (const Point & master_point, Point * figure_point)
{
	// translation formulas
	// specific for each figure
}

void Basic_Geometry::Basic_Figure::triangulate (int * n, int (*nodes)[3])
{
}

Basic_Geometry::Point Basic_Geometry::Basic_Figure::get_point (int k)
{
	return base_points[k];
}

void Basic_Geometry::Basic_Figure::get_length_modifiers (double * x_step, double * y_step)
{
	// returns length modifiers
	// specific for each figure
}

void Basic_Geometry::Basic_Figure::sort_points ()
{
	// sorts points in comvenient order
	// specific for each figure
}

bool Basic_Geometry::Basic_Figure::point_inside (const Point & p)
{
	return false;
}

int Basic_Geometry::Basic_Figure::get_type ()
{
	return BASIC_FIGURE_TYPE_NO;
}

int Basic_Geometry::Basic_Figure::get_material ()
{
	return material;
}

int Basic_Geometry::Basic_Figure::get_amount_of_sides ()
{
	return 0;
}

void Basic_Geometry::Basic_Figure::get_side (int k, Point * p0, Point * p1)
{
}

bool Basic_Geometry::Basic_Figure::has_a_side (const Point & p0, const Point & p1)
{
	return false;
}

Basic_Geometry::Rectangle::Rectangle ()
{
	base_points_amount = 4;
	base_points = new Point[base_points_amount];
}

Basic_Geometry::Rectangle::~Rectangle ()
{

}

void Basic_Geometry::Rectangle::triangulate (int * n, int (*nodes)[3])
{
	// two triangles
	*n = 2;
	// local numbers of nodes
	for (int i = 0; i < 3; i++)
	{
		nodes[0][i] = i;
		nodes[1][i] = i + 1;
	}
}

void Basic_Geometry::Rectangle::get_point (const Point & master_point, Point * figure_point)
{
	figure_point->x = master_point.x * (base_points[1].x - base_points[0].x) + base_points[0].x;
	figure_point->y = master_point.y * (base_points[2].y - base_points[0].y) + base_points[0].y;
}

void Basic_Geometry::Rectangle::get_length_modifiers (double * x_step, double * y_step)
{
	// return inverted sides lengths
	*x_step = 1.0 / (base_points[1].x - base_points[0].x); 
	*y_step = 1.0 / (base_points[2].y - base_points[0].y);
}

void Basic_Geometry::Rectangle::sort_points ()
{
	// points order from min x to max x, from min y to max y 
	Point point = { 1e+10, 1e+10 };
	int pos = -1;
	// find min x min y
	for (int i = 0; i < 4; i++)
	{
		if (base_points[i] < point)
		{
			pos = i;
			point = base_points[i];
		}
	}
	point = base_points[0];
	base_points[0] = base_points[pos];
	base_points[pos] = point;
	
	// find max x max y
	point = { -1e+10, -1e+10 };
	for (int i = 0; i < 4; i++)
	{
		if (base_points[i] > point)
		{
			pos = i;
			point = base_points[i];
		}
	}
	point = base_points[3];
	base_points[3] = base_points[pos];
	base_points[pos] = point;
	// sort 1 and 2 point
	if (base_points[1].y > base_points[2].y)
	{
		point = base_points[1];
		base_points[1] = base_points[2];
		base_points[2] = point;
	}
}

bool Basic_Geometry::Rectangle::point_inside (const Point & p)
{
	if (((base_points[0].x - PRECIS < p.x) && (p.x < base_points[3].x + PRECIS))
		&& ((base_points[0].y - PRECIS < p.y) && (p.y < base_points[3].y + PRECIS)))
		return true;
	return false;
}

int Basic_Geometry::Rectangle::get_type ()
{
	return BASIC_FIGURE_TYPE_RECTANGLE;
}

int Basic_Geometry::Rectangle::get_amount_of_sides ()
{
	return 4;
}

void Basic_Geometry::Rectangle::get_side (int k, Point * p0, Point * p1)
{
	switch (k)
	{
	case 0:
	{
		*p0 = base_points[0];
		*p1 = base_points[1];
		break;
	}
	case 1:
	{
		*p0 = base_points[2];
		*p1 = base_points[3];
		break;
	}
	case 2:
	{
		*p0 = base_points[0];
		*p1 = base_points[2];
		break;
	}
	case 3:
	{
		*p0 = base_points[1];
		*p1 = base_points[3];
		break;
	}
	}


}

bool Basic_Geometry::Rectangle::has_a_side (const Point & p0, const Point & p1)
{
	if (p0 == base_points[0] && p1 == base_points[1])
		return true;
	if (p0 == base_points[2] && p1 == base_points[3])
		return true;
	if (p0 == base_points[0] && p1 == base_points[2])
		return true;
	if (p0 == base_points[1] && p1 == base_points[3])
		return true;
	return false;
}

Basic_Geometry::Convex_quadrilateral::Convex_quadrilateral ()
{
	base_points_amount = 4;
	base_points = new Point[base_points_amount];
}

Basic_Geometry::Convex_quadrilateral::~Convex_quadrilateral ()
{
}

void Basic_Geometry::Convex_quadrilateral::triangulate (int * n, int (*nodes)[3])
{
	// two triangles
	*n = 2;
	// local numbers of nodes
	for (int i = 0; i < 3; i++)
	{
		nodes[0][i] = i;
		nodes[1][i] = i + 1;
	}
}

void Basic_Geometry::Convex_quadrilateral::get_point (const Point & master_point, Point * figure_point)
{
	figure_point->x = (1.0 - master_point.x) * (1.0 - master_point.y) * base_points[0].x +
		master_point.x * (1.0 - master_point.y) * base_points[1].x +
		(1.0 - master_point.x) * master_point.y * base_points[2].x +
		master_point.x * master_point.y * base_points[3].x;
	figure_point->y = (1.0 - master_point.x) * (1.0 - master_point.y) * base_points[0].y +
		master_point.x * (1.0 - master_point.y) * base_points[1].y +
		(1.0 - master_point.x) * master_point.y * base_points[2].y +
		master_point.x * master_point.y * base_points[3].y;
}

void Basic_Geometry::Convex_quadrilateral::get_length_modifiers (double * x_step, double * y_step)
{
	// return inverted sides lengths
	double l[4];
	l[0] = sqrt (pow (base_points[0].x - base_points[1].x, 2.0) + pow (base_points[0].y - base_points[1].y, 2.0));
	l[1] = sqrt (pow (base_points[2].x - base_points[3].x, 2.0) + pow (base_points[2].y - base_points[3].y, 2.0));

	l[2] = sqrt (pow (base_points[0].x - base_points[2].x, 2.0) + pow (base_points[0].y - base_points[2].y, 2.0));
	l[3] = sqrt (pow (base_points[1].x - base_points[3].x, 2.0) + pow (base_points[1].y - base_points[3].y, 2.0));
	*x_step = l[0] > l[1] ? 1.0 / l[0] : 1.0 / l[1];
	*y_step = l[2] > l[3] ? 1.0 / l[2] : 1.0 / l[3];
}

void Basic_Geometry::Convex_quadrilateral::sort_points ()
{
	// sort by y increasing order
	std::sort (base_points, base_points + base_points_amount, 
		[](Point a, Point b) { return a.y < b.y; });
	Point p;
	if (base_points[0].x > base_points[1].x)
	{
		p = base_points[1];
		base_points[1] = base_points[0];
		base_points[0] = p;
	}
	if (base_points[2].x > base_points[3].x)
	{
		p = base_points[3];
		base_points[3] = base_points[3];
		base_points[2] = p;
	}
	
}

bool Basic_Geometry::Convex_quadrilateral::point_inside (const Point & p)
{
	// calc sq area of two triangles
	double cqsa = fabs ((base_points[0].x - base_points[2].x) * (base_points[1].y - base_points[2].y) - 
		(base_points[1].x - base_points[2].x) * (base_points[0].y - base_points[2].y)) / 2.0 +
		fabs ((base_points[3].x - base_points[2].x) * (base_points[1].y - base_points[2].y) -
		(base_points[1].x - base_points[2].x) * (base_points[3].y - base_points[2].y)) / 2.0;

	// cacl sq area of 4 triangles
	double psa = 0.0;
	for (int i = 0; i < 4; i++)
	{
		psa += fabs ((base_points[(i / 2) * 3].x - p.x) * (base_points[i % 2 + 1].y - p.y) -
			(base_points[i % 2 + 1].x - p.x) * (base_points[(i / 2) * 3].y - p.y)) / 2.0;
	}
	if (fabs (psa - cqsa) < PRECIS)
		return true;
	return false;
}

int Basic_Geometry::Convex_quadrilateral::get_type ()
{
	return BASIC_FIGURE_TYPE_CONVEX_QUADRILATERAL;
}

int Basic_Geometry::Convex_quadrilateral::get_amount_of_sides ()
{
	return 4;
}

void Basic_Geometry::Convex_quadrilateral::get_side (int k, Point * p0, Point * p1)
{
	switch (k)
	{
	case 0:
	{
		*p0 = base_points[0];
		*p1 = base_points[1];
		break;
	}
	case 1:
	{
		*p0 = base_points[2];
		*p1 = base_points[3];
		break;
	}
	case 2:
	{
		*p0 = base_points[0];
		*p1 = base_points[2];
		break;
	}
	case 3:
	{
		*p0 = base_points[1];
		*p1 = base_points[3];
		break;
	}
	}
}

bool Basic_Geometry::Convex_quadrilateral::has_a_side (const Point & p0, const Point & p1)
{
	if (p0 == base_points[0] && p1 == base_points[1])
		return true;
	if (p0 == base_points[2] && p1 == base_points[3])
		return true;
	if (p0 == base_points[0] && p1 == base_points[2])
		return true;
	if (p0 == base_points[1] && p1 == base_points[3])
		return true;
	return false;
}
