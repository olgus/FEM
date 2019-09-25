#pragma once
#include <stdio.h>
#include <math.h>
#include <algorithm>
// basic figures that can be translated onto master element
// things to add
// different step length by x and y
// translation formulas
// 

namespace Basic_Geometry
{
#define PRECIS 1e-10
	// types of areas
	enum Figure_Type {
		BASIC_FIGURE_TYPE_NO, BASIC_FIGURE_TYPE_RECTANGLE, BASIC_FIGURE_TYPE_CONVEX_QUADRILATERAL
	};
	
	struct Point
	{
		double x, y;

		Point & operator= (const Point & p)
		{
			x = p.x;
			y = p.y;
			return *this;
		}

		bool operator<(const Point & p) const
		{
			if (fabs (x - p.x) < PRECIS)
				return (y < p.y);
			else
				return (x < p.x);

			//return ((x < p.x) || (fabs (x - p.x) < PRECIS) && (y < p.y));
		}

		bool operator>(const Point & p) const
		{
			return ((x > p.x) || (fabs (x - p.x) < PRECIS) && (y > p.y));
		}

		bool operator==(const Point & p) const
		{
			return ((fabs (x - p.x) < PRECIS) && (fabs (y - p.y) < PRECIS));
		}

		double distance (const Point & p) const
		{
			return (sqrt (pow (x - p.x, 2.0) + pow (y - p.y, 2.0)));
		}
	};

	class Basic_Figure
	{
	protected:
		int base_points_amount;
		int material;
		Point * base_points;
	public:	
		Basic_Figure ();
		Basic_Figure (const Basic_Figure & bf);
		~Basic_Figure ();
		// set point of a figure
		void set_point (int k_point, const Point & point);
		// translation
		void get_point (int k_point, Point * figure_point);
		// get amount of base points
		int get_amount_of_base_points ();

		virtual void set_material (int Material);
		virtual void get_point (const Point & master_point, Point * figure_point);
		virtual void triangulate (int * n, int (*nodes)[3]);
		virtual Point get_point (int k);
		// step lenght
		virtual void get_length_modifiers (double * x_step, double * y_step);
		// sort the points for convenience
		virtual void sort_points ();
		// check if the point is inside the figure
		virtual bool point_inside (const Point & p);
		virtual int get_type ();
		virtual int get_material ();
		virtual int get_amount_of_sides ();
		virtual void get_side (int k, Point * p0, Point * p1);
		virtual bool has_a_side (const Point & p0, const Point & p1);
	};
	// rectangle
	class Rectangle : public Basic_Figure
	{
	private:
	public:
		Rectangle ();
		~Rectangle ();

		void triangulate (int * n, int (*nodes)[3]) override;
		void get_point (const Point & master_point, Point * figure_point) override;
		// step lenght
		void get_length_modifiers (double * x_step, double * y_step) override;
		// sort the points for convenience
		void sort_points () override;
		// check if the point is inside the figure
		bool point_inside (const Point & p)override;
		int get_type () override;
		int get_amount_of_sides () override;
		void get_side (int k, Point * p0, Point * p1) override;
		bool has_a_side (const Point & p0, const Point & p1) override;
	};
	// convex quadrilateral
	class Convex_quadrilateral : public Basic_Figure
	{
	private:
	public:
		Convex_quadrilateral ();
		~Convex_quadrilateral ();

		void triangulate (int * n, int (*nodes)[3]) override;
		void get_point (const Point & master_point, Point * figure_point) override;
		// step lenght
		void get_length_modifiers (double * x_step, double * y_step) override;
		// sort the points for convenience
		void sort_points () override;
		// check if the point is inside the figure
		bool point_inside (const Point & p)override;
		int get_type () override;
		int get_amount_of_sides () override;
		void get_side (int k, Point * p0, Point * p1) override;
		bool has_a_side (const Point & p0, const Point & p1) override;
	};
	// triangle

	// figure factory
	inline Basic_Figure * make_basic_figure (int figure_type)
	{
		switch (figure_type)
		{
		default:
			return NULL;
		case BASIC_FIGURE_TYPE_RECTANGLE:
			return new Rectangle ();
		case BASIC_FIGURE_TYPE_CONVEX_QUADRILATERAL:
			return new Convex_quadrilateral ();
		}
	}
}