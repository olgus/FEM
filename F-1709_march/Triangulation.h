#pragma once
#include "Basic_Figure.h"
#include <queue>
#include <vector>
#include <time.h>
#include <random>
#include <algorithm>

using namespace Basic_Geometry;

namespace Triangulation
{
#define LVL_PRECIS 1e-1
	struct Dense_Area
	{
		int type; // type
		Basic_Figure * bf; // the figure
		double lvl; // lvl
	};

	struct Side_point
	{
		Point p;
		double r;
	};

	struct Side
	{
		Point p0, p1;
		std::vector<Side_point> side_points;
	};

	class Master_Element
	{
	protected:
		double x_step, y_step;
		int n, m; // amount of steps
		double ** lvls;
		Basic_Figure * bf; // the figure
	public:
		Master_Element ();
		Master_Element (Basic_Figure * BF);
		~Master_Element ();
		virtual void set_steps (double base_x_step, double base_y_step);
		virtual void set_lvls (int n_dense_areas, Dense_Area * dense_areas);
		virtual void get_length_modifiers (double * x_step, double * y_step);
		virtual double get_max_step ();
		virtual void generate_points (std::vector<Point> * gen_points, std::vector<Side> * sides);
	};

	struct Master_Element_Square : public Master_Element
	{
		Master_Element_Square ();
		Master_Element_Square (Basic_Figure * BF);
		~Master_Element_Square ();

		void set_steps (double base_x_step, double base_y_step) override;
		void set_lvls (int n_dense_areas, Dense_Area * dense_areas) override;
		void generate_points (std::vector<Point> * gen_points, std::vector<Side> * sides) override;
	};

	// master element factory
	inline Master_Element * make_master_element (Basic_Figure * bf)
	{
		switch (bf->get_type ())
		{
		default:
			return NULL;
		case BASIC_FIGURE_TYPE_RECTANGLE:
		case BASIC_FIGURE_TYPE_CONVEX_QUADRILATERAL:
		{
			return new Master_Element_Square (bf);
		}
		}
	}

	// generate a point
	void generate_a_point (double * x_bound, double * y_bound, double * point);

	class Point_Generator
	{
		// generates points 
		// duh

	private:
		// points
		int n_base_points;
		Point * base_points;
		double * base_steps;
		int n_dense_areas;
		Dense_Area * dense_areas;
		std::vector<Side> sides;
		double base_lvl;
		// TODO separate inside points and boundary points ?
		// save somehow outside structures
		// boundary points of basic figures? 
		//std::vector <std::pair <Basic_Geometry::Point, Basic_Geometry::Point>> sides;

		void read_data (char * file_input);
		void output (char * file_output);
		// generate on figure's unique sides
		//void generate_on_a_side (int k_bf, const Master_Element & me);
		void generate_on_sides ();
		// generate points by generating them in master elements
		void generate_on_basic_figures ();
		// check that sides are on one side
		bool one_side (const Side & s1, const Side & s2);
	public:
		Point_Generator ();
		Point_Generator (const Point_Generator & pg);
		~Point_Generator ();

		Basic_Figure ** basic_figures;
		int n_basic_figures;

		// read area details
		void generate_points (char * file_input, char * file_output);
		void generate_points (char * file_input);

		// if you don't want to generate again
		void read_points (char * file_input);

		std::vector<Point> gen_points;
	};
}