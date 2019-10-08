#include "Triangulation.h"

namespace Triangulation
{
	void Point_Generator::read_data (char * file_input)
	{
		FILE * file = fopen (file_input, "r");
		// amount of mesh base points
		fscanf (file, "%i", &n_base_points);
		base_points = new Point[n_base_points];
		// level
		fscanf (file, "%lf", &base_lvl);
		// mesh base points
		for (int i = 0; i < n_base_points; i++)
		{
			fscanf (file, "%lf %lf", &base_points[i].x, &base_points[i].y);
		}
		// amount of basic figures
		fscanf (file, "%i", &n_basic_figures);
		basic_figures = new Basic_Figure *[n_basic_figures];
		base_steps = new double[n_basic_figures];
		{
			int type;
			int material;
			for (int i = 0; i < n_basic_figures; i++)
			{
				// figure type
				fscanf (file, "%i", &type);
				// material
				fscanf (file, "%i", &material);
				// base step length
				fscanf (file, "%lf", &base_steps[i]);
				base_steps[i] *= pow (2.0, -base_lvl);
				// respective number of points
				basic_figures[i] = make_basic_figure (type);
				basic_figures[i]->set_material (material);
				int point;
				for (int p = 0, p_end = basic_figures[i]->get_amount_of_base_points (); p < p_end; p++)
				{
					fscanf (file, "%i", &point);
					basic_figures[i]->set_point (p, base_points[point]);
				}
				basic_figures[i]->sort_points ();
			}
		}

		// amount of areas with deeper density
		fscanf (file, "%i", &n_dense_areas);
		// areas with deeper density
		if (n_dense_areas > 0)
		{
			dense_areas = new Dense_Area[n_dense_areas];
			{
				Point point;
				for (int i = 0; i < n_dense_areas; i++)
				{
					fscanf (file, "%i %lf", &dense_areas[i].type, &dense_areas[i].lvl);
					dense_areas[i].bf = make_basic_figure (dense_areas[i].type);
					for (int p = 0, p_end = dense_areas[i].bf->get_amount_of_base_points (); p < p_end; p++)
					{
						fscanf (file, "%lf %lf", &point.x, &point.y);
						dense_areas[i].bf->set_point (p, point);
					}
					dense_areas[i].bf->sort_points ();
				}
			}
		}
	}

	void Point_Generator::output (char * file_output)
	{
		FILE * file = fopen (file_output, "w");
		for (int i = 0, i_end = (int)gen_points.size (); i < i_end; i++)
		{
			fprintf (file, "%.15lf %.15lf\n", gen_points[i].x, gen_points[i].y);
		}
		fclose (file);
	}

	Point_Generator::Point_Generator ()
	{
		n_base_points = 0;
		n_basic_figures = 0;
		n_dense_areas = 0;

		base_points = NULL;
		basic_figures = NULL;
		dense_areas = NULL;
		base_steps = NULL;
		gen_points.clear ();
	}

	Point_Generator::Point_Generator (const Point_Generator & pg)
	{
	}

	Point_Generator::~Point_Generator ()
	{
		if (base_points != NULL)
			delete[] base_points;
		if (basic_figures != NULL)
		{
			for (int i = 0;i < n_basic_figures; i++)
			{
				delete basic_figures[i];
			}
			delete[] basic_figures;
		}
		if (dense_areas != NULL)
			delete[] dense_areas;
		if (base_steps != NULL)
			delete[] base_steps;
	}

	void Point_Generator::generate_points (char * file_input, char * file_output)
	{
		read_data (file_input);
		generate_on_basic_figures ();
		output (file_output);
	}

	void Point_Generator::generate_points (char * file_input)
	{
		read_data (file_input);
		generate_on_basic_figures ();
	}

	void Point_Generator::read_points (char * file_input, char * file_points_input)
	{
		FILE * file = fopen (file_input, "r");
		// amount of mesh base points
		fscanf (file, "%i", &n_base_points);
		base_points = new Point[n_base_points];
		// level
		fscanf (file, "%lf", &base_lvl);
		// mesh base points
		for (int i = 0; i < n_base_points; i++)
		{
			fscanf (file, "%lf %lf", &base_points[i].x, &base_points[i].y);
		}
		// amount of basic figures
		fscanf (file, "%i", &n_basic_figures);
		basic_figures = new Basic_Figure *[n_basic_figures];
		base_steps = new double[n_basic_figures];
		{
			int type;
			int material;
			for (int i = 0; i < n_basic_figures; i++)
			{
				// figure type
				fscanf (file, "%i", &type);
				// material
				fscanf (file, "%i", &material);
				// base step length
				fscanf (file, "%lf", &base_steps[i]);
				base_steps[i] *= pow (2.0, -base_lvl);
				// respective number of points
				basic_figures[i] = make_basic_figure (type);
				basic_figures[i]->set_material (material);
				int point;
				for (int p = 0, p_end = basic_figures[i]->get_amount_of_base_points (); p < p_end; p++)
				{
					fscanf (file, "%i", &point);
					basic_figures[i]->set_point (p, base_points[point]);
				}
				basic_figures[i]->sort_points ();
			}
		}
		fclose (file);

		file = fopen (file_points_input, "r");
		// read points
		{
			double x, y;
			while (!feof (file))
			{
				fscanf (file, "%lf %lf", &x, &y);
				Point p = { x, y };
				gen_points.push_back (p);
			}
			fclose (file);
			gen_points.pop_back ();
		}
	}
		
	void Point_Generator::generate_on_sides ()
	{
		std::vector <std::vector <int>> seq_sides;
		int same_side;
		// sort sides into groups that make one side
		for (int i = 0, i_end = (int)sides.size (); i < i_end; i++)
		{
			same_side = -1;
			for (int j = 0, j_end = (int)seq_sides.size (); j < j_end && same_side < 0; j++)
			{
				// check if current sides is one the same side as all others in a group
				if (one_side (sides[seq_sides[j][0]], sides[i]))
				{
					same_side = j;
				}
			} 
			if (same_side != -1)
			{
				seq_sides[same_side].push_back (i);
			}
			else
			{
				std::vector <int> s = {i};
				seq_sides.push_back (s);
			}
		}
		// go by groups
		for (int i = 0, i_end = (int)seq_sides.size(); i < i_end; i++)
		{
			struct Point_on_a_side
			{
				Side_point sp;
				int side;
			};
			std::vector <Point_on_a_side> points_on_a_side;
			for (int j = 0, j_end = (int)seq_sides[i].size (); j < j_end; j++)
			{
				for (int k = 0, k_end = (int)sides[seq_sides[i][j]].side_points.size (); k < k_end; k++)
				{
					points_on_a_side.push_back ({ sides[seq_sides[i][j]].side_points[k], seq_sides[i][j] });
				}
			}
			// sort points in each group
			std::sort (points_on_a_side.begin (), points_on_a_side.end (),
				[](Point_on_a_side a, Point_on_a_side b) {return a.sp.p < b.sp.p; });
			// go by points
			bool add, prev;
			prev = true;
			if (points_on_a_side.size () > 2)
			{
				gen_points.push_back (points_on_a_side[0].sp.p);
				for (int j = 1, j_end = (int)points_on_a_side.size () - 1; j < j_end; j++)
				{
					add = false;
					if ((points_on_a_side[j - 1].side == points_on_a_side[j + 1].side))
					{
						if ((points_on_a_side[j].side != points_on_a_side[j + 1].side) &&
							(points_on_a_side[j].side != points_on_a_side[j - 1].side))
						{
							if ((points_on_a_side[j].sp.r < points_on_a_side[j + 1].sp.r - PRECIS))
								add = true;
							if (!prev)
								add = true;
						}
						else
						{
							add = true;
						}
					}
					else
					{
						if (prev)
						{
							double dist1 = points_on_a_side[j].sp.p.distance (points_on_a_side[j + 1].sp.p);
							double dist2 = points_on_a_side[j].sp.p.distance (points_on_a_side[j - 1].sp.p);
							double dist = dist1 > dist2 ? dist2 : dist1;
							double r = points_on_a_side[j].sp.r > points_on_a_side[j + 1].sp.r ?
								(points_on_a_side[j + 1].sp.r > points_on_a_side[j - 1].sp.r ? points_on_a_side[j - 1].sp.r : points_on_a_side[j + 1].sp.r) :
								(points_on_a_side[j].sp.r > points_on_a_side[j - 1].sp.r ? points_on_a_side[j - 1].sp.r : points_on_a_side[j].sp.r);
							if (dist > r / 2.0)
								add = true;
						}
						else
							add = true;
					}

					if (add)
					{
						gen_points.push_back (points_on_a_side[j].sp.p);
						prev = true;
					}
					else
						prev = false;
				}
				gen_points.push_back (points_on_a_side[(int)points_on_a_side.size () - 1].sp.p);
			}
			else
			{
				for (int j = 0, j_end = (int)points_on_a_side.size (); j < j_end; j++)
					gen_points.push_back (points_on_a_side[j].sp.p);
			}
		}
	}

	void Point_Generator::generate_on_basic_figures ()
	{
		double xlm, ylm;
		// add all base points
		//for (int i = 0; i < n_base_points; i++)
		//{
		//	gen_points.push_back (base_points[i]);
		//}
		Master_Element ** mes;
		mes = new Master_Element*[n_basic_figures];
		for (int k = 0; k < n_basic_figures; k++)
		{
			// create a master element for bf
			mes[k] = make_master_element (basic_figures[k]);
			// define base steps
			mes[k]->get_length_modifiers (&xlm, &ylm);
			mes[k]->set_steps (base_steps[k] * xlm, base_steps[k] * ylm);
			// set lvls for said master element
			mes[k]->set_lvls (n_dense_areas, dense_areas);
		}
		//// sort basic figures by max(lvl)
		//{
		//	std::sort (mes, mes + n_basic_figures,
		//		[](Master_Element* a, Master_Element* b) {return a->get_max_step () < b->get_max_step (); });
		//}
		// go by basic figures
		for (int k = 0; k < n_basic_figures; k++)
		{
			mes[k]->generate_points (&gen_points, &sides);
		}
		generate_on_sides ();
		std::reverse (gen_points.begin (), gen_points.end ());

		for (int k = 0; k < n_basic_figures; k++)
			delete mes[k];
		delete[] mes;
	}

	bool Point_Generator::one_side (const Side & s1, const Side & s2)
	{
		double sq1 = ((s1.p0.x - s2.p0.x) * (s1.p1.y - s2.p0.y) - 
			(s1.p1.x - s2.p0.x) * (s1.p0.y - s2.p0.y)) / 2.0;
		double sq2 = ((s1.p0.x - s2.p1.x) * (s1.p1.y - s2.p1.y) -
			(s1.p1.x - s2.p1.x) * (s1.p0.y - s2.p1.y)) / 2.0;
		if ((fabs (sq1) < PRECIS) && fabs (sq2) < PRECIS)
			return true;
		return false;
	}

	Master_Element::Master_Element ()
	{
		lvls = NULL;
		bf = NULL;
	}

	Master_Element::Master_Element (Basic_Figure * BF)
	{
		lvls = NULL;
		bf = BF;
	}

	Master_Element::~Master_Element ()
	{
		if (lvls != NULL)
		{
			for (int i = 0; i < n; i++)
				delete[] lvls[i];
			delete[] lvls;
		}
	}

	void Master_Element::set_steps (double base_x_step, double base_y_step)
	{
		
	}

	void Master_Element::set_lvls (int n_dense_areas, Dense_Area * dense_areas)
	{
	
	}

	void Master_Element::get_length_modifiers (double * x_step, double * y_step)
	{
		if (bf != NULL)
			bf->get_length_modifiers (x_step, y_step);
	}

	double Master_Element::get_max_step ()
	{
		return x_step > y_step ? x_step : y_step;
	}

	void Master_Element::generate_points (std::vector<Point>* gen_points, std::vector<Side>* sides)
	{
	
	}

	Master_Element_Square::Master_Element_Square ()
	{
	}

	Master_Element_Square::Master_Element_Square (Basic_Figure * BF) : Master_Element (BF)
	{
	}
		
	Master_Element_Square::~Master_Element_Square ()
	{
	}

	void Master_Element_Square::set_steps (double base_x_step, double base_y_step)
	{
		if (base_x_step > 1.0 - PRECIS)
			m = 1;
		else
		{
			m = (int)(1.0 / base_x_step);
			double leftover = 1.0 - m * base_x_step;
			if (leftover > 0.5 * base_x_step)
				m++;
		}
		x_step = 1.0 / (double)m;

		if (base_y_step > 1.0 - PRECIS)
			n = 1;
		else
		{
			n = (int)(1.0 / base_y_step);
			double leftover = 1.0 - n * base_y_step;
			if (leftover > 0.5 * base_y_step)
				n++;
		}
		y_step = 1.0 / (double)n;
	}

	void Master_Element_Square::set_lvls (int n_dense_areas, Dense_Area * dense_areas)
	{
		lvls = new double *[n];
		for (int i = 0; i < n; i++)
			lvls[i] = new double[m];

		// init lvls
		double lvl;
		int count;
		double x[3], y[3];
		Point mp, fp;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				// base line is 0
				lvl = 0.0;
				// check if rectangle falls into the area
				if (n_dense_areas > 0)
				{
					x[0] = x_step * j;
					x[1] = (x_step * (2 * j + 1)) / 2.0;
					x[2] = x_step * (j + 1);
					y[0] = y_step * i;
					y[1] = (y_step * (2 * i + 1)) / 2.0;
					y[2] = y_step * (i + 1);

					for (int k = 0; k < n_dense_areas; k++)
					{
						count = 0;
						for (int p = 0; p < 9 && count < 5; p++)
						{
							mp.x = x[p % 3];
							mp.y = y[p / 3];
							bf->get_point (mp, &fp);
							if (dense_areas[k].bf->point_inside (fp))
							{
								// set the lvl as the lvl of the area
								count++;
							}
							if (count > 4)
							{
								if (dense_areas[k].lvl > lvl)
									lvl = dense_areas[k].lvl;
							}
						}
					}
				}
				lvls[i][j] = lvl;
			}
		}

		// smooth
		// make a queue
		std::queue<std::pair<int, int>> Q;
		// start the queue as full of areas with max lvl > 0
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				if (lvls[i][j] > PRECIS)
					Q.push (std::pair <int, int> (i, j));
			}
		}
		std::pair<int, int> p;
		// while queue is not empty
		while (!Q.empty ())
		{
			p = Q.front ();
			Q.pop ();
			// check areas on the sides 
			// if the difference is > 1
			// change the lvl of the area to +1
			// push it into the queue

			if (p.first > 0)
			{
				if (lvls[p.first][p.second] - lvls[p.first - 1][p.second] > 1.0)
				{
					lvls[p.first - 1][p.second] = lvls[p.first][p.second] - 1;
					Q.push (std::pair <int, int> (p.first - 1, p.second));
				}

				if (p.second > 0)
				{
					if (lvls[p.first][p.second] - lvls[p.first - 1][p.second - 1] > 1.0)
					{
						lvls[p.first - 1][p.second - 1] = lvls[p.first][p.second] - 1;
						Q.push (std::pair <int, int> (p.first - 1, p.second - 1));
					}
				}
				if (p.second < m - 1)
				{
					if (lvls[p.first][p.second] - lvls[p.first - 1][p.second + 1] > 1.0)
					{
						lvls[p.first - 1][p.second + 1] = lvls[p.first][p.second] - 1;
						Q.push (std::pair <int, int> (p.first - 1, p.second + 1));
					}
				}
			}
			if (p.first < n - 1)
			{
				if (lvls[p.first][p.second] - lvls[p.first + 1][p.second] > 1.0)
				{
					lvls[p.first + 1][p.second] = lvls[p.first][p.second] - 1;
					Q.push (std::pair <int, int> (p.first + 1, p.second));
				}

				if (p.second > 0)
				{
					if (lvls[p.first][p.second] - lvls[p.first + 1][p.second - 1] > 1.0)
					{
						lvls[p.first + 1][p.second - 1] = lvls[p.first][p.second] - 1;
						Q.push (std::pair <int, int> (p.first + 1, p.second - 1));
					}
				}
				if (p.second < m - 1)
				{
					if (lvls[p.first][p.second] - lvls[p.first + 1][p.second + 1] > 1.0)
					{
						lvls[p.first + 1][p.second + 1] = lvls[p.first][p.second] - 1;
						Q.push (std::pair <int, int> (p.first + 1, p.second + 1));
					}
				}
			}
			if (p.second > 0)
			{
				if (lvls[p.first][p.second] - lvls[p.first][p.second - 1] > 1.0)
				{
					lvls[p.first][p.second - 1] = lvls[p.first][p.second] - 1;
					Q.push (std::pair <int, int> (p.first, p.second - 1));
				}
			}
			if (p.second < m - 1)
			{
				if (lvls[p.first][p.second] - lvls[p.first][p.second + 1] > 1.0)
				{
					lvls[p.first][p.second + 1] = lvls[p.first][p.second] - 1;
					Q.push (std::pair <int, int> (p.first, p.second + 1));
				}
			}
		}
	}
	
	void Master_Element_Square::generate_points (std::vector<Point>* gen_points, std::vector<Side>* sides)
	{
		// add sides
		{
			Point p0, p1;
			for (int i = 0, i_end = bf->get_amount_of_sides (); i < i_end; i++)
			{
				bf->get_side (i, &p0, &p1);
				sides->push_back ({ p0, p1 });
			}

		}
		// add corner points into gen_points
		for (int k = 0, k_end = bf->get_amount_of_base_points (); k < k_end; k++)
			gen_points->push_back (bf->get_point (k));
		// go by m-e sections
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				// if section is lvl 0, set the point randomly in center area
				if (lvls[i][j] < LVL_PRECIS)
				{
					double x_bound[2];
					double y_bound[2];
					double point[2];
					{
						switch (bf->get_type ())
						{
						case BASIC_FIGURE_TYPE_RECTANGLE:
						{
							Point mp = { x_step * j,  y_step * i };
							Point fp0, fpN;
							bf->get_point (mp, &fp0);
							mp = { x_step * (j + 1), y_step * (i + 1) };
							bf->get_point (mp, &fpN);

							x_bound[0] = fp0.x + (fpN.x - fp0.x) * 0.2;
							y_bound[0] = fp0.y + (fpN.y - fp0.y) * 0.2;
							x_bound[1] = fp0.x + (fpN.x - fp0.x) * 0.8;
							y_bound[1] = fp0.y + (fpN.y - fp0.y) * 0.8;

							// keep the side
							{
								double xlm, ylm;
								bf->get_length_modifiers (&xlm, &ylm);
								if (i == 0 && j != m - 1)
								{
									sides->operator[]((int)sides->size () - 4).side_points.push_back ({ { fpN.x, fp0.y }, x_step / xlm });
								}
								if (i == n - 1 && j != m - 1)
								{
									sides->operator[]((int)sides->size () - 3).side_points.push_back ({ { fpN.x, fpN.y }, x_step / xlm });
								}
								if (j == 0 && i != n - 1)
								{
									sides->operator[]((int)sides->size () - 2).side_points.push_back ({ { fp0.x, fpN.y }, y_step / ylm });
								}
								if (j == m - 1 && i != n - 1)
								{
									sides->operator[]((int)sides->size () - 1).side_points.push_back ({ { fpN.x, fpN.y }, y_step / ylm });
								}
							}
							break;
						}
						}

					}
					Triangulation::generate_a_point (x_bound, y_bound, point);
					gen_points->push_back ({ point[0], point[1] });

				}
				else
				{
					// if section is bigger lvl
					// define steps for it
					double sec_x_step = pow (2.0, -lvls[i][j]) * x_step;
					double sec_y_step = pow (2.0, -lvls[i][j]) * y_step;
					// define amounts
					int sec_n, sec_m;
					sec_m = (int)(x_step / sec_x_step);
					sec_n = (int)(y_step / sec_y_step);
					if (x_step - sec_x_step * sec_m > 0.5 * sec_x_step)
					{
						sec_m++;
					}
					sec_x_step = x_step / (double)(sec_m);
					if (y_step - sec_y_step * sec_n > 0.5 * sec_y_step)
					{
						sec_n++;
					}
					sec_y_step = y_step / (double)(sec_n);

					for (int sec_i = 0; sec_i < sec_n; sec_i++)
					{
						for (int sec_j = 0; sec_j < sec_m; sec_j++)
						{
							{
								double x_bound[2];
								double y_bound[2];
								double point[2];
								{
									switch (bf->get_type ())
									{
									case BASIC_FIGURE_TYPE_RECTANGLE:
									{
										Point mp = { x_step * j + sec_j * sec_x_step,  y_step * i + sec_i * sec_y_step };
										Point fp0, fpN;
										bf->get_point (mp, &fp0);
										mp = { x_step * j + (sec_j + 1) * sec_x_step,  y_step * i + (sec_i + 1) * sec_y_step };
										bf->get_point (mp, &fpN);

										x_bound[0] = fp0.x + (fpN.x - fp0.x) * 0.2;
										y_bound[0] = fp0.y + (fpN.y - fp0.y) * 0.2;
										x_bound[1] = fp0.x + (fpN.x - fp0.x) * 0.8;
										y_bound[1] = fp0.y + (fpN.y - fp0.y) * 0.8;

										// keep the side
										{
											double xlm, ylm;
											bf->get_length_modifiers (&xlm, &ylm);
											if (i == 0 && sec_i == 0 /*&& sec_j != sec_m - 1*/)
											{
												sides->operator[]((int)sides->size () - 4).side_points.push_back ({ { fpN.x, fp0.y }, sec_x_step / xlm });
											}
											if (i == n - 1 && sec_i == sec_n - 1 && j != m - 1)
											{
												sides->operator[]((int)sides->size () - 3).side_points.push_back ({ { fpN.x, fpN.y }, sec_x_step / xlm });
											}
											if (j == 0 && sec_j == 0 && i != n - 1)
											{
												sides->operator[]((int)sides->size () - 2).side_points.push_back ({ { fp0.x, fpN.y }, sec_y_step / ylm });
											}
											if (j == m - 1 && sec_j == sec_m - 1 /*&& sec_i != sec_n - 1*/)
											{
												sides->operator[]((int)sides->size () - 1).side_points.push_back ({ { fpN.x, fpN.y }, sec_y_step / ylm });
											}
										}

										break;
									}
									}

								}
								Triangulation::generate_a_point (x_bound, y_bound, point);
								gen_points->push_back ({ point[0], point[1] });
							}
						}
					}
				}
			}
		}
	}

	void Triangulation::generate_a_point (double * x_bound, double * y_bound, double * point)
	{
		static std::random_device rd;
		std::mt19937 rng (rd ());
		std::uniform_real_distribution<double> gen_x (x_bound[0], x_bound[1]); // uniform, unbiased
		point[0] = gen_x (rng);
		std::uniform_real_distribution<double> gen_y (y_bound[0], y_bound[1]); // uniform, unbiased
		point[1] = gen_y (rng);
		//point[0] = (x_bound[1] + x_bound[0]) / 2.0;
		//point[1] = (y_bound[1] + y_bound[0]) / 2.0; 
	}
}