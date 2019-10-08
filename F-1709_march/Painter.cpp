#include "Painter.h"

Painter * Painter::currentPainter = NULL;

void Painter::initialize (int Argc, char * Argv[])
{
	argc = Argc;
	argv = Argv;
}

void Painter::resize ()
{
	double L[2];
	L[0] = (cN[0] - c0[0]) * scale_x;
	L[1] = (cN[1] - c0[1]) * scale_y;

	double coef;
	L[0] < L[1] ? coef = L[0] / L[1] : coef = L[1] / L[0];

	pic_height = pic_width = width;

	if (L[0] < L[1])
		pic_width = width = (int)(height * coef);
	else
		pic_height = height = (int)(width * coef);

	width += 2 * Frame;
	height += 2 * Frame;

	if (legend_blocks > 0)
	{
		// there is legend block
		width += W_legend + Frame;

		picture.c0[0] = 2.0 * (double)(Frame) / (double)width - 1.0;
		picture.c0[1] = 2.0 * (double)(Frame) / (double)height - 1.0;

		picture.cN[0] = 2.0 * (double)(Frame + pic_width) / (double)width - 1.0;
		picture.cN[1] = 2.0 * (double)(Frame + pic_height) / (double)height - 1.0;

		legend.c0[0] = 2.0 * (double)(2.0 * Frame + pic_width) / (double)width - 1.0;
		legend.c0[1] = 2.0 * (double)(Frame) / (double)height - 1.0;

		legend.cN[0] = 1.0 - (double)(2.0 * Frame) / width;
		legend.cN[1] = 2.0 * (double)(Frame + pic_height) / (double)height - 1.0;
	}
	else
	{
		// no legend 
		picture.c0[0] = 2.0 * (double)(Frame) / (double)width - 1.0;
		picture.c0[1] = 2.0 * (double)(Frame) / (double)height - 1.0;

		picture.cN[0] = 2.0 * (double)(Frame + pic_width) / (double)width - 1.0;
		picture.cN[1] = 2.0 * (double)(Frame + pic_height) / (double)height - 1.0;
	}
}

void Painter::draw_blocks ()
{
	glClear (GL_COLOR_BUFFER_BIT);

	glLoadIdentity ();
	glMatrixMode (GL_MODELVIEW);

	if ((draw_setup & DRAW_SETUP_FIELD) == DRAW_SETUP_FIELD)
		draw_block_field ();

	if ((draw_setup & DRAW_SETUP_ISO) == DRAW_SETUP_ISO)
		draw_block_iso ();

	if ((draw_setup & DRAW_SETUP_MESH) == DRAW_SETUP_MESH)
		draw_block_mesh ();

	draw_block_stroke_lines ();

	if ((draw_setup & DRAW_SETUP_AREA_DIVIDERS) == DRAW_SETUP_AREA_DIVIDERS)
		draw_block_areas ();

	if ((draw_setup & DRAW_SETUP_PLOT) == DRAW_SETUP_PLOT)
		draw_block_plot ();

	if ((draw_setup & DRAW_SETUP_VECTOR_FIELD) == DRAW_SETUP_VECTOR_FIELD)
		draw_block_vector_field ();

	draw_block_white_frame ();
	draw_block_frame ();
	draw_block_strokes ();

	draw_block_legend ();

	draw_name ();

	glutSwapBuffers ();

	// save picture
	CaptureScreenShot ();

	glutLeaveMainLoop ();
}

void Painter::draw_block_field ()
{
	int nodes[3];
	double point[2];

	glPushMatrix ();

	glTranslatef ((GLfloat)((picture.cN[0] + picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] + picture.c0[1]) / 2.0), (GLfloat)0.0);
	glScalef ((GLfloat)((picture.cN[0] - picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] - picture.c0[1]) / 2.0), (GLfloat)1.0);

	glScalef ((GLfloat)(2.0 / (cN[0] - c0[0])), (GLfloat)(2.0 / (cN[1] - c0[1])), (GLfloat)1.0);
	glTranslatef ((GLfloat)(-(cN[0] + c0[0]) / 2.0), (GLfloat)(-(cN[1] + c0[1]) / 2.0), (GLfloat)0.0);

	double value;
	double rgb[3];
	bool in_pic;
	// go through elements 
	for (int i = 0, i_end = mesh->get_n_elements (); i < i_end; i++)
	{
		in_pic = true;
		if (area_to_draw != -1)
		{
			if (mesh->get_area (i) != area_to_draw)
			{
				in_pic = false;
			}
		}
		//if (custom_selection != 0)
		{
			in_pic = false;
			for (int k = 0; k < 3 && !in_pic; k++)
			{
				mesh->get_node_coordinates (mesh->get_node_number (i, k), point);
				if (((c0[0] - 1e-10 < point[0]) && (point[0] < cN[0] + 1e-10)) || 
					((c0[1] - 1e-10 < point[1]) && (point[1] < cN[1] + 1e-10)))
					in_pic = true;
			}

		}
		if (in_pic)
		{
			mesh->get_base_nodes (i, nodes);
			glBegin (GL_TRIANGLES);
			for (int k = 0; k < 3; k++)
			{
				mesh->get_node_coordinates (nodes[k], point);
				{
					task->get_solution_in_point (k_system, point, &value);
					get_color (DRAW_SETUP_FIELD, value, rgb);
					glColor3f ((GLfloat)rgb[0], (GLfloat)rgb[1], (GLfloat)rgb[2]);
				}
				glVertex3f ((GLfloat)point[0], (GLfloat)point[1], 0);
			}
			glEnd ();

			if (task->is_symmetrical () != 0)
			{
				glBegin (GL_TRIANGLES);
				for (int k = 0; k < 3; k++)
				{
					mesh->get_node_coordinates (nodes[k], point);
					point[0] = -point[0];
					{
						task->get_solution_in_point (k_system, point, &value);
						get_color (DRAW_SETUP_FIELD, value, rgb);
						glColor3f ((GLfloat)rgb[0], (GLfloat)rgb[1], (GLfloat)rgb[2]);
					}
					glVertex3f ((GLfloat)point[0], (GLfloat)point[1], 0);
				}
				glEnd ();
			}
		}
	}

	glPopMatrix ();
}

void Painter::draw_block_legend ()
{
	// move into legend block
	glPushMatrix ();
	glTranslatef ((GLfloat)((legend.cN[0] + legend.c0[0]) / 2.0), (GLfloat)((legend.cN[1] + legend.c0[1]) / 2.0), (GLfloat)0.0);
	glScalef ((GLfloat)((legend.cN[0] - legend.c0[0]) / 2.0), (GLfloat)((legend.cN[1] - legend.c0[1]) / 2.0), (GLfloat)1.0);

	double taken_height = 0.0;
	glLineWidth ((GLfloat)1.0);

	if ((draw_setup & DRAW_SETUP_MESH) == DRAW_SETUP_MESH)
	{
		int n_elements = mesh->get_n_elements ();
		int n_nodes = mesh->get_n_nodes ();

		// min, max
		glColor3f ((GLfloat)0.0, (GLfloat)0.0, (GLfloat)0.0);
		char n_el[128];

		sprintf (n_el, "%i nodes, %i elements", n_nodes, n_elements);
		int lenght = (int)strlen (n_el); // amount of characters
		int str_pel_width = 0;
		for (int j = 0; j < lenght; j++)
			str_pel_width += glutBitmapWidth (painter_font, n_el[j]);

		// start position
		double x = -1.0;
		double y = 0.95;

		glColor3f ((GLfloat)0.0, (GLfloat)0.0, (GLfloat)0.0);
		glRasterPos2f ((GLfloat)(x), (GLfloat)(y));
		
		for (int j = 0; j < strlen (n_el); j++)
		{
			glutBitmapCharacter (painter_font, n_el[j]);
		}

	}
	if ((legend_blocks & DATA_BLOCK_MIN_MAX) == DATA_BLOCK_MIN_MAX)
	{
		if ((max_value - min_value > 1e-15) && (
			((draw_setup & DRAW_SETUP_FIELD) == DRAW_SETUP_FIELD) || 
			(((draw_setup & DRAW_SETUP_ISO) == DRAW_SETUP_ISO)/* && (color_scale_iso != -1)*/)))
		{
			double scale_height = 2.0 / (double)height * 25; // last is amount of pixels that i want scale to take
			// one unit block	
			glPushMatrix ();
			glTranslatef ((GLfloat)(0.0), (GLfloat)(1.0 - scale_height / 2.0), (GLfloat)0.0);
			glScalef ((GLfloat)(1.0), (GLfloat)(scale_height / 2.0), (GLfloat)1.0);

			double line_height;
			line_height = (double)(glutBitmapHeight (painter_font) - 4.0);
			line_height *= 2.0 / (double)height;
			line_height /= ((legend.cN[1] - legend.c0[1]) / 2.0);
			line_height /= (scale_height / 2.0);

			int step = 10;
			double f = min_value;
			double x = -0.8;
			double hx = 1.6 / (double)step;
			double hf = (max_value - min_value) / (double)step;
			double rgb[3];
			get_color (DRAW_SETUP_FIELD, f, rgb);

			glBegin (GL_TRIANGLE_STRIP);
			glColor3f ((GLfloat)rgb[0], (GLfloat)rgb[1], (GLfloat)rgb[2]);
			glVertex3f ((GLfloat)(x), (GLfloat)(-1.0), 0);
			glVertex3f ((GLfloat)(x), (GLfloat)(1.0), 0);

			for (int i = 1; i <= step; i++)
			{
				x = -0.8 + hx * i;
				f = min_value + hf * i;

				get_color (DRAW_SETUP_FIELD, f, rgb);
				glColor3f ((GLfloat)rgb[0], (GLfloat)rgb[1], (GLfloat)rgb[2]);
				glVertex3f ((GLfloat)(x), (GLfloat)(-1.0), 0);
				glVertex3f ((GLfloat)(x), (GLfloat)(1.0), 0);
			}
			glEnd ();

			// frame
			glColor3f ((GLfloat)0.0, (GLfloat)0.0, (GLfloat)0.0);
			glBegin (GL_LINE_LOOP);
			glVertex3f ((GLfloat)(-0.8), (GLfloat)(-1.0), 0);
			glVertex3f ((GLfloat)(0.8), (GLfloat)(-1.0), 0);
			glVertex3f ((GLfloat)(0.8), (GLfloat)(1.0), 0);
			glVertex3f ((GLfloat)(-0.8), (GLfloat)(1.0), 0);
			glEnd ();

			// min, max
			glColor3f ((GLfloat)0.0, (GLfloat)0.0, (GLfloat)0.0);
			char n_el[16];

			if (fabs (min_value) < 0.09 && fabs (min_value) > 1e-10)
				sprintf (n_el, "%.1e", min_value);
			else
				sprintf (n_el, "%.2lf", min_value);
			int lenght = (int)strlen (n_el); // amount of characters
			int str_pel_width = 0;
			for (int j = 0; j < lenght; j++)
				str_pel_width += glutBitmapWidth (painter_font, n_el[j]);

			glRasterPos2f ((GLfloat)(-0.8 - str_pel_width / (double)width / ((legend.cN[0] - legend.c0[0]) / 2.0)),
				(GLfloat)(-1.0 - line_height));
			for (int j = 0; j < strlen (n_el); j++)
			{
				glutBitmapCharacter (painter_font, n_el[j]);
			}

			if (fabs (max_value) < 0.09 && fabs (max_value) > 1e-10)
				sprintf (n_el, "%.1e", max_value);
			else
				sprintf (n_el, "%.2lf", max_value);
			lenght = (int)strlen (n_el); // amount of characters
			str_pel_width = 0;
			for (int j = 0; j < lenght; j++)
				str_pel_width += glutBitmapWidth (painter_font, n_el[j]);

			glRasterPos2f ((GLfloat)(0.8 - str_pel_width / (double)width / ((legend.cN[0] - legend.c0[0]) / 2.0)),
				(GLfloat)(-1.0 - line_height));
			for (int j = 0; j < strlen (n_el); j++)
			{
				glutBitmapCharacter (painter_font, n_el[j]);
			}

			line_height = (double)(glutBitmapHeight (painter_font));
			line_height *= 2.0 / (double)height;
			taken_height += scale_height + line_height + 10.0 * 2.0 / (double)height;
			glPopMatrix ();
		}
	}
	if ((legend_blocks & DATA_BLOCK_FULL_ISO_LEGEND) == DATA_BLOCK_FULL_ISO_LEGEND)
	{
		double * isovalues = new double[N_isovalues];
		double dif = (max_value - min_value) / max_value;
		if (fabs (dif) > 1e-25)
		{
			double step = (max_value - min_value) / (N_isovalues + 1);

			for (int i = 0; i < N_isovalues; i++)
			{
				isovalues[i] = min_value + (i + 1) * step;
			}
		}

		double scale_height = 1.0 / ((legend.cN[1] - legend.c0[1]) / 2.0);
		scale_height *= 2.0 / (double)height;
		double line_pix_height = 16.0; // the height of one section
		double line_height = line_pix_height; // the height of one section
		scale_height *= (int)((N_isovalues + 1) / 2) * line_height;

		glPushMatrix ();
		glTranslatef ((GLfloat)(0.0), (GLfloat)(1.0 - taken_height - scale_height / 2.0), (GLfloat)0.0);
		glScalef ((GLfloat)(1.0), (GLfloat)(scale_height / 2.0), (GLfloat)1.0);

		line_height *= 2.0 / (double)height;
		line_height /= ((legend.cN[1] - legend.c0[1]) / 2.0);
		line_height /= (scale_height / 2.0);

		double rgb[3];
		//if (color_scale_iso != COLOR_SCALE_BLACK)
		{
			for (int i = 0; i < N_isovalues; i++)
			{
				get_color (DRAW_SETUP_ISO, isovalues[i], rgb);
				glColor3f ((GLfloat)rgb[0], (GLfloat)rgb[1], (GLfloat)rgb[2]);
				glLineWidth ((GLfloat)iso_line_width);
				glBegin (GL_LINES);
				if (i < ((N_isovalues + 1) / 2))
				{
					glVertex3f ((GLfloat)(-0.9), (GLfloat)(1.0 - line_height / 2.0 - line_height * i), 0);
					glVertex3f ((GLfloat)(-0.6), (GLfloat)(1.0 - line_height / 2.0 - line_height * i), 0);
				}
				else
				{
					glVertex3f ((GLfloat)(0.1),
						(GLfloat)(1.0 - line_height / 2.0 - line_height * (i % ((N_isovalues + 1) / 2))), 0);
					glVertex3f ((GLfloat)(0.4),
						(GLfloat)(1.0 - line_height / 2.0 - line_height * (i % ((N_isovalues + 1) / 2))), 0);
				}
				glEnd ();

				glColor3f ((GLfloat)0.0, (GLfloat)0.0, (GLfloat)0.0);
				char n_el[16];

				if (fabs (isovalues[i]) < 0.09 && fabs (isovalues[i]) > 1e-10)
					sprintf (n_el, "%.1e", isovalues[i]);
				else
					sprintf (n_el, "%.2lf", isovalues[i]);

				double half_height;
				half_height = ((double)glutBitmapHeight (painter_font) / 2.0 - 4.0);
				half_height *= 2.0 / (double)height;
				half_height /= ((legend.cN[1] - legend.c0[1]) / 2.0);
				half_height /= (scale_height / 2.0);

				if (i < ((N_isovalues + 1) / 2))
				{
					glRasterPos2f ((GLfloat)(-0.55), (GLfloat)(1.0 - line_height / 2.0 - line_height * i - half_height));
				}
				else
				{
					glRasterPos2f ((GLfloat)(0.45),
						(GLfloat)(1.0 - line_height / 2.0 - line_height * (i % ((N_isovalues + 1) / 2)) - half_height));
				}
				for (int j = 0; j < strlen (n_el); j++)
				{
					glutBitmapCharacter (painter_font, n_el[j]);
				}
			}
		}

		delete[] isovalues;

		glPopMatrix ();
		taken_height += (int)((N_isovalues + 10) / 2) * line_pix_height * 2.0 / (double)height;
		{
			int N_spec_iso = (int)spec_iso.size ();
			double scale_height = 1.0 / ((legend.cN[1] - legend.c0[1]) / 2.0);
			scale_height *= 2.0 / (double)height;
			double line_pix_height = 16.0; // the height of one section
			double line_height = line_pix_height; // the height of one section
			scale_height *= (int)((N_spec_iso + 1) / 2) * line_height;

			glPushMatrix ();
			glTranslatef ((GLfloat)(0.0), (GLfloat)(1.0 - taken_height - scale_height / 2.0), (GLfloat)0.0);
			glScalef ((GLfloat)(1.0), (GLfloat)(scale_height / 2.0), (GLfloat)1.0);

			line_height *= 2.0 / (double)height;
			line_height /= ((legend.cN[1] - legend.c0[1]) / 2.0);
			line_height /= (scale_height / 2.0);

			double rgb[3];
			for (int i = 0; i < N_spec_iso; i++)
			{
				get_color (DRAW_SETUP_ISO, spec_iso[i], rgb);
				glColor3f ((GLfloat)rgb[0], (GLfloat)rgb[1], (GLfloat)rgb[2]);
				glLineWidth ((GLfloat)(iso_line_width + 2.0));
				glBegin (GL_LINES);
				if (i < ((N_spec_iso + 1) / 2))
				{
					glVertex3f ((GLfloat)(-0.9), (GLfloat)(1.0 - line_height / 2.0 - line_height * i), 0);
					glVertex3f ((GLfloat)(-0.6), (GLfloat)(1.0 - line_height / 2.0 - line_height * i), 0);
				}
				else
				{
					glVertex3f ((GLfloat)(0.1),
						(GLfloat)(1.0 - line_height / 2.0 - line_height * (i % ((N_spec_iso + 1) / 2))), 0);
					glVertex3f ((GLfloat)(0.4),
						(GLfloat)(1.0 - line_height / 2.0 - line_height * (i % ((N_spec_iso + 1) / 2))), 0);
				}
				glEnd ();

				glColor3f ((GLfloat)0.0, (GLfloat)0.0, (GLfloat)0.0);
				char n_el[16];

				if (fabs (spec_iso[i]) < 0.09 && fabs (spec_iso[i]) > 1e-10)
					sprintf (n_el, "%.2e", spec_iso[i]);
				else
					sprintf (n_el, "%.2lf", spec_iso[i]);

				double half_height;
				half_height = ((double)glutBitmapHeight (painter_font) / 2.0 - 4.0);
				half_height *= 2.0 / (double)height;
				half_height /= ((legend.cN[1] - legend.c0[1]) / 2.0);
				half_height /= (scale_height / 2.0);

				if (i < ((N_spec_iso + 1) / 2))
				{
					glRasterPos2f ((GLfloat)(-0.55), (GLfloat)(1.0 - line_height / 2.0 - line_height * i - half_height));
				}
				else
				{
					glRasterPos2f ((GLfloat)(0.45),
						(GLfloat)(1.0 - line_height / 2.0 - line_height * (i % ((N_spec_iso + 1) / 2)) - half_height));
				}
				for (int j = 0; j < strlen (n_el); j++)
				{
					glutBitmapCharacter (painter_font, n_el[j]);
				}
			}

			glPopMatrix ();
			taken_height += (int)((N_spec_iso + 5) / 2) * line_pix_height * 2.0 / (double)height;
		}
	}
	if ((legend_blocks & DATA_BLOCK_PLOT_LEGEND) == DATA_BLOCK_PLOT_LEGEND)
	{
		int N_graphs = (int)plot_data.size () + (int)splines.size ();
		
		double scale_height = 1.0 / ((legend.cN[1] - legend.c0[1]) / 2.0);
		scale_height *= 2.0 / (double)height;
		double line_pix_height = 16.0; // the height of one section
		double line_height = line_pix_height; // the height of one section
		scale_height *= (int)((N_graphs + 1) / 2) * line_pix_height;
		glLineWidth ((GLfloat)(1.0));

		glPushMatrix ();
		glTranslatef ((GLfloat)(0.0), (GLfloat)(1.0 - taken_height - scale_height / 2.0), (GLfloat)0.0);
		glScalef ((GLfloat)(1.0), (GLfloat)(scale_height / 2.0), (GLfloat)1.0);

		line_height *= 2.0 / (double)height;
		line_height /= ((legend.cN[1] - legend.c0[1]) / 2.0);
		line_height /= (scale_height / 2.0);
		glLineWidth (plot_line_width);

		double point[2];
		for (int i = 0; i < N_graphs; i++)
		{
			double rgb[3];
			if (flag_plot_rainbow)
			{
				get_color (i, rgb);
				glColor3f (rgb[0], rgb[1], rgb[2]);
			}
			else
				glColor3f (0.0, 0.0, 0.0);

			glBegin (GL_LINES);
			if (i < ((N_graphs + 1) /*/ 2*/))
			{
				glVertex3f ((GLfloat)(-0.9), (GLfloat)(1.0 - line_height / 2.0 - line_height * i), 0);
				glVertex3f ((GLfloat)(-0.3), (GLfloat)(1.0 - line_height / 2.0 - line_height * i), 0);
			}
			else
			{
				glVertex3f ((GLfloat)(0.1),
					(GLfloat)(1.0 - line_height / 2.0 - line_height * (i % ((N_graphs + 1) / 2))), 0);
				glVertex3f ((GLfloat)(0.7),
					(GLfloat)(1.0 - line_height / 2.0 - line_height * (i % ((N_graphs + 1) / 2))), 0);
			}
			glEnd ();

			glColor3f ((GLfloat)0.0, (GLfloat)0.0, (GLfloat)0.0);
			char n_el[16];

			if (use_handlers)
				sprintf (n_el, "%s", handlers[i]);
			else
				sprintf (n_el, "%i", i + 1);

			double half_height;
			half_height = ((double)glutBitmapHeight (painter_font) / 2.0 - 4.0);
			half_height *= 2.0 / (double)height;
			half_height /= ((legend.cN[1] - legend.c0[1]) / 2.0);
			half_height /= (scale_height / 2.0);

			if (i < ((N_graphs + 1) /*/ 2*/))
			{
				glRasterPos2f ((GLfloat)(-0.25), (GLfloat)(1.0 - line_height / 2.0 - line_height * i - half_height));
			}
			else
			{
				glRasterPos2f ((GLfloat)(0.75),
					(GLfloat)(1.0 - line_height / 2.0 - line_height * (i % ((N_graphs + 1) / 2)) - half_height));
			}
			for (int j = 0; j < strlen (n_el); j++)
			{
				glutBitmapCharacter (painter_font, n_el[j]);
			}
		}

		double reverse_scale[2];
		reverse_scale[0] = 1.0 / ((legend.cN[0] - legend.c0[0]) / 2.0);
		reverse_scale[1] = 1.0 / ((scale_height / 2.0) * ((legend.cN[1] - legend.c0[1]) / 2.0));
		// add figures
		if (!flag_plot_rainbow)
		{
			for (int i = 0; i < N_graphs; i++)
			{
				if (i < ((N_graphs + 1)/* / 2*/))
				{
					point[0] = -0.6;
					point[1] = 1.0 - line_height / 2.0 - line_height * i;
				}
				else
				{
					point[0] = 0.4;
					point[1] = 1.0 - line_height / 2.0 - line_height * (i % ((N_graphs + 1) / 2));
				}
				add_figure (point, reverse_scale, i);
			}
		}
		glPopMatrix ();

		taken_height += (int)((N_graphs + 5) / 2) * line_pix_height * 2.0 / (double)height;
	}
	if ((legend_blocks & DATA_BLOCK_VECTOR_MIN_MAX) == DATA_BLOCK_VECTOR_MIN_MAX)
	{
		double scale_height = 1.0 / ((legend.cN[1] - legend.c0[1]) / 2.0);
		scale_height *= 2.0 / (double)height;
		scale_height *= max_V_pix_length * 2.0;

		double min_length = min_V_pix_length * 2.0 / (double)width;
		min_length /= ((legend.cN[0] - legend.c0[0]) / 2.0);
		double max_length = max_V_pix_length * 2.0 / (double)width;
		max_length /= ((legend.cN[0] - legend.c0[0]) / 2.0);

		double scale_length;
		double s_x, s_y;
		s_x = 2.0 * (double)figure_size / (double)width;
		s_x /= ((legend.cN[0] - legend.c0[0]) / 2.0);

		s_y = 2.0 * (double)figure_size / (double)height;
		s_y /= ((legend.cN[1] - legend.c0[1]) / 2.0);
		s_y /= (scale_height / 2.0);

		char n_el[16];


		double half_height;
		half_height = ((double)glutBitmapHeight (painter_font) / 2.0 - 4.0);
		half_height *= 2.0 / (double)height;
		half_height /= ((legend.cN[1] - legend.c0[1]) / 2.0);
		half_height /= (scale_height / 2.0);

		glLineWidth (1.0);
		glPushMatrix ();
		glTranslatef ((GLfloat)(0.0), (GLfloat)(1.0 - taken_height - scale_height / 2.0), (GLfloat)0.0);
		glScalef ((GLfloat)(1.0), (GLfloat)(scale_height / 2.0), (GLfloat)1.0);
		
		if (fabs (max_V_length - min_V_length) < 1e-7)
		{
			glPushMatrix ();
			glTranslatef ((GLfloat)(-0.1), (GLfloat)(0.0), (GLfloat)(0.0));
			scale_length = (max_length + min_length) / 2.0;
			glScalef ((GLfloat)(scale_length), (GLfloat)(scale_length), 1.0);

			// line
			glBegin (GL_LINES);
			glVertex3f ((GLfloat)(-1.0), (GLfloat)(0.0), 0);
			glVertex3f ((GLfloat)(1.0), (GLfloat)(0.0), 0);
			glEnd ();

			// triangle
			glPushMatrix ();
			glTranslatef ((GLfloat)(1.0), (GLfloat)(0.0), (GLfloat)(0.0));
			glScalef ((GLfloat)(s_x / scale_length), (GLfloat)(s_y / scale_length), 1.0);

			glColor3f (0.0, 0.0, 0.0);
			glBegin (GL_TRIANGLES);
			glVertex3f ((GLfloat)(0.0), (GLfloat)(-1.0), 0);
			glVertex3f ((GLfloat)(2.5), (GLfloat)(0.0), 0);
			glVertex3f ((GLfloat)(0.0), (GLfloat)(1.0), 0);
			glEnd ();
			glPopMatrix ();

			glPopMatrix ();
			
			if (fabs (max_V_length) < 0.09 && fabs (max_V_length) > 1e-10)
				sprintf (n_el, "%.1e", max_V_length);
			else
				sprintf (n_el, "%.1lf", max_V_length);

			glRasterPos2f ((GLfloat)(-0.25), (GLfloat)(-1.0));
			for (int j = 0; j < strlen (n_el); j++)
			{
				glutBitmapCharacter (painter_font, n_el[j]);
			}


		}
		else
		{
			glPushMatrix ();
			glTranslatef ((GLfloat)(-0.8), (GLfloat)(0.0), (GLfloat)(0.0));
			scale_length = (max_length);
			glScalef ((GLfloat)(scale_length), (GLfloat)(scale_length), 1.0);

			// line
			glBegin (GL_LINES);
			glVertex3f ((GLfloat)(-1.0), (GLfloat)(0.0), 0);
			glVertex3f ((GLfloat)(1.0), (GLfloat)(0.0), 0);
			glEnd ();

			// triangle
			glPushMatrix ();
			glTranslatef ((GLfloat)(1.0), (GLfloat)(0.0), (GLfloat)(0.0));
			glScalef ((GLfloat)(s_x / scale_length), (GLfloat)(s_y / scale_length), 1.0);

			glColor3f (0.0, 0.0, 0.0);
			glBegin (GL_TRIANGLES);
			glVertex3f ((GLfloat)(0.0), (GLfloat)(-1.0), 0);
			glVertex3f ((GLfloat)(2.5), (GLfloat)(0.0), 0);
			glVertex3f ((GLfloat)(0.0), (GLfloat)(1.0), 0);
			glEnd ();
			glPopMatrix ();

			glPopMatrix ();
			
			if (fabs (max_V_length) < 0.09 && fabs (max_V_length) > 1e-10)
				sprintf (n_el, "%.1e", max_V_length);
			else
				sprintf (n_el, "%.1lf", max_V_length);

			glRasterPos2f ((GLfloat)(-1.0), (GLfloat)(-1.0));
			for (int j = 0; j < strlen (n_el); j++)
			{
				glutBitmapCharacter (painter_font, n_el[j]);
			}


			glPushMatrix ();
			glTranslatef ((GLfloat)(0.2), (GLfloat)(0.0), (GLfloat)(0.0));
			scale_length = (min_length);
			glScalef ((GLfloat)(scale_length), (GLfloat)(scale_length), 1.0);

			// line
			glBegin (GL_LINES);
			glVertex3f ((GLfloat)(-1.0), (GLfloat)(0.0), 0);
			glVertex3f ((GLfloat)(1.0), (GLfloat)(0.0), 0);
			glEnd ();

			// triangle
			glPushMatrix ();
			glTranslatef ((GLfloat)(1.0), (GLfloat)(0.0), (GLfloat)(0.0));
			glScalef ((GLfloat)(s_x / scale_length), (GLfloat)(s_y / scale_length), 1.0);

			glColor3f (0.0, 0.0, 0.0);
			glBegin (GL_TRIANGLES);
			glVertex3f ((GLfloat)(0.0), (GLfloat)(-1.0), 0);
			glVertex3f ((GLfloat)(2.5), (GLfloat)(0.0), 0);
			glVertex3f ((GLfloat)(0.0), (GLfloat)(1.0), 0);
			glEnd ();
			glPopMatrix ();

			glPopMatrix ();

			if (fabs (min_V_length) < 0.09 && fabs (min_V_length) > 1e-10)
				sprintf (n_el, "%.1e", min_V_length);
			else
				sprintf (n_el, "%.1lf", min_V_length);

			glRasterPos2f ((GLfloat)(0.1), (GLfloat)(-1.0));
			for (int j = 0; j < strlen (n_el); j++)
			{
				glutBitmapCharacter (painter_font, n_el[j]);
			}
		}
		
		glPopMatrix ();
	}
	glPopMatrix ();
}

void Painter::draw_block_vector_field ()
{
	// net 
	glPushMatrix ();
	glTranslatef ((GLfloat)((picture.cN[0] + picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] + picture.c0[1]) / 2.0), (GLfloat)0.0);
	glScalef ((GLfloat)((picture.cN[0] - picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] - picture.c0[1]) / 2.0), (GLfloat)1.0);

	glScalef ((GLfloat)(2.0 / (cN[0] - c0[0])), (GLfloat)(2.0 / (cN[1] - c0[1])), (GLfloat)1.0);
	glTranslatef ((GLfloat)(-(cN[0] + c0[0]) / 2.0), (GLfloat)(-(cN[1] + c0[1]) / 2.0), (GLfloat)0.0);

	double min_length = min_V_pix_length * 2.0 / (double)width;
	min_length /= (2.0 / (cN[0] - c0[0]));
	min_length /= ((picture.cN[0] - picture.c0[0]) / 2.0);
	double max_length = max_V_pix_length * 2.0 / (double)width;
	max_length /= (2.0 / (cN[0] - c0[0]));
	max_length /= ((picture.cN[0] - picture.c0[0]) / 2.0);

	double point[2], step[2], g[2];
	step[0] = (cN[0] - c0[0]) / (double)(vect_n_x + 1);
	step[1] = (cN[1] - c0[1]) / (double)(vect_n_y + 1);
	double V_length, scale_length;

	for (int i = 1; i < vect_n_x + 1; i++)
	{
		point[0] = c0[0] + i * step[0];
		for (int j = 1; j < vect_n_y + 1; j++)
		{
			point[1] = c0[1] + j * step[1];
			
			task->get_derivative (0, 0, point, &(g[0]));
			task->get_derivative (0, 1, point, &(g[1]));
			V_length = sqrt (g[0] * g[0] + g[1] * g[1]);

			if (V_length > max_V_length)
				max_V_length = V_length;
			if (V_length < min_V_length)
				min_V_length = V_length;
		}
	}

	double s_x, s_y;
	s_x = 2.0 * (double)figure_size / (double)width;
	s_x /= ((picture.cN[0] - picture.c0[0]) / 2.0);
	s_x /= (GLfloat)(2.0 / (cN[0] - c0[0]));
	
	s_y = 2.0 * (double)figure_size / (double)height;
	s_y /= ((picture.cN[1] - picture.c0[1]) / 2.0);
	s_y /= (GLfloat)(2.0 / (cN[1] - c0[1]));

	glColor3f ((GLfloat)0.0, (GLfloat)0.0, (GLfloat)0.0);

	double angle;
	for (int i = 1; i < vect_n_x + 1; i++)
	{
		point[0] = c0[0] + i * step[0];
		for (int j = 1; j < vect_n_y + 1; j++)
		{
			point[1] = c0[1] + j * step[1];
			
			task->get_derivative (0, 0, point, &(g[0]));
			task->get_derivative (0, 1, point, &(g[1]));
			V_length = sqrt (g[0] * g[0] + g[1] * g[1]);

			if (anti_gradient == 1)
			{
				g[0] = -g[0];
				g[1] = -g[1];
			}
			
			if (fabs (max_V_length - min_V_length) < 1e-7)
				scale_length = (max_length + min_length) / 2.0;
			else
				scale_length = min_length + (max_length - min_length) * ((V_length - min_V_length) / (max_V_length - min_V_length));			
			
			angle = acos (g[0] / V_length) * 180.0 / PI;
			if ((g[1] / V_length) < 0)
				angle = -angle; // IF ERROR LOOK HERE

			if (V_length > 1.0)
			{
				glPushMatrix ();
				glTranslatef ((GLfloat)(point[0]), (GLfloat)(point[1]), (GLfloat)0.0);
				glScalef ((GLfloat)(scale_length), (GLfloat)(scale_length), (GLfloat)1.0);
				glRotatef ((GLfloat)(angle), (GLfloat)0.0, (GLfloat)0.0, (GLfloat)1.0);

				// line
				glBegin (GL_LINES);
				glVertex3f ((GLfloat)(-1.0), (GLfloat)(0.0), 0);
				glVertex3f ((GLfloat)(1.0), (GLfloat)(0.0), 0);
				glEnd ();

				// triangle
				glPushMatrix ();
				glTranslatef ((GLfloat)(1.0), (GLfloat)(0.0), (GLfloat)(0.0));
				glScalef ((GLfloat)(s_x / scale_length * 0.7), (GLfloat)(s_y / scale_length * 0.7), 1.0);

				glColor3f (0.0, 0.0, 0.0);
				glBegin (GL_TRIANGLES);
				glVertex3f ((GLfloat)(-1.0), (GLfloat)(-1.0), 0);
				glVertex3f ((GLfloat)(1.5), (GLfloat)(0.0), 0);
				glVertex3f ((GLfloat)(-1.0), (GLfloat)(1.0), 0);
				glEnd ();

				glPopMatrix ();
				glPopMatrix ();

				if (task->is_symmetrical () != 0)
				{
					glPushMatrix ();
					glTranslatef ((GLfloat)(-point[0]), (GLfloat)(point[1]), (GLfloat)0.0);
					glScalef ((GLfloat)(scale_length), (GLfloat)(scale_length), (GLfloat)1.0);

					angle = -angle + 180;

					glRotatef ((GLfloat)(angle), (GLfloat)0.0, (GLfloat)0.0, (GLfloat)1.0);

					// line
					glBegin (GL_LINES);
					glVertex3f ((GLfloat)(-1.0), (GLfloat)(0.0), 0);
					glVertex3f ((GLfloat)(1.0), (GLfloat)(0.0), 0);
					glEnd ();

					// triangle
					glPushMatrix ();
					glTranslatef ((GLfloat)(1.0), (GLfloat)(0.0), (GLfloat)(0.0));
					glScalef ((GLfloat)(s_x / scale_length * 0.7), (GLfloat)(s_y / scale_length * 0.7), 1.0);

					glColor3f (0.0, 0.0, 0.0);
					glBegin (GL_TRIANGLES);
					glVertex3f ((GLfloat)(-1.0), (GLfloat)(-1.0), 0);
					glVertex3f ((GLfloat)(1.5), (GLfloat)(0.0), 0);
					glVertex3f ((GLfloat)(-1.0), (GLfloat)(1.0), 0);
					glEnd ();
					glPopMatrix ();
					glPopMatrix ();
				}
			}
		}
	}

	if (task->is_symmetrical () != 0)
	{
		point[0] = 0.0;
		for (int j = vect_n_y / 2; j < vect_n_y + 1; j++)
		{
			point[1] = c0[1] + j * step[1];

			task->get_derivative (0, 0, point, &(g[0]));
			task->get_derivative (0, 1, point, &(g[1]));
			V_length = sqrt (g[0] * g[0] + g[1] * g[1]);

			if (anti_gradient == 1)
			{
				g[0] = -g[0];
				g[1] = -g[1];
			}

			if (fabs (max_V_length - min_V_length) < 1e-7)
				scale_length = (max_length + min_length) / 2.0;
			else
				scale_length = min_length + (max_length - min_length) * ((V_length - min_V_length) / (max_V_length - min_V_length));

			angle = acos (g[0] / V_length) * 180.0 / PI;
			if ((g[1] / V_length) < 0)
				angle = -angle; // IF ERROR LOOK HERE

			if (V_length > 1.0)
			{
				glPushMatrix ();
				glTranslatef ((GLfloat)(point[0]), (GLfloat)(point[1]), (GLfloat)0.0);
				glScalef ((GLfloat)(scale_length), (GLfloat)(scale_length), (GLfloat)1.0);
				glRotatef ((GLfloat)(angle), (GLfloat)0.0, (GLfloat)0.0, (GLfloat)1.0);

				// line
				glBegin (GL_LINES);
				glVertex3f ((GLfloat)(-1.0), (GLfloat)(0.0), 0);
				glVertex3f ((GLfloat)(1.0), (GLfloat)(0.0), 0);
				glEnd ();

				// triangle
				glPushMatrix ();
				glTranslatef ((GLfloat)(1.0), (GLfloat)(0.0), (GLfloat)(0.0));
				glScalef ((GLfloat)(s_x / scale_length * 0.7), (GLfloat)(s_y / scale_length * 0.7), 1.0);

				glColor3f (0.0, 0.0, 0.0);
				glBegin (GL_TRIANGLES);
				glVertex3f ((GLfloat)(-1.0), (GLfloat)(-1.0), 0);
				glVertex3f ((GLfloat)(1.5), (GLfloat)(0.0), 0);
				glVertex3f ((GLfloat)(-1.0), (GLfloat)(1.0), 0);
				glEnd ();

				glPopMatrix ();
				glPopMatrix ();
			}
		}
	}
	glPopMatrix ();
}

void Painter::draw_block_iso ()
{
	double * isovalues = new double[N_isovalues + (int)spec_iso.size()];
	std::vector<IsovalueSection> isolines_data;

	// set isovalues
	double dif = (max_value - min_value);
	if (fabs (max_value) > 1e-10)
		dif /= max_value;

	if (fabs (dif) > 1e-25)
	{
		double step = (max_value - min_value) / (N_isovalues + 1);

		for (int i = 0; i < N_isovalues; i++)
		{
			isovalues[i] = min_value + (i + 1) * step;
		}
		for (int i = N_isovalues, j = 0; j < (int)spec_iso.size (); i++, j++)
		{
			isovalues[i] = spec_iso[j];
		}

		double localMin, localMax;
		double v;
		double point[2];
		int nodes[3];
		int k;
		IsovalueSection iS;
		bool in_pic;
		int counter; 

		// go through triangles
		for (int i = 0, i_end = mesh->get_n_elements (); i < i_end; i++)
		{
			in_pic = false;

			// check if all element nodes are in

			for (int k = 0; k < 3 && !in_pic; k++)
			{
				mesh->get_base_nodes (i, nodes);					
				mesh->get_node_coordinates (nodes[k], point);
				if ( ( (c0[0] - 1e-7 < point[0]) && (point[0] < cN[0] + 1e-7) ) ||
					((c0[1] - 1e-7 < point[1]) && (point[1] < cN[1] + 1e-7)) )
				{
					in_pic = true;
				}
			}				

			if (in_pic)
			{
				// get min/max values in triangles' base nodes
				mesh->get_base_nodes (i, nodes);
				localMin = 1e+20;
				localMax = -1e+20;
				for (int k = 0; k < 3; k++)
				{
					mesh->get_node_coordinates (nodes[k], point);
					{
						if ((draw_setup & DRAW_SETUP_DERIVATIVE) == DRAW_SETUP_DERIVATIVE)
							task->get_derivative (k_system, der_var, point, &v);
						else
							task->get_solution_in_point (k_system, point, &v);

						if (v > localMax)
							localMax = v;
						if (v < localMin)
							localMin = v;
					}
				}
				// check if any isovalue fits into the range
				for (k = 0; k < N_isovalues + (int)spec_iso.size(); k++)
				{
					// if any fits, add it into isovalue_data
					if (localMin < isovalues[k] && isovalues[k] < localMax)
					{
						iS.isovalue = k;
						iS.triangle = i;
						isolines_data.push_back (iS);
					}
				}
			}
		}
		double c1[2], c2[2];
		{
			for (size_t i = 0, i_end = isolines_data.size (); i < i_end; i++)
			{
				bool isoline_section = task->get_isoline_section (k_system, der_var, isolines_data[i].triangle, isovalues[isolines_data[i].isovalue], c1, c2);
				if (isoline_section)
				{
					isolines_data[i].point[0].set_dim (2);
					isolines_data[i].point[0].set_point (c1);

					isolines_data[i].point[1].set_dim (2);
					isolines_data[i].point[1].set_point (c2);
				}
				else
				{
					isolines_data[i].isovalue = -1;
				}
			}
		}
	}

	glPushMatrix ();
	glTranslatef ((GLfloat)((picture.cN[0] + picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] + picture.c0[1]) / 2.0), (GLfloat)0.0);
	glScalef ((GLfloat)((picture.cN[0] - picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] - picture.c0[1]) / 2.0), (GLfloat)1.0);

	glScalef ((GLfloat)(2.0 / (cN[0] - c0[0])), (GLfloat)(2.0 / (cN[1] - c0[1])), (GLfloat)1.0);
	glTranslatef ((GLfloat)(-(cN[0] + c0[0]) / 2.0), (GLfloat)(-(cN[1] + c0[1]) / 2.0), (GLfloat)0.0);

	// draw each section of the same isovalue
	for (size_t i = 0, i_end = isolines_data.size (); i < i_end; i++)
	{
		if (isolines_data[i].isovalue != -1)
		{
			double rgb[3];
			get_color (DRAW_SETUP_ISO, isovalues[isolines_data[i].isovalue], rgb);
		//	if (isolines_data[i].isovalue < N_isovalues)
				glColor3f ((GLfloat)rgb[0], (GLfloat)rgb[1], (GLfloat)rgb[2]);
			//else
			//	glColor3f (0.0, 0.0, 0.0);
			
			
			 if (isolines_data[i].isovalue < N_isovalues)
				glLineWidth ((GLfloat)iso_line_width);
			else
				glLineWidth ((GLfloat)(iso_line_width + 1.0));

			glBegin (GL_LINES);
			glVertex3f ((GLfloat)isolines_data[i].point[0].X (), (GLfloat)isolines_data[i].point[0].Y (), 0);
			glVertex3f ((GLfloat)isolines_data[i].point[1].X (), (GLfloat)isolines_data[i].point[1].Y (), 0);
			glEnd ();

			if (task->is_symmetrical () != 0)
			{
				if (task->is_symmetrical () == 2)
				{
					get_color (DRAW_SETUP_ISO, -isovalues[isolines_data[i].isovalue], rgb);
					glColor3f ((GLfloat)rgb[0], (GLfloat)rgb[1], (GLfloat)rgb[2]);
				}

				glBegin (GL_LINES);
				glVertex3f ((GLfloat)-isolines_data[i].point[0].X (), (GLfloat)isolines_data[i].point[0].Y (), 0);
				glVertex3f ((GLfloat)-isolines_data[i].point[1].X (), (GLfloat)isolines_data[i].point[1].Y (), 0);
				glEnd ();
			}
		}
	}
	glPopMatrix ();

	// clean up
	isolines_data.clear ();

	delete[] isovalues;
}

void Painter::draw_block_white_frame ()
{
	glPushMatrix ();

	glTranslatef ((GLfloat)(0.0), (GLfloat)((-1.0 + picture.c0[1]) / 2.0), (GLfloat)0.0);
	glScalef ((GLfloat)(1.0), (GLfloat)((-1.0 - picture.c0[1]) / 2.0), (GLfloat)1.0);

	glColor3f ((GLfloat)(1.0), (GLfloat)(1.0), (GLfloat)(1.0));
	glBegin (GL_TRIANGLE_FAN);
	glVertex3f ((GLfloat)(-1.0), (GLfloat)(1.0), 0);
	glVertex3f ((GLfloat)(-1.0), (GLfloat)(-1.0), 0);
	glVertex3f ((GLfloat)(1.0), (GLfloat)(-1.0), 0);
	glVertex3f ((GLfloat)(1.0), (GLfloat)(1.0), 0);
	glEnd ();

	glPopMatrix ();

	glPushMatrix ();

	glTranslatef ((GLfloat)((-1.0 + picture.c0[0]) / 2.0), (GLfloat)((1.0 + picture.c0[1]) / 2.0), (GLfloat)0.0);
	glScalef ((GLfloat)((-1.0 - picture.c0[0]) / 2.0), (GLfloat)((1.0 - picture.c0[1]) / 2.0), (GLfloat)1.0);

	glColor3f ((GLfloat)(1.0), (GLfloat)(1.0), (GLfloat)(1.0));
	glBegin (GL_TRIANGLE_FAN);
	glVertex3f ((GLfloat)(-1.0), (GLfloat)(1.0), 0);
	glVertex3f ((GLfloat)(-1.0), (GLfloat)(-1.0), 0);
	glVertex3f ((GLfloat)(1.0), (GLfloat)(-1.0), 0);
	glVertex3f ((GLfloat)(1.0), (GLfloat)(1.0), 0);
	glEnd ();

	glPopMatrix ();


	glPushMatrix ();

	glTranslatef ((GLfloat)((1.0 + picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] + 1.0) / 2.0), (GLfloat)0.0);
	glScalef ((GLfloat)((1.0 - picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] - 1.0) / 2.0), (GLfloat)1.0);

	glColor3f ((GLfloat)(1.0), (GLfloat)(1.0), (GLfloat)(1.0));
	glBegin (GL_TRIANGLE_FAN);
	glVertex3f ((GLfloat)(-1.0), (GLfloat)(1.0), 0);
	glVertex3f ((GLfloat)(-1.0), (GLfloat)(-1.0), 0);
	glVertex3f ((GLfloat)(1.0), (GLfloat)(-1.0), 0);
	glVertex3f ((GLfloat)(1.0), (GLfloat)(1.0), 0);
	glEnd ();

	glPopMatrix ();

	glPushMatrix ();

	glTranslatef ((GLfloat)((picture.cN[0] + 1.0) / 2.0), (GLfloat)((picture.cN[1] + picture.c0[1]) / 2.0), (GLfloat)0.0);
	glScalef ((GLfloat)((picture.cN[0] - 1.0) / 2.0), (GLfloat)((picture.cN[1] - picture.c0[1]) / 2.0), (GLfloat)1.0);

	glColor3f ((GLfloat)(1.0), (GLfloat)(1.0), (GLfloat)(1.0));
	glBegin (GL_TRIANGLE_FAN);
	glVertex3f ((GLfloat)(-1.0), (GLfloat)(1.0), 0);
	glVertex3f ((GLfloat)(-1.0), (GLfloat)(-1.0), 0);
	glVertex3f ((GLfloat)(1.0), (GLfloat)(-1.0), 0);
	glVertex3f ((GLfloat)(1.0), (GLfloat)(1.0), 0);
	glEnd ();

	glPopMatrix ();
}

void Painter::draw_block_mesh ()
{
	int nodes[3];
	double point[2];

	glPushMatrix ();
	glTranslatef ((GLfloat)((picture.cN[0] + picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] + picture.c0[1]) / 2.0), (GLfloat)0.0);
	glScalef ((GLfloat)((picture.cN[0] - picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] - picture.c0[1]) / 2.0), (GLfloat)1.0);

	glScalef ((GLfloat)(2.0 / (cN[0] - c0[0])), (GLfloat)(2.0 / (cN[1] - c0[1])), (GLfloat)1.0);
	glTranslatef ((GLfloat)(-(cN[0] + c0[0]) / 2.0), (GLfloat)(-(cN[1] + c0[1]) / 2.0), (GLfloat)0.0);

	glLineWidth ((GLfloat)mesh_line_width);

	bool in_pic;
	int order;
	// go through elements 
	for (int i = 0, i_end = mesh->get_n_elements (); i < i_end; i++)
	{
		in_pic = true;
		if (area_to_draw != -1)
		{
			if (mesh->get_area (i) != area_to_draw)
			{
				in_pic = false;
			}
		}
		if (custom_selection != 0)
		{
			// if element node are not in the selection
			in_pic = false;
			for (int k = 0; k < 3 && !in_pic; k++)
			{
				mesh->get_node_coordinates (mesh->get_node_number (i, k), point);
				if ((c0[0] - 1e-10 < point[0]) && (point[0] < cN[0] + 1e-10) && (c0[1] - 1e-10 < point[1]) && (point[1] < cN[1] + 1e-10))
					in_pic = true;
			}				

		}
		if (in_pic)
		{
			mesh->get_base_nodes (i, nodes);

			glColor3f ((GLfloat)1.0, (GLfloat)1.0, (GLfloat)1.0);
			order = mesh->get_area (i);
			if (order == 2)
				glColor3f ((GLfloat)1.0, (GLfloat)0.82, (GLfloat)0.88);
			
			// fill element
			double center[2] = {0.0, 0.0};
			{
				glBegin (GL_TRIANGLES);
				for (int k = 0; k < 3; k++)
				{
					mesh->get_node_coordinates (nodes[k], point);
					glVertex3f ((GLfloat)(point[0]), (GLfloat)(point[1]), 0);
					for (int j = 0; j < 2; j++)
						center[j] += point[j];
				}
				glEnd ();
			}

			// outline it
			glColor3f ((GLfloat)0.0, (GLfloat)0.0, (GLfloat)0.0);
			glBegin (GL_LINE_LOOP);
			for (int k = 0; k < 3; k++)
			{
				mesh->get_node_coordinates (nodes[k], point);
				glVertex3f ((GLfloat)(point[0]), (GLfloat)(point[1]), 0);
			}
			glEnd ();

			// number 
			if ((mesh_data & DRAW_SETUP_MESH_MARKERS) == DRAW_SETUP_MESH_MARKERS)
			{
				for (int j = 0; j < 2; j++)
					center[j] /= 3.0;

				char n_el[16];
				sprintf (n_el, "%i", i);
				glColor3f ((GLfloat)0.2, (GLfloat)0.25, (GLfloat)0.77);
				glRasterPos2f ((GLfloat)(center[0]), (GLfloat)(center[1]));
				for (int j = 0; j < strlen (n_el); j++)
				{
					glutBitmapCharacter (painter_font, n_el[j]);
				}
			}
		}
	}

	// nodes
	if ((mesh_data & DRAW_SETUP_MESH_MARKERS) == DRAW_SETUP_MESH_MARKERS)
	{
		double point[2];
		for (int i = 0, i_end = mesh->get_n_nodes (); i < i_end; i++)
		{
			mesh->get_node_coordinates (i, point);

			glColor3f ((GLfloat)0.0, (GLfloat)0.0, (GLfloat)0.0);
			glPointSize (point_size);
			glBegin (GL_POINTS);
			glVertex3f (point[0], point[1], 0.0);
			glEnd ();

			char n_el[16];
			sprintf (n_el, "%i", i);
			glColor3f ((GLfloat)0.73, (GLfloat)0.24, (GLfloat)0.68);
			glRasterPos2f ((GLfloat)(point[0]), (GLfloat)(point[1]));
			for (int j = 0; j < strlen (n_el); j++)
			{
				glutBitmapCharacter (painter_font, n_el[j]);
			}
		}
	}
	glPopMatrix ();
}

void Painter::draw_block_plot ()
{
	glPushMatrix ();
	glTranslatef ((GLfloat)((picture.cN[0] + picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] + picture.c0[1]) / 2.0), (GLfloat)0.0);
	glScalef ((GLfloat)((picture.cN[0] - picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] - picture.c0[1]) / 2.0), (GLfloat)1.0);

	glScalef ((GLfloat)(2.0 / (cN[0] - c0[0])), (GLfloat)(2.0 / (cN[1] - c0[1])), (GLfloat)1.0);
	glTranslatef ((GLfloat)(-(cN[0] + c0[0]) / 2.0), (GLfloat)(-(cN[1] + c0[1]) / 2.0), (GLfloat)0.0);

	// draw pointed
	glLineWidth (plot_line_width);
	for (size_t i = 0, i_end = plot_data.size (); i < i_end; i++)
	{
		double rgb[3];
		if (flag_plot_rainbow)
		{
			get_color (i, rgb);
			glColor3f (rgb[0], rgb[1], rgb[2]);
		}
		else
			glColor3f (0.0, 0.0, 0.0);
		{
			glBegin (GL_LINE_STRIP);
			// draw plot from point to the point 
			for (size_t j = 0, j_end = plot_data[i].size (); j < j_end; j++)
			{
				glVertex3f ((GLfloat)(plot_data[i][j].first), (GLfloat)(plot_data[i][j].second), 0);
			}
			glEnd ();
		}
	}
	// draw splines
	for (size_t i = 0, i_end = splines.size (); i < i_end; i++)
	{
		double rgb[3];
		if (flag_plot_rainbow)
		{
			get_color (i, rgb);
			glColor3f (rgb[0], rgb[1], rgb[2]);
		}
		else
			glColor3f (0.0, 0.0, 0.0);

		double h = (cN[0] - c0[0]) / (0.7 * width);
		double x[1];
		double value = 0.0;

		glBegin (GL_LINE_STRIP);
		// draw plot from point to the point 
		for (int j = 0, j_end = (int)((cN[0] - c0[0]) / h) + 1; j < j_end; j++)
		{
			x[0] = c0[0] + j * h;
			if (splines[i]->get_solution_in_point (x, &value))
				glVertex3f ((GLfloat)(x[0]), (GLfloat)(value), 0);
		}
		glEnd ();
	}

	// add markers
	int counter = 0;
	double reverse_scale[2];
	reverse_scale[0] = 1.0 / (((picture.cN[0] - picture.c0[0]) / 2.0) * (2.0 / (cN[0] - c0[0])));
	reverse_scale[1] = 1.0 / (((picture.cN[1] - picture.c0[1]) / 2.0) * (2.0 / (cN[1] - c0[1])));

	// for pointed
	if (!flag_plot_rainbow)
	{
		for (size_t i = 0, i_end = plot_data.size (); i < i_end; i++)
		{
			double figure_point[2];
			// go by points
			for (size_t j = 0, j_end = plot_data[i].size (); j < j_end; j++)
			{
				figure_point[0] = plot_data[i][j].first;
				figure_point[1] = plot_data[i][j].second;
				add_figure (figure_point, reverse_scale, counter);
			}
			counter++;
		}
	}
	double x_axis_step, y_axis_step;
	scale_axes (cN[0], cN[1], c0[0], c0[1], &x_axis_step, &y_axis_step);
	// for splines
	double margin = 0.05;
	if (!flag_plot_rainbow)
	{
		for (size_t i = 0, i_end = splines.size (); i < i_end; i++)
		{
			// go by x
			double c[2];
			double figure_point[2];
			// start point
			splines[i]->get_boundaries (figure_point, c);
			if (splines[i]->get_solution_in_point (figure_point, &figure_point[1]))
			{
				add_figure (figure_point, reverse_scale, counter);
			}

			// points between
			figure_point[0] = c0[0];
			for (int j = 0; figure_point[0] < cN[0] + 1e-5; j++)
			{
				figure_point[0] = c0[0] + j * x_axis_step;
				// get y value
				if (splines[i]->get_solution_in_point (figure_point, &figure_point[1]))
				{
					add_figure (figure_point, reverse_scale, counter);
				}
			}

			// end point
			splines[i]->get_boundaries (c, figure_point);
			if (splines[i]->get_solution_in_point (figure_point, &figure_point[1]))
			{
				add_figure (figure_point, reverse_scale, counter);
			}
			counter++;
		}
	}
	glPopMatrix ();
}

void Painter::draw_block_areas ()
{
	std::vector <std::pair<int, int>> dividers_nodes;
	mesh->get_area_dividers (&dividers_nodes);

	double c[2];
	glPushMatrix ();
	glTranslatef ((GLfloat)((picture.cN[0] + picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] + picture.c0[1]) / 2.0), (GLfloat)0.0);
	glScalef ((GLfloat)((picture.cN[0] - picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] - picture.c0[1]) / 2.0), (GLfloat)1.0);

	glScalef ((GLfloat)(2.0 / (cN[0] - c0[0])), (GLfloat)(2.0 / (cN[1] - c0[1])), (GLfloat)1.0);
	glTranslatef ((GLfloat)(-(cN[0] + c0[0]) / 2.0), (GLfloat)(-(cN[1] + c0[1]) / 2.0), (GLfloat)0.0);
	
	glLineWidth ((GLfloat)(1.5));
	glColor3f ((GLfloat)0.0, (GLfloat)0.0, (GLfloat)0.0);
	for (int i = 0; i < (int)dividers_nodes.size (); i++)
	{
		glBegin (GL_LINES);
		mesh->get_node_coordinates (dividers_nodes[i].first, c);
		if ((c0[0] - 1e-7 < c[0] && c[0] < cN[0] + 1e-7) &&
			(c0[1] - 1e-7 < c[1] && c[1] < cN[1] + 1e-7))
		{
			glVertex3f ((GLfloat)(c[0]), (GLfloat)(c[1]), 0);
			if (task->is_symmetrical () != 0)
			{
				glVertex3f ((GLfloat)(-c[0]), (GLfloat)(c[1]), 0);
			}
		}
		mesh->get_node_coordinates (dividers_nodes[i].second, c);
		if ((c0[0] - 1e-7 < c[0] && c[0] < cN[0] + 1e-7) &&
			(c0[1] - 1e-7 < c[1] && c[1] < cN[1] + 1e-7))
		{
			glVertex3f ((GLfloat)(c[0]), (GLfloat)(c[1]), 0);
			if (task->is_symmetrical () != 0)
			{
				glVertex3f ((GLfloat)(-c[0]), (GLfloat)(c[1]), 0);
			}
		}
		glEnd ();
	}
	glPopMatrix ();

	glLineWidth ((GLfloat)1.0);
}

void Painter::draw_block_frame ()
{
	// frame
	glLineWidth ((GLfloat)(1.0));

	glPushMatrix ();
	glTranslatef ((GLfloat)((picture.cN[0] + picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] + picture.c0[1]) / 2.0), (GLfloat)0.0);
	glScalef ((GLfloat)((picture.cN[0] - picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] - picture.c0[1]) / 2.0), (GLfloat)1.0);

	glColor3f ((GLfloat)0.0, (GLfloat)0.0, (GLfloat)0.0);
	glBegin (GL_LINE_LOOP);
	glVertex3f ((GLfloat)(-1.0), (GLfloat)(-1.0), 0);
	glVertex3f ((GLfloat)(1.0), (GLfloat)(-1.0), 0);
	glVertex3f ((GLfloat)(1.0), (GLfloat)(1.0), 0);
	glVertex3f ((GLfloat)(-1.0), (GLfloat)(1.0), 0);
	glEnd ();

	glPopMatrix ();
}

void Painter::draw_name ()
{
	// text
	int lenght = (int)strlen (picture_name); // amount of characters in the name
	double x, y; // starting text position

	int str_pel_width = 0;
	for (int j = 0; j < lenght; j++)
		str_pel_width += glutBitmapWidth (painter_font, picture_name[j]);

	// start position
	x = 2.0 * (double)(Frame + pic_width - str_pel_width) / (double)width - 1.0;
	y = picture.cN[1] + 2.0 * 4.0 / (double)height;

	glColor3f ((GLfloat)0.0, (GLfloat)0.0, (GLfloat)0.0);
	glRasterPos2f ((GLfloat)(x), (GLfloat)(y));

	for (int j = 0; j < lenght; j++)
	{
		glutBitmapCharacter (painter_font, picture_name[j]);
	}
}

void Painter::draw_block_strokes ()
{
	double x_axis_step, y_axis_step;
	scale_axes (cN[0], cN[1], c0[0], c0[1], &x_axis_step, &y_axis_step);
	glLineWidth ((GLfloat)(1.0));

	// sections on axises and writings
	{
		// by x
		{
			glPushMatrix ();
			glTranslatef ((GLfloat)((picture.cN[0] + picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] + picture.c0[1]) / 2.0), (GLfloat)0.0);
			glScalef ((GLfloat)((picture.cN[0] - picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] - picture.c0[1]) / 2.0), (GLfloat)1.0);

			glScalef ((GLfloat)(2.0 / (cN[0] - c0[0])), (GLfloat)(1.0), (GLfloat)1.0);
			glTranslatef ((GLfloat)(-(cN[0] + c0[0]) / 2.0), (GLfloat)(0.0), (GLfloat)0.0);

			int n = (int)(log10 (fabs (x_axis_step)));
			double x_start = (double)((int)round (c0[0] / pow (10.0, n)) * x_axis_step);
			x_start = c0[0];
			if ((draw_setup & DRAW_SETUP_PLOT) != DRAW_SETUP_PLOT)
				if (task != NULL)
					if (task->is_symmetrical ())
						x_start = -cN[0];
			double x = x_start;
			for (int i = 0; x < cN[0] + 1e-10; i++)
			{

				// sections
				glColor3f (0.0, 0.0, 0.0);
				glBegin (GL_LINES);
				glVertex3f ((GLfloat)(x), (GLfloat)(-1.0 - 3.0 * 2.0 / (double)height), 0);
				glVertex3f ((GLfloat)(x), (GLfloat)(-1.0 + 3.0 * 2.0 / (double)height), 0);
				glEnd ();
				// text under it
				glColor3f ((GLfloat)0.0, (GLfloat)0.0, (GLfloat)0.0);
				char n_el[16];

				if ((fabs (x_axis_step) < 0.009) || (fabs (x_axis_step) > 9e+3))
					sprintf (n_el, "%.1e", x);
				else
					sprintf (n_el, "%.2lf", x);

				int lenght = (int)strlen (n_el); // amount of characters
				int str_pel_width = 0;
				for (int j = 0; j < lenght; j++)
					str_pel_width += glutBitmapWidth (painter_font, n_el[j]);

				double half_length;
				half_length = 1.0 / (2.0 / (cN[0] - c0[0]));
				half_length /= ((picture.cN[0] - picture.c0[0]) / 2.0);
				half_length *= 2.0 / (double)width;
				half_length *= (double)(str_pel_width) / 2.0;

				double height_gap;
				height_gap = 1.0 / ((picture.cN[1] - picture.c0[1]) / 2.0);
				height_gap *= 2.0 / (double)height;
				height_gap *= glutBitmapHeight (painter_font);

				glRasterPos2f ((GLfloat)(x - half_length),
					(GLfloat)(-1.0 - height_gap));

				for (int j = 0; j < strlen (n_el); j++)
				{
					glutBitmapCharacter (painter_font, n_el[j]);
				}
				x = x_start + i * x_axis_step;
			}

			glPopMatrix ();

			glRasterPos2f ((GLfloat)(picture.cN[0] + 10.0 / (double)(width)), (GLfloat)(picture.c0[1] + 2.0 / (double)(height)));
			for (int j = 0; j < strlen (axis_name_x); j++)
			{
				glutBitmapCharacter (painter_font, axis_name_x[j]);
			}

		}

		// by y
		{
			glPushMatrix ();
			glTranslatef ((GLfloat)((picture.cN[0] + picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] + picture.c0[1]) / 2.0), (GLfloat)0.0);
			glScalef ((GLfloat)((picture.cN[0] - picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] - picture.c0[1]) / 2.0), (GLfloat)1.0);
			glScalef ((GLfloat)(1.0), (GLfloat)(2.0 / (cN[1] - c0[1])), (GLfloat)1.0);
			glTranslatef ((GLfloat)(0.0), (GLfloat)(-(cN[1] + c0[1]) / 2.0), (GLfloat)0.0);

			// get log 10 (a[i]), rounded up
			int n = (int)(log10 (fabs (y_axis_step)));
			double y_start = (double)((int)round (c0[1] / pow (10.0, n)) * y_axis_step);
			y_start = c0[1];
			double y = y_start;
			for (int i = 0; y < cN[1] + 1e-10; i++)
			{
				// sections
				glColor3f (0.0, 0.0, 0.0);
				glBegin (GL_LINES);
				glVertex3f ((GLfloat)(-1.0 - 3.0 * 2.0 / (double)width), (GLfloat)(y), 0);
				glVertex3f ((GLfloat)(-1.0 + 3.0 * 2.0 / (double)width), (GLfloat)(y), 0);
				glEnd ();
				// text under it
				glColor3f (0.0, 0.0, 0.0);
				char n_el[16];

				if ((fabs (y_axis_step) < 0.009) || (fabs (y_axis_step) > 9e+3))
					sprintf (n_el, "%.2e", y);
				else
					sprintf (n_el, "%.2lf", y);

				int lenght = (int)strlen (n_el); // amount of characters
				int str_pel_width = 0;
				for (int j = 0; j < lenght; j++)
					str_pel_width += glutBitmapWidth (painter_font, n_el[j]);

				double line_width;
				line_width = 1.0 / ((picture.cN[0] - picture.c0[0]) / 2.0);
				line_width *= 2.0 / (double)width;
				line_width *= (double)(str_pel_width + 4);

				double half_height;
				half_height = ((double)glutBitmapHeight (painter_font) / 2.0 - 4.0);
				half_height *= 2.0 / (double)height;
				half_height /= (2.0 / (cN[1] - c0[1]));
				half_height /= ((picture.cN[1] - picture.c0[1]) / 2.0);

				glRasterPos2f ((GLfloat)(-1.0 - line_width), (GLfloat)(y - half_height));

				for (int j = 0; j < strlen (n_el); j++)
				{
					glutBitmapCharacter (painter_font, n_el[j]);
				}
				y = y_start + i * y_axis_step;
			}
			glPopMatrix ();


			glRasterPos2f ((GLfloat)(picture.c0[0] - 10.0 / (double)(width)), (GLfloat)(picture.cN[1] + 10.0 / (double)(height)));
			for (int j = 0; j < strlen (axis_name_y); j++)
			{
				glutBitmapCharacter (painter_font, axis_name_y[j]);
			}
		}
	}
}

void Painter::draw_block_stroke_lines ()
{
	double x_axis_step, y_axis_step;
	scale_axes (cN[0], cN[1], c0[0], c0[1], &x_axis_step, &y_axis_step);
	glLineWidth ((GLfloat)(1.0));

	// sections on axises and writings
	{
		// by x
		{
			glPushMatrix ();
			glTranslatef ((GLfloat)((picture.cN[0] + picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] + picture.c0[1]) / 2.0), (GLfloat)0.0);
			glScalef ((GLfloat)((picture.cN[0] - picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] - picture.c0[1]) / 2.0), (GLfloat)1.0);

			glScalef ((GLfloat)(2.0 / (cN[0] - c0[0])), (GLfloat)(1.0), (GLfloat)1.0);
			glTranslatef ((GLfloat)(-(cN[0] + c0[0]) / 2.0), (GLfloat)(0.0), (GLfloat)0.0);

			int n = (int)(log10 (fabs (x_axis_step)));
			double x_start = (double)((int)round (c0[0] / pow (10.0, n)) * x_axis_step);
			x_start = c0[0];
			if ((draw_setup & DRAW_SETUP_PLOT) != DRAW_SETUP_PLOT)
				if (task != NULL)
					if (task->is_symmetrical ())
						x_start = -cN[0];
			double x = x_start;
			for (int i = 0; x < cN[0] + 1e-10; i++)
			{
				// lines
				if ((draw_setup & DRAW_SETUP_PLOT) == DRAW_SETUP_PLOT)
				{
					glColor3f ((GLfloat)0.9, (GLfloat)0.9, (GLfloat)0.9);
					glBegin (GL_LINES);
					glVertex3f ((GLfloat)(x), (GLfloat)(1.0), 0);
					glVertex3f ((GLfloat)(x), (GLfloat)(-1.0), 0);
					glEnd ();
				}

				x = x_start + i * x_axis_step;
			}

			glPopMatrix ();

			glRasterPos2f ((GLfloat)(picture.cN[0] + 10.0 / (double)(width)), (GLfloat)(picture.c0[1] + 2.0 / (double)(height)));
			for (int j = 0; j < strlen (axis_name_x); j++)
			{
				glutBitmapCharacter (painter_font, axis_name_x[j]);
			}

		}

		// by y
		{
			glPushMatrix ();
			glTranslatef ((GLfloat)((picture.cN[0] + picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] + picture.c0[1]) / 2.0), (GLfloat)0.0);
			glScalef ((GLfloat)((picture.cN[0] - picture.c0[0]) / 2.0), (GLfloat)((picture.cN[1] - picture.c0[1]) / 2.0), (GLfloat)1.0);
			glScalef ((GLfloat)(1.0), (GLfloat)(2.0 / (cN[1] - c0[1])), (GLfloat)1.0);
			glTranslatef ((GLfloat)(0.0), (GLfloat)(-(cN[1] + c0[1]) / 2.0), (GLfloat)0.0);

			// get log 10 (a[i]), rounded up
			int n = (int)(log10 (fabs (y_axis_step)));
			double y_start = (double)((int)round (c0[1] / pow (10.0, n)) * y_axis_step);
			y_start = c0[1];
			double y = y_start;
			for (int i = 0; y < cN[1] + 1e-10; i++)
			{
				// lines
				if ((draw_setup & DRAW_SETUP_PLOT) == DRAW_SETUP_PLOT)
				{
					glColor3f ((GLfloat)0.9, (GLfloat)0.9, (GLfloat)0.9);
					glBegin (GL_LINES);
					glVertex3f ((GLfloat)(1.0), (GLfloat)(y), 0);
					glVertex3f ((GLfloat)(-1.0), (GLfloat)(y), 0);
					glEnd ();
				}
				
				y = y_start + i * y_axis_step;
			}
			glPopMatrix ();


			glRasterPos2f ((GLfloat)(picture.c0[0] - 10.0 / (double)(width)), (GLfloat)(picture.cN[1] + 10.0 / (double)(height)));
			for (int j = 0; j < strlen (axis_name_y); j++)
			{
				glutBitmapCharacter (painter_font, axis_name_y[j]);
			}
		}
	}
}

void Painter::scale_axes (double pxmax, double pymax, double pxmin, double pymin, double * x_axis_step, double * y_axis_step)
{
	// sets pxmin, pxmax, pymin, pymax, x_axis_step, y_axis_step to make a pretty picture
	if (custom_scale == 0)
	{
		double a[] = { pxmin, pxmax, pymin, pymax };
		double dif[] = { pxmax - pxmin, pymax - pymin };
		int i_dif = 0;
		for (int i = 0; i < 4; i++)
		{
			if (i == 2)
				i_dif = 1;
			// do nothing with 0
			if (fabs (a[i]) > 1e-10)
			{
				// get log 10 (a[i]), rounded up
				int n = (int)(log10 (fabs (a[i])));
				// get 10 ^ n
				double b = pow (10, n);
				b *= (int)(round)(a[i] / b);

				// adjust it to cover a[i] properly
				// if dif is bigger than half the 10^n than n
				double step_adj;
				if (dif[i_dif] / pow (10, n) < 5 - 1e-5)
				{
					step_adj = pow (10, n - 1);
				}
				else
				{
					step_adj = pow (10, n);
				}
				// else n - 1
				double bs = b;
				if (b > a[i])
				{
					for (int j = 0; j < 10 && bs > a[i]; j++)
					{
						bs = b - j * step_adj;
					}
				}
				else
				{
					for (int j = 0; j < 10 && bs < a[i]; j++)
					{
						bs = b + j * step_adj;
					}
				}
				a[i] = bs;
			}
		}
		//pxmin = a[0];
		//pxmax = a[1];
		//pymin = a[2];
		//pymax = a[3];

		dif[0] = a[1] - a[0];
		dif[1] = a[3] - a[2];
		// gets steps
		// get difference
		int	n = (int)round (log10 (fabs (dif[0]))) - 1;

		if (fabs (pxmax) > 99.0 || fabs (pxmin) > 99.0)
			*x_axis_step = pow (10, n) * 2.0;
		else
			*x_axis_step = pow (10, n);
		n = (int)round (log10 (fabs (dif[1]))) - 1;
		*y_axis_step = pow (10, n);

		//*y_axis_step = 0.05;
	}
	else
	{
		*x_axis_step = (pxmax - pxmin) / n_x_strokes;
		*y_axis_step = (pymax - pymin) / n_y_strokes;
	}
}

Painter::Painter ()
{
	pic_counter = 0;
	width = height = 450;
	scale_x = scale_y = 1.0;
	k_system = 0;
	min_value = -1e+10;
	max_value = 1e+10;
	custom_min_max = 0;
	painter_font = GLUT_BITMAP_HELVETICA_12;

	wchar_t pic_name[128];
	swprintf (pic_name, 128, L"Pictures//img%i.png", pic_counter);
	szDestFile = pic_name;
	axis_name_x = "";
	axis_name_y = "";
	picture_name = " ";
	legend_blocks = DATA_BLOCK_MIN_MAX;
	draw_setup = 0;
	der_var = -1;

	vect_n_x = vect_n_y = 10;
	anti_gradient = 0;
	min_V_pix_length = 4;
	max_V_pix_length = 8;
	max_V_length = -1e+20;
	min_V_length = 1e+20;

	Frame = 60;
	//W_legend = 150;
	W_legend = 300;

	custom_scale = 0;
	custom_selection = 0;
	n_x_strokes = n_y_strokes = 1;

	task = NULL;
	mesh = NULL;
	color_scale_iso = COLOR_SCALE_BLACK;
	color_scale_field = COLOR_SCALE_RAINBOW;
	area_to_draw = -1;
	mesh_data = 0;
	N_isovalues = 0;
	mesh_line_width = 1.0;
	iso_line_width = 1.0;
	spec_iso.clear ();
	use_handlers = false;
	flag_plot_rainbow = false;
	figure_size = 5;
	point_size = 1.0;

	plot_line_width = 1.0;
}

Painter::~Painter ()
{
	if (task != NULL)
		task = NULL;
	if (mesh != NULL)
		mesh = NULL;
}

void Painter::display ()
{
	currentPainter->draw_blocks ();
}

void Painter::set_max_resolution (int max_res)
{
	width = height = max_res;
}

void Painter::set_scale (double Scale_x, double Scale_y)
{
	scale_x = Scale_x;
	scale_y = Scale_y;
}

void Painter::set_mesh_line_width (double Mesh_line_width)
{
	mesh_line_width = Mesh_line_width;
}

void Painter::set_iso_line_width (double Iso_line_width)
{
	iso_line_width = Iso_line_width;
}

void Painter::get_color (int setup, double value, double * rgb)
{
	if (value > max_value)
		value = max_value;
	if (value < min_value)
		value = min_value;

	int color_scale = -1;

	if (setup == DRAW_SETUP_ISO)
	{
		color_scale = color_scale_iso;
	}
	if (setup == DRAW_SETUP_FIELD)
	{
		color_scale = color_scale_field;
	}

	// switch by color scheme
	switch (color_scale)
	{
	case COLOR_SCALE_RED_YELLOW:
	{
		rgb[0] = 1;
		rgb[1] = (max_value - value) / (max_value - min_value);
		rgb[2] = 0;
		break;
	}
	case COLOR_SCALE_RAINBOW:
	{
		int cmin, addcmax;
		cmin = 0;
		addcmax = 255;

		double nrv = 1.0 - (max_value - value) / (max_value - min_value);
		if (fabs (max_value - min_value) < 1e-10)
			nrv = 1.0;
		double hue = 240.0 * nrv;
		double saturation = 1.0;
		double luminocity = 0.5;
		double c = (1.0 - fabs (2.0 * luminocity - 1.0)) * saturation;
		double xc = c * (1.0 - fabs ((hue / 60.0) - floor (hue / 120.0) * 2.0 - 1.0));
		double m = luminocity - c / 2.0;
		double cr, cg, cb;

		if (hue < 60)
		{
			cr = 0;
			cg = xc;
			cb = c;

		}
		if (60 <= hue && hue < 120)
		{
			cr = 0;
			cg = c;
			cb = xc;

		}
		if (120 <= hue && hue < 180)
		{
			cr = xc;
			cg = c;
			cb = 0;
		}
		if (180 <= hue)
		{
			cr = c;
			cg = xc;
			cb = 0;
		}
		if (cr < 0.0) cr = 0.0;
		if (cr > 1.0) cr = 1.0;
		if (cg < 0.0) cg = 0.0;
		if (cg > 1.0) cg = 1.0;
		if (cb < 0.0) cb = 0.0;
		if (cb > 1.0) cb = 1.0;

		rgb[0] = cr + m;
		rgb[1] = cg + m;
		rgb[2] = cb + m;
		break;
	}
	case COLOR_SCALE_BLACK_WHITE:
	{
		rgb[0] = rgb[1] = rgb[2] = (max_value - value) / (max_value - min_value);
		break;
	}
	case COLOR_SCALE_GREY:
	{
		rgb[0] = rgb[1] = rgb[2] = 0.3;
		break;
	}
	default:
	{
		rgb[0] = rgb[1] = rgb[2] = 0.0;
		break;
	}
	}
}

void Painter::get_color (int number, double * rgb)
{
	switch (number)
	{
	case 1:
		rgb[0] = 0.0;
		rgb[1] = 166.0;
		rgb[2] = 0.0;
		break;
	case 2:
		rgb[0] = 184.0;
		rgb[1] = 3.0;
		rgb[2] = 3.0;
		break;
	case 3:
		rgb[0] = 0.0;
		rgb[1] = 121.0;
		rgb[2] = 242.0;
		break;
	case 4:
		rgb[0] = 242.0;
		rgb[1] = 97.0;
		rgb[2] = 0.0;
		break;
	case 5:
		rgb[0] = 97.0;
		rgb[1] = 30.0;
		rgb[2] = 111.0;
		break;
	case 7:
		rgb[0] = 128.0;
		rgb[1] = 255.0;
		rgb[2] = 0.0;
		break;
	case 6:
		rgb[0] = 254.0;
		rgb[1] = 5.0;
		rgb[2] = 142.0;
		break;
	case 8:
		rgb[0] = 97.0;
		rgb[1] = 30.0;
		rgb[2] = 111.0;
		break;
	default:
		rgb[0] = rgb[1] = rgb[2] = 0.0;
	}

	rgb[0] /= 255.0;
	rgb[1] /= 255.0;
	rgb[2] /= 255.0;
}

void Painter::add_figure (double * point, double * reverse_scale, int counter)
{
	double s_x, s_y;
	s_x = 2.0 * (double)figure_size / (double)width;
	s_x *= reverse_scale[0];
	s_y = 2.0 * (double)figure_size / (double)height;
	s_y *= reverse_scale[1];

	glPushMatrix ();

	glTranslatef ((GLfloat)(point[0]), (GLfloat)(point[1]), (GLfloat)(0.0));
	glScalef ((GLfloat)s_x, (GLfloat)s_y, 1.0);

	switch (counter)
	{
	case 0: // black rectangle
	{
		glColor3f (0.0, 0.0, 0.0);
		glBegin (GL_POLYGON);
		glVertex3f ((GLfloat)(-1.0), (GLfloat)(-1.0), 0);
		glVertex3f ((GLfloat)(1.0), (GLfloat)(-1.0), 0);
		glVertex3f ((GLfloat)(1.0), (GLfloat)(1.0), 0);
		glVertex3f ((GLfloat)(-1.0), (GLfloat)(1.0), 0);
		glEnd ();
		break;
	}
	case 1: // white rectangle
	{
		glColor3f ((GLfloat)1.0, (GLfloat)1.0, (GLfloat)1.0);
		glBegin (GL_POLYGON);
		glVertex3f ((GLfloat)(-1.0), (GLfloat)(-1.0), 0);
		glVertex3f ((GLfloat)(1.0), (GLfloat)(-1.0), 0);
		glVertex3f ((GLfloat)(1.0), (GLfloat)(1.0), 0);
		glVertex3f ((GLfloat)(-1.0), (GLfloat)(1.0), 0);
		glEnd ();

		glColor3f (0.0, 0.0, 0.0);
		glBegin (GL_LINE_LOOP);
		glVertex3f ((GLfloat)(-1.0), (GLfloat)(-1.0), 0);
		glVertex3f ((GLfloat)(1.0), (GLfloat)(-1.0), 0);
		glVertex3f ((GLfloat)(1.0), (GLfloat)(1.0), 0);
		glVertex3f ((GLfloat)(-1.0), (GLfloat)(1.0), 0);
		glEnd ();
		break;
	}
	case 2: // black triangle
	{
		glColor3f (0.0, 0.0, 0.0);
		glBegin (GL_TRIANGLES);
		glVertex3f ((GLfloat)(-1.0), (GLfloat)(-1.0), 0);
		glVertex3f ((GLfloat)(0.0), (GLfloat)(1.5), 0);
		glVertex3f ((GLfloat)(1.0), (GLfloat)(-1.0), 0);
		glEnd ();
		break;
	}
	case 3: // white triangle
	{
		glColor3f ((GLfloat)1.0, (GLfloat)1.0, (GLfloat)1.0);
		glBegin (GL_TRIANGLES);
		glVertex3f ((GLfloat)(-1.0), (GLfloat)(-1.0), 0);
		glVertex3f ((GLfloat)(0.0), (GLfloat)(1.5), 0);
		glVertex3f ((GLfloat)(1.0), (GLfloat)(-1.0), 0);
		glEnd ();

		glColor3f (0.0, 0.0, 0.0);
		glBegin (GL_LINE_LOOP);
		glVertex3f ((GLfloat)(-1.0), (GLfloat)(-1.0), 0);
		glVertex3f ((GLfloat)(0.0), (GLfloat)(1.5), 0);
		glVertex3f ((GLfloat)(1.0), (GLfloat)(-1.0), 0);
		glEnd ();
		break;
	}
	case 4: // black circle
	{
		glColor3f (0.0, 0.0, 0.0);
		double DEG2RAD = 3.14159 / 180.0;
		glBegin (GL_TRIANGLE_FAN);
		glVertex3f (0.0, 0.0, 0.0);
		for (int i = 0; i < 361; i++)
		{
			double degInRad = i*DEG2RAD;
			glVertex3f ((GLfloat)(cos (degInRad) * 1.0), (GLfloat)(sin (degInRad) * 1.0), 0.0);
		}
		glEnd ();
		break;
	}
	case 5: // white circle
	{
		glColor3f ((GLfloat)1.0, (GLfloat)1.0, (GLfloat)1.0);
		double DEG2RAD = 3.14159 / 180.0;
		glBegin (GL_TRIANGLE_FAN);
		glVertex3f (0.0, 0.0, 0.0);

		for (int i = 0; i < 361; i++)
		{
			double degInRad = i*DEG2RAD;
			glVertex3f ((GLfloat)(cos (degInRad) * 1.0), (GLfloat)(sin (degInRad) * 1.0), 0.0);
		}
		glEnd ();

		glColor3f (0.0, 0.0, 0.0);
		glBegin (GL_LINE_LOOP);
		for (int i = 0; i < 361; i++)
		{
			double degInRad = i*DEG2RAD;
			glVertex3f ((GLfloat)(cos (degInRad) * 1.0), (GLfloat)(sin (degInRad) * 1.0), 0.0);
		}
		glEnd ();
		break;
	}
	default:
	{
		// case 6 - inf
		// stars with different amount of corners
		int angles = (int)((double)counter / 2.0) + 2;

		glRotatef ((GLfloat)(54.0), (GLfloat)(0.0), (GLfloat)(0.0), (GLfloat)(1.0));

		if (counter % 2 == 0)
			glColor3f (0.0, 0.0, 0.0);
		else
			glColor3f ((GLfloat)1.0, (GLfloat)1.0, (GLfloat)1.0);

		double DEG2RAD = 3.14159 / 180.0;
		glBegin (GL_TRIANGLE_FAN);
		glVertex3f (0.0, 0.0, 0.0);
		int step = (int)(360.0 / (double)(2.0 * angles));
		for (int j = 0, i = 0; i < 361; i += step, j++)
		{
			double degInRad = i*DEG2RAD;
			double rad = j % 2 == 0 ? 2.0 : 1.0;
			glVertex3f ((GLfloat)(cos (degInRad) * rad), (GLfloat)(sin (degInRad) * rad), 0.0);
		}
		glEnd ();

		if (counter % 2 != 0)
		{
			glColor3f (0.0, 0.0, 0.0);
			glBegin (GL_LINE_LOOP);
			for (int j = 0, i = 0; i < 361; i += step, j++)
			{
				double degInRad = i*DEG2RAD;
				double rad = j % 2 == 0 ? 2.0 : 1.0;
				glVertex3f ((GLfloat)(cos (degInRad) * rad), (GLfloat)(sin (degInRad) * rad), 0.0);
			}
			glEnd ();
		}
		break;
	}
	}
	glPopMatrix ();
}

int Painter::GetEncoderClsid (const WCHAR* format, CLSID* pClsid)
{
	UINT  num = 0;          // number of image encoders
	UINT  size = 0;         // size of the image encoder array in bytes
	Gdiplus::ImageCodecInfo* pImageCodecInfo = NULL;

	Gdiplus::GetImageEncodersSize (&num, &size);
	if (size == 0)
		return -1;  // Failure

	pImageCodecInfo = (Gdiplus::ImageCodecInfo*)(malloc (size));
	if (pImageCodecInfo == NULL)
		return -1;  // Failure

	GetImageEncoders (num, size, pImageCodecInfo);

	for (UINT j = 0; j < num; ++j)
	{
		if (wcscmp (pImageCodecInfo[j].MimeType, format) == 0)
		{
			*pClsid = pImageCodecInfo[j].Clsid;
			free (pImageCodecInfo);
			return j;  // Success
		}
	}

	free (pImageCodecInfo);
	return -1;  // Failure
}

bool Painter::CaptureScreenShot ()
{
	int size = width * height;
	UINT * pixels = new UINT[width * height];
	memset (pixels, 0, sizeof (UINT)* width * height);

	glFlush (); 
	glFinish ();

	glReadPixels (0, 0, width, height, GL_BGRA_EXT, GL_UNSIGNED_BYTE, pixels);

	if (NULL == pixels)
		return false;

	// Initialize GDI+
	Gdiplus::GdiplusStartupInput gdiplusStartupInput;
	ULONG_PTR gdiplusToken;
	Gdiplus::GdiplusStartup (&gdiplusToken, &gdiplusStartupInput, NULL);

	{
		// Create the dest image
		Gdiplus::Bitmap DestBmp (width, height, PixelFormat32bppARGB);

		Gdiplus::Rect rect1 (0, 0, width, height);

		Gdiplus::BitmapData bitmapData;
		memset (&bitmapData, 0, sizeof (bitmapData));
		DestBmp.LockBits (
			&rect1,
			Gdiplus::ImageLockModeRead,
			PixelFormat32bppARGB,
			&bitmapData);

		int nStride1 = bitmapData.Stride;
		if (nStride1 < 0)
			nStride1 = -nStride1;

		UINT* DestPixels = (UINT*)bitmapData.Scan0;

		if (!DestPixels)
		{
			delete[] pixels;
			return false;
		}

		for (UINT row = 0; row < bitmapData.Height; ++row)
		{
			for (UINT col = 0; col < bitmapData.Width; ++col)
			{
				DestPixels[row * nStride1 / 4 + col] = pixels[row * width + col];
			}
		}

		DestBmp.UnlockBits (
			&bitmapData);

		delete[] pixels;
		pixels = NULL;

		DestBmp.RotateFlip (Gdiplus::RotateNoneFlipY);

		CLSID Clsid;
		int result = GetEncoderClsid (L"image/png", &Clsid);

		if (result < 0)
			return false;

		Gdiplus::Status status = DestBmp.Save (szDestFile.c_str (), &Clsid);
	}
	// Shutdown GDI+
	Gdiplus::GdiplusShutdown (gdiplusToken);

	return true;
}

void Painter::draw_to_file (std::wstring SzDestFile)
{
	// field boundareis
	if ((draw_setup & DRAW_SETUP_PLOT) != DRAW_SETUP_PLOT)
	{
		if (custom_selection == 0)
		{
			if (area_to_draw == -1)
			{
				mesh->get_0_boundaries (c0);
				mesh->get_N_boundaries (cN);
			}
			else
			{
				mesh->get_area_0_boundaries (area_to_draw, c0);
				mesh->get_area_N_boundaries (area_to_draw, cN);
			}
		}
		if (task != NULL)
			if (task->is_symmetrical ())
				c0[0] = -cN[0];
	}

	// resize
	resize ();

	// min max
	if (custom_min_max == 0 &&
		((draw_setup & DRAW_SETUP_FIELD) == DRAW_SETUP_FIELD || (draw_setup & DRAW_SETUP_ISO) == DRAW_SETUP_ISO))
	{
		double * point_coord = new double[2];

		task->get_min_max (k_system, &min_value, &max_value);
		if (task->is_symmetrical () == 2)
		{
			abs (max_value) > abs (min_value) ? min_value = -max_value : max_value = -min_value;
		}

		//min_value = optimization::find_minimum_HJ (task, k_system, point_coord, c0, cN);
		//max_value = optimization::find_maximum_HJ (task, k_system, point_coord, c0, cN);
	}

	if (custom_min_max == 0 &&
		((draw_setup & DRAW_SETUP_DERIVATIVE) == DRAW_SETUP_DERIVATIVE))
	{
		double * point_coord = new double[2];

		min_value = optimization::find_minimum_der_HJ (task, k_system, der_var, point_coord, c0, cN);
		max_value = optimization::find_maximum_der_HJ (task, k_system, der_var, point_coord, c0, cN);
	}

	glutInit (&argc, argv);
	glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize (width, height);

	glutSetOption (GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);

	glutCreateWindow ("Graphics");

	glClearColor (1.0f, 1.0f, 1.0f, 1.0f);
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	currentPainter = this;
	glutDisplayFunc (Painter::display);

	szDestFile = SzDestFile;

	// draw
	glutMainLoop ();
}

void Painter::reset ()
{
	width = height = 450;
	scale_x = scale_y = 1.0;
	k_system = 0;
	min_value = -1e+10;
	max_value = 1e+10;
	custom_min_max = 0;
	painter_font = GLUT_BITMAP_HELVETICA_12;

	pic_counter++;
	wchar_t pic_name[128];
	swprintf (pic_name, 128, L"Pictures//img%i.png", pic_counter);
	szDestFile = pic_name;

	axis_name_x = "";
	axis_name_y = "";
	picture_name = " ";
	legend_blocks = DATA_BLOCK_MIN_MAX;
	draw_setup = 0;
	der_var = -1;

	vect_n_x = vect_n_y = 10;
	anti_gradient = 0;
	min_V_pix_length = 4;
	max_V_pix_length = 8;
	max_V_length = -1e+20;
	min_V_length = 1e+20;

	Frame = 60;
	W_legend = 300;

	custom_scale = 0;
	custom_selection = 0;
	n_x_strokes = n_y_strokes = 1;

	color_scale_iso = COLOR_SCALE_BLACK;
	color_scale_field = COLOR_SCALE_RAINBOW;
	area_to_draw = -1;
	mesh_data = 0;
	N_isovalues = 0;
	mesh_line_width = 1.0;
	iso_line_width = 1.0;
	spec_iso.clear ();
	use_handlers = false;
	flag_plot_rainbow = false;
	figure_size = 5;
	point_size = 1.0;

	plot_line_width = 1.0;
}

void Painter::set_axis_names (char * Axis_name_x, char * Axis_name_y)
{
	axis_name_x = Axis_name_x;
	axis_name_y = Axis_name_y;
}

void Painter::set_picture_name (char * Picture_name)
{
	picture_name = Picture_name;
}

void Painter::no_legend ()
{
	legend_blocks = 0;
	W_legend = 0;
}

void Painter::set_axes_scale (int N_x_strokes, int N_y_strokes)
{
	if (N_x_strokes > 0 && N_y_strokes > 0)
	{
		custom_scale = 1;
		n_x_strokes = N_x_strokes;
		n_y_strokes = N_y_strokes;
	}
}

void Painter::use_bigger_font ()
{
	painter_font = GLUT_BITMAP_HELVETICA_18;
}

void Painter::set_task (Task_pointer * Task)
{
	task = Task;
	mesh = task->get_mesh_pointer ();
}

void Painter::set_mesh (Mesh_Prototype * Mesh)
{
	mesh = Mesh;
}

void Painter::draw_field (int K_system)
{
	k_system = K_system;
	draw_setup = draw_setup | DRAW_SETUP_FIELD;
}

void Painter::draw_field (int K_system, int color_scale)
{
	k_system = K_system;
	draw_setup = draw_setup | DRAW_SETUP_FIELD;
	color_scale_field = color_scale;
}

void Painter::select_section (double * C0, double * CN)
{
	for (int i = 0; i < 2; i++)
	{
		c0[i] = C0[i];
		cN[i] = CN[i];
	}
	custom_selection = 1;
}

void Painter::draw_mesh ()
{
	draw_setup = draw_setup | DRAW_SETUP_MESH;
}

void Painter::draw_derivative (int K_system, int Der_var)
{
	k_system = K_system;
	der_var = Der_var;
	draw_setup = draw_setup | DRAW_SETUP_DERIVATIVE;
}

void Painter::draw_contour_lines (int K_system)
{
	k_system = K_system;
	draw_setup = draw_setup | DRAW_SETUP_ISO;
	N_isovalues = 10;
}

void Painter::draw_contour_lines (int K_system, int n_amount)
{
	k_system = K_system;
	draw_setup = draw_setup | DRAW_SETUP_ISO;
	N_isovalues = n_amount;
}

void Painter::draw_contour_lines (int K_system, int n_amount, int color_scale)
{
	k_system = K_system;
	draw_setup = draw_setup | DRAW_SETUP_ISO;
	color_scale_iso = color_scale;
	N_isovalues = n_amount;
}

void Painter::draw_field_and_contour_lines (int K_system)
{
	k_system = K_system;
	draw_setup = draw_setup | DRAW_SETUP_ISO | DRAW_SETUP_FIELD;
}

void Painter::set_area (int Area_to_draw)
{
	area_to_draw = Area_to_draw;
}

void Painter::set_min_max (double Min_value, double Max_value)
{
	min_value = Min_value;
	max_value = Max_value;
	custom_min_max = 1;
}

void Painter::set_legend_full_iso ()
{
	legend_blocks = legend_blocks | DATA_BLOCK_FULL_ISO_LEGEND;
}

void Painter::add_isoline (double sp_iso_value)
{
	spec_iso.push_back (sp_iso_value);

	draw_setup = draw_setup | DRAW_SETUP_ISO;
	legend_blocks = legend_blocks | DATA_BLOCK_FULL_ISO_LEGEND;
}

void Painter::add_area_dividers ()
{
	draw_setup = draw_setup | DRAW_SETUP_AREA_DIVIDERS;
}

void Painter::set_point_size (double Point_size)
{
	point_size = Point_size;
}

void Painter::add_mesh_markers ()
{
	mesh_data = mesh_data | DRAW_SETUP_MESH_MARKERS;
}

void Painter::draw_gradient_field (int K_system)
{
	draw_setup = draw_setup | DRAW_SETUP_VECTOR_FIELD;
	legend_blocks = legend_blocks | DATA_BLOCK_VECTOR_MIN_MAX;
}

void Painter::draw_antigradient_field (int K_system)
{
	draw_setup = draw_setup | DRAW_SETUP_VECTOR_FIELD;
	anti_gradient = 1;
	legend_blocks = legend_blocks | DATA_BLOCK_VECTOR_MIN_MAX;
}

void Painter::set_vector_net (int Vect_n_x, int Vect_n_y)
{
	vect_n_x = Vect_n_x;
	vect_n_y = Vect_n_y;
}

void Painter::set_plot_line_width (double Plot_line_width)
{
	plot_line_width = Plot_line_width;
}

void Painter::plot_rainbow ()
{
	flag_plot_rainbow = true;
}

void Painter::draw_plot (char * scale_file_name, std::vector<char*> Plot_file_names)
{
}

void Painter::draw_plot (std::vector<char*> Plot_file_names, std::vector<Spline_Pointer*> Splines)
{
	legend_blocks = DATA_BLOCK_PLOT_LEGEND;
	draw_setup = draw_setup | DRAW_SETUP_PLOT;

	splines = Splines;
	std::vector<std::pair<double, double>> current_data;
	double x, y;
	c0[0] = c0[1] = 1e+20;
	cN[0] = cN[1] = -1e+20;

	// read everything into vectors
	for (size_t i = 0; i < Plot_file_names.size (); i++)
	{
		current_data.clear ();
		FILE * file = fopen (Plot_file_names[i], "r");
		while (!feof (file))
		{
			fscanf (file, "%lf %lf", &x, &y);
			// save min/max of x
			if (x < c0[0])
				c0[0] = x;
			if (x > cN[0])
				cN[0] = x;

			// save min/max of y
			if (y < c0[1])
				c0[1] = y;
			if (y > cN[1])
				cN[1] = y;
			current_data.push_back (std::pair <double, double> (x, y));
		}
		plot_data.push_back (current_data);
		fclose (file);
	}

	for (size_t i = 0, i_end = Splines.size (); i < i_end; i++)
	{
		y = optimization::min_Spline_GS (Splines[i], &x);
		if (x < c0[0])
			c0[0] = x;
		if (y < c0[1])
			c0[1] = y;

		y = optimization::max_Spline_GS (Splines[i], &x);
		if (x > cN[0])
			cN[0] = x;
		if (y > cN[1])
			cN[1] = y;
	}
}

void Painter::add_names (std::vector<char*> plot_names)
{
	use_handlers = true;
	handlers.insert (handlers.end(), plot_names.begin (), plot_names.end());
}

void Painter::set_plot_boundaries (double * C0, double * CN)
{
	c0[0] = C0[0];
	c0[1] = C0[1];

	cN[0] = CN[0];
	cN[1] = CN[1];
}

void Painter::set_figure_size (int Figure_size)
{
	figure_size = Figure_size;
}
