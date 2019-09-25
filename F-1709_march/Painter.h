#pragma once
#include <stdio.h>
#include <Windows.h>
#include "GL\glut.h"
#include "GL\freeglut.h"
#include <vector>
#include <gdiplus.h>
#include <windows.h>  
#include "Painter_def.h"

#include "Mesh_Pointer.h"
#include "Task_Pointer.h"
#include "Spline_1D_Pointer.h"

#include "Optimization.h"
#include "Triangular_Mesh.h"
#include "Task_magn_expl.h"

struct IsovalueSection
{
public:
	int triangle;
	int isovalue;
	Point_2D point[2];
};

struct Block
{
public:
	double c0[2], cN[2];
};

/*
painter that draws:
+ 2D solution 
- mesh (colored by order, colored by area)
- solution as field, b/w, color scale rainbox, color scale red-yellow
- axis data included
- min/max fixed or found as an option
- legend included (min/max)
- solution as contour lines, b/w, color scale rainbox, color scale red-yellow
- optional vertical symmetry 
- add contour lines with fixed isovalues
- draw only a subsection by area (rescale) 
- derivative

+ vector field
- for gradient of the task

+ 1D graphs
- from multiple files (different scales)
- x from one file and multiple y files
- axis data and legend are included
*/

class Painter
{
private:

	static Painter * currentPainter;
	int argc;
	char ** argv;
	int pic_counter;
	int k_system;

	// universal variables
	int width, height; // size of the picture itself
	int pic_width, pic_height; // size of the picture itself
	double scale_x, scale_y; // scaling options
	std::wstring szDestFile; // file to write in
	char * axis_name_x, * axis_name_y; // axis names
	char * picture_name; // name of the picture 
	double c0[2], cN[2]; // field boundaries
	void * painter_font;

	int legend_blocks; // specifies which blocks are to put in legend section
	int draw_setup; // specifies what to draw
	Block picture; // picture block
	Block legend; // legend block

	int Frame; // width of a frame
	int W_legend; // width of a legend 

	int custom_scale; // custom scale for axes strokes
	int custom_selection;
	int n_x_strokes, n_y_strokes; // amount of strokes

	// 2D solution section
	Task_pointer * task; // task to draw
	Mesh_Prototype * mesh; // mesh to draw on
	int color_scale_iso, color_scale_field; // color scales
	int area_to_draw; // number of the area to draw, if -1 then everything
	int mesh_data; // defines extra stuff for mesh data
	int N_isovalues; // amount of contour lines to draw
	double mesh_line_width; // line width for mesh
	double iso_line_width; // line width for contour lines
	int custom_min_max; // search min/max flag
	double min_value, max_value; // min and max values of the solution
	int der_var; // derivative variable, -1 for solution
	std::vector<double> spec_iso;
	int figure_size; // size of the figure for plots
	double point_size; // for nodes

	// vector field
	int vect_n_x, vect_n_y;
	int anti_gradient;
	double min_V_pix_length, max_V_pix_length;
	double max_V_length, min_V_length;

	// 1D plot
	double plot_line_width; // line width for plot
	std::vector <std::vector<std::pair<double, double>>> plot_data; // vector with plot points
	std::vector<Spline_Pointer *> splines; // spline plots
	std::vector <char *> handlers;
	bool use_handlers;
	bool flag_plot_rainbow;
protected:
	void resize (); // sets width and height according to picture
	void draw_blocks (); // display functions by blocks
	void draw_block_field (); // field
	void draw_block_legend ();// legend

	void draw_block_vector_field (); // vector field
	void draw_block_iso (); // isolines
	void draw_block_white_frame (); // adds frame to cover pieces of triangles
	void draw_block_mesh (); // mesh
	void draw_block_plot (); // plots
	void draw_block_areas (); // plots
	// data
	void draw_block_frame (); // axes, frame
	void draw_name (); // adds a name of the picture
	void draw_block_strokes (); // adds strokes and their values
	void draw_block_stroke_lines (); // adds net for plot
	void scale_axes (double pxmax, double pymax, double pxmin, double pymin, double * x_axis_step, double * y_axis_step);

	// color scales
	void get_color (int setup, double value, double * rgb);
	void get_color (int number, double * rgb);

	void add_figure (double * point, double * reverse_scale, int counter);

	// file section
	int Painter::GetEncoderClsid (const WCHAR* format, CLSID* pClsid);
	bool Painter::CaptureScreenShot ();
public:

	Painter ();
	~Painter ();

	static void display (); // display function for opengl

	void initialize (int argc, char * argv[]);
	void draw_to_file (std::wstring SzDestFile); // draw into file current picture 
	void reset (); // dump all the options

	void set_max_resolution (int max_res); // set width/ height
	void set_scale (double Scale_x, double Scale_y); // set scale
	void set_mesh_line_width (double Mesh_line_width); // set mesh line width
	void set_iso_line_width (double Iso_line_width); // set contour line width
	void set_axis_names (char * Axis_name_x, char * Axis_name_y); // set axis names
	void set_picture_name (char * Picture_name); // set axis names
	void no_legend (); // no legend
	void set_axes_scale (int N_x_strokes, int N_y_strokes); // sets custom amount of strokes for axes
	void use_bigger_font ();

	// 2D solution section
	void set_task (Task_pointer * Task); // set task
	void set_mesh (Mesh_Prototype * Mesh); // set mesh, if different from task's mesh
	void draw_field (int K_system); // setup = draw field
	void draw_field (int K_system, int color_scale); // setup = draw field with specific color scale
	void select_section (double * C0, double * CN); // select specific section on the field to draw
	void draw_mesh (); // setup = draw mesh
	void draw_derivative (int K_system, int Der_var); // setup = draw derivative, still needs field and iso options!
	void draw_contour_lines (int K_system); // setup = draw contour lines
	void draw_contour_lines (int K_system, int n_amount); // setup = draw contour lines
	void draw_contour_lines (int K_system, int n_amount, int color_scale); // setup = draw contour lines
	void draw_field_and_contour_lines (int K_system); // setup = field and contour lines
	void set_area (int Area_to_draw); // set specific area (material) to draw solution there
	void set_min_max (double Min_value, double Max_value); // set custom min max
	void set_legend_full_iso ();
	void add_isoline (double sp_iso_value);
	void add_area_dividers ();
	void set_point_size (double Point_size);

	// vector field
	void draw_gradient_field (int K_system); // draw gradient field
	void draw_antigradient_field (int K_system); // draw antigradient field
	void set_vector_net (int Vect_n_x, int Vect_n_y); // sets amounts of vectors to display

	// 1D plot
	void set_plot_line_width (double Plot_line_width); // set line width
	void plot_rainbow ();
	void draw_plot (char * scale_file_name, std::vector<char *> Plot_file_names); // draw plots for single scale in different file
	void draw_plot (std::vector<char *> Plot_file_names, std::vector<Spline_Pointer *> Splines); // draw plots from multiple files with possible different scales
	void add_names (std::vector<char *> plot_names); // draw plots from multiple files with possible different scales
	void set_plot_boundaries (double * C0, double * CN);
	void set_figure_size (int Figure_size);
};