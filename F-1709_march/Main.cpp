#include "Task_Types.h"
#include "Integration.h"
#include "Free_meshing.h"
#include "MS_comparison.h"
#include "Sparse_Matrixes.h"

#include "Triangulation.h"
#include "Find_theta.h"

#include <string> 

#include <ctime> 

int main (int argc, char *argv[])
{
	{
		//Task_Triangular_Hierarchical pt;
		//pt.prepare ("Source Files//Mesh_lin.txt", "Source Files//Conditions.txt", false);

		//int solver_param[] = {SOLVER_METHOD_PARDISO, SOLVER_DECOMP_TYPE_LU, 5, 3, SOLVER_MKL_NO};
		//pt.solve_task (solver_param);

		//Painter * painter = new Painter ();
		//painter->initialize (argc, argv);
		//painter->set_task (&pt);
		//painter->draw_field (0, COLOR_SCALE_RAINBOW);
		//painter->draw_contour_lines (0, 20);
		//painter->draw_to_file (L"Pictures//test.png");
	}
	{
		//double h[2];
		//h[0] = h[1] = 0.01;
		//int n;

		//Free_mesh mesh;
		//mesh.build_Mesh ("Source Files//Simple//Mesh_free1.txt", "Source Files//Simple//density.txt", "Source Files//Simple//order.txt", "MeshData//Nodes.txt", "MeshData//Triangles.txt");
		////Triangular_Mesh tm;
		////tm.build_Mesh ("Source Files//Mesh_lin.txt", true, &n);

		//Task_Triangular_Hierarchical pr;
		//Cartesian_Test_Task pr;
		//pr.prepare ("MeshData//Nodes.txt", "MeshData//Triangles.txt", "Source Files//conditions.txt");
		////pr.prepare ("Source Files//Mesh_lin.txt", "Source Files//conditions.txt", true);
		//pr.solve_task (0, 1, 30);
		//pr.calc_gradient (0);

		//Painter painter;
		//painter.initialize (argc, argv);
		//painter.set_task (&pr);
		//painter.set_picture_name ("Pic");
		//painter.draw_field (0, COLOR_SCALE_RED_YELLOW);
		//painter.draw_antigradient_field (0);
		//painter.draw_contour_lines (0);
		//painter.set_max_resolution (600);
		//painter.set_scale (1.0, 1.0);
		//painter.add_isoline (100.0);
		//painter.add_isoline (200.0);
		//painter.add_isoline (300.0);
		//painter.draw_to_file (L"Pictures//mesh1.png");


		//pr.save_slice (0, "Source Files//slice2.txt", "Result//slice_y.txt", true);
		//pr.calc_gradient (0);
		//
		//pr.save_slice_first_derivative (0, "Source Files//slice1.txt", "Result//derivative_x.txt", true);
		//pr.save_slice_first_derivative (0, "Source Files//slice2.txt", "Result//derivative_y.txt", true);
		//pr.save_slice_first_derivative (0, "Source Files//slice2.txt", "Result//derivative_y_out.txt", false);
		//pr.calc_gradient (0);
		//pr.save_slice_first_derivative (0, "Source Files//slice2.txt", "Result//derivative_y_sr.txt", false);
		//////double value;
		//////double coord[] = {0.0, 1.9};
		////////////pr.get_solution_in_point (0, coord, &value);
		//Painter painter;
		//painter.initialize (argc, argv);
		//painter.set_task (&pr);
		//painter.set_picture_name ("Mesh");
		//painter.draw_mesh ();
		//painter.draw_to_file (L"Pictures//mesh.png");

		//double a[] = {1e-3, 1e-2};
		//Spline_1D_Smooth sp;
		//sp.prepare ("Result//derivative_y.txt", a);
		//sp.prepare ("Result//slice_y.txt", a);
		//sp.solve_task ();
		//sp.save_slice (0, "Source Files//slice_1D.txt", "Result//slice_1D.txt", false);
		//sp.save_slice_first_derivative (0, "Source Files//slice_1D.txt", "Result//Spline_1D.txt", false);
	}
	{
		//Free_mesh mesh;
		//mesh.build_Mesh ("Source Files//Comparison//Mesh_free2.txt", "Source Files//Comparison//density.txt", "Source Files//Comparison//order.txt", "MeshData//Nodes.txt", "MeshData//Triangles.txt");
		//////Triangular_Mesh tm;
		//////tm.build_Mesh ("Source Files//Mesh_lin.txt", true, &n);

		////Task_Triangular_Hierarchical pr;
		//Cartesian_Test_Task pr;
		//pr.prepare ("MeshData//Nodes.txt", "MeshData//Triangles.txt", "Source Files//conditions.txt");
		//////pr.prepare ("Source Files//Mesh_lin.txt", "Source Files//conditions.txt", true);
		//pr.solve_task (0, 30);

		//Painter painter;
		//painter.initialize (argc, argv);
		//painter.set_task (&pr);
		//painter.set_picture_name ("PicPicPic");
		//painter.set_axes_scale (5, 5);
		//painter.draw_gradient_field (0);
		//painter.set_axis_names ("x", "Vy");
		//painter.draw_field_and_contour_lines (0);
		//painter.draw_to_file (L"Pictures//mesh2.png");
	}
	{
		//Free_mesh mesh;
		//mesh.build_Mesh ("Source Files//Comparison//Mesh_free3.txt", "Source Files//Comparison//density.txt", "Source Files//Comparison//order.txt", "MeshData//Nodes.txt", "MeshData//Triangles.txt");
		//////Triangular_Mesh tm;
		//////tm.build_Mesh ("Source Files//Mesh_lin.txt", true, &n);

		////Task_Triangular_Hierarchical pr;
		//Cartesian_Test_Task pr;
		//pr.prepare ("MeshData//Nodes.txt", "MeshData//Triangles.txt", "Source Files//conditions.txt");
		//////pr.prepare ("Source Files//Mesh_lin.txt", "Source Files//conditions.txt", true);
		//pr.solve_task (0, 30);

		//Painter painter;
		//painter.initialize (argc, argv);
		//painter.set_task (&pr);
		//painter.set_picture_name ("PicPicPicPicPicPic");
		//painter.draw_mesh ();
		//painter.draw_field (0, COLOR_SCALE_BLACK_WHITE);
		//painter.set_axes_scale (2, 20);
		//painter.draw_gradient_field (0);
		//painter.set_iso_line_width (3.0);
		//painter.set_axis_names ("X", "Y");
		//painter.draw_contour_lines (0, 20, COLOR_SCALE_RAINBOW);
		//painter.set_legend_full_iso ();
		//painter.add_isoline (100.0);
		//painter.add_isoline (200.0);
		//painter.add_isoline (300.0);
		//painter.add_isoline (400.0);
		//painter.draw_to_file (L"Pictures//mesh3.png");
	}
	{
		//Task_Vector_2D pr;
		//pr.prepare ("MeshData//Nodes.txt", "MeshData//Triangles.txt", "Source Files//conditions.txt");
		////pr.prepare ("Source Files//Mesh_lin.txt", "Source Files//conditions.txt", true);
		//pr.solve_task (0, 10);
		//pr.fprint_edges ("Result Files Extra//edges.txt");
		//pr.fprint_functions_by_elements ("Result Files Extra//functions.txt");
		//double value[2];
		//double coord[] = { 0.0, 0.0 };
		//pr.get_solution_in_point (0, coord, value);
		//printf ("%lf %lf\n", value[0], value[1]);
		//coord[0] = coord[1] = 0.5;
		//pr.get_solution_in_point (0, coord, value);
		//printf ("%lf %lf\n", value[0], value[1]);
		//coord[1] = 0.8;
		//pr.get_solution_in_point (0, coord, value);
		//printf ("%lf %lf\n", value[0], value[1]);
		//coord[0] = coord[1] = 1.0;
		//pr.get_solution_in_point (0, coord, value);
		//printf ("%lf %lf\n", value[0], value[1]);
		//pr.save_evenly_2D (0, h, "Result//slice_even.txt", false);

		//Painter painter1;
		//painter = &painter1;
		//painter->set_painter (&pr, 1000, 0);
		//painter->Initialize (argc, argv);
		//painter->draw_mesh (L"Pictures//mesh.png");
		//double h_vf = 0.1;
		//painter->draw_vector_field (0, h_vf, L"Pictures//v_field.png");
		//painter = NULL;
	}
	{
	//	double extra[2];
	//	extra[0] = 10.0;
	//	extra[1] = 1e-10;
	//	 
	//	Spline_1D_Smooth sp1;
	//	sp1.prepare ("Result//TC//water_bottom_cryst//det//detgraphs//cell_Vx_1_na_0.07300.txt", extra, 1.0 / 1.7);
	//	sp1.solve_task ();
	//	extra[0] = 10.0;
	//	extra[1] = 1e-10;
	//	Spline_1D_Smooth sp2;
	//	sp2.prepare ("Result//TC//water_bottom_cryst//det//detgraphs//cell_Vx_5_na_0.07300.txt", extra, 1.0 / 1.7);
	//	sp2.solve_task ();
	//	Spline_1D_Smooth sp3;
	//	sp3.prepare ("Result//TC//water_bottom_cryst//det//detgraphs//cell_Vx_8_na_0.07300.txt", extra, 1.0 / 1.7);
	//	sp3.solve_task ();
	//	//Spline_1D_Smooth sp4;
	//	//sp4.prepare ("Result//TC//water_bottom_cryst//det//detgraphs//T_vert_0.12000.txt", extra, 1.0 / 2.5);
	//	//sp4.solve_task ();
	//	//Spline_1D_Smooth sp5;
	//	//sp5.prepare ("Result//TC//water_bottom_cryst//det//detgraphs//T_vert_0.14000.txt", extra, 1.0 / 2.5);
	//	//sp5.solve_task ();
	//	//Spline_1D_Smooth sp6;
	//	//sp6.prepare ("Result//TC//water_bottom_cryst//det//detgraphs//T_vert_0.16000.txt", extra, 1.0 / 2.5);
	//	//sp6.solve_task ();
	//	//Spline_1D_Smooth sp7;
	//	//sp7.prepare ("Result//TC//water_bottom_cryst//det//detgraphs//T_vert_0.17500.txt", extra, 1.0 / 2.5);
	//	//sp7.solve_task ();

	//	std::vector<char *> plots = { /*"Result//TC//water_bottom_cryst//det//detgraphs//cell_T1_0.07300.txt",
	//		"Result//TC//water_bottom_cryst//det//detgraphs//cell_T5_0.07300.txt",
	//		"Result//TC//water_bottom_cryst//det//detgraphs//cell_T8_0.07300.txt"
	//		"Result//Melt_SF//temp//data8.txt",
	//		"Result//Melt_SF//temp//data10.txt", "Result//Melt_SF//temp//data11.txt"*/};
	//	//std::vector<char *> plots = { };
	//	std::vector<Spline_Pointer *> splines;
	//	splines.push_back (&sp1);
	//	splines.push_back (&sp2);
	//	splines.push_back (&sp3);
	//	//splines.push_back (&sp4);
	//	//splines.push_back (&sp5);
	//	//splines.push_back (&sp6);
	//	//splines.push_back (&sp7);
	//	////
	//	std::vector<char *> names = { "49.35", "53.55", "56.7", "738", "861", "984", "1076.25" };
	//	Painter painter;
	//	painter.initialize (argc, argv);
	//	//painter.plot_rainbow ();
	//	painter.add_names (names);
	//	painter.set_figure_size (4);
	//	painter.set_max_resolution (700);
	//	painter.set_plot_line_width (3.0);
	//	painter.set_plot_line_width (1.0);
	//	painter.set_scale (1.0, 10.0);
	//	painter.set_axes_scale (14, 13);
	//	painter.draw_plot (plots, splines);
	//	painter.set_axis_names ("Y", "|Vx|");
	//	double c0[] = {0.0, -0.8};
	//	double cN[] = {35.0, 0.5};
	//	painter.set_plot_boundaries (c0, cN);
	//	painter.set_picture_name (" ");
	//	painter.draw_to_file (L"Pictures//TC//water_bottom_cryst//graphs//Vx_na_cell_hor_cu_bw.png");

	/////*	FILE * file_q_Table = fopen ("Result//Melt_SF//q_table.txt", "w");

	////	double ht = 0.05;
	////	double t = 0.05;
	////	double value;
	////	double c[] = { 0.0 };
	////	char name[128];
	////	
	////	for (int i = 1; t < 10.0; i++)
	////	{
	////		t = ht * i;
	////		Spline_1D_Smooth sp1;
	////		sprintf (name, "Result//Melt_SF//50 mm//hfb_t_%0.5lf.txt", t);
	////		sp1.prepare (name, extra);
	////		sp1.solve_task ();

	////		sp1.get_solution_in_point (c, &value);
	////		fprintf (file_q_Table, "%0.2lf %.16lf\n", t , value);
	////	}

	////	fclose (file_q_Table);*/
	}
	{
		//Free_mesh mesh;
		//mesh.build_Mesh ("Source Files//Iron//Mesh_free.txt", "Source Files//Iron//density.txt", "Source Files//Iron//order.txt", "MeshData//Nodes.txt", "MeshData//Triangles.txt");
		//Cartesian_Test_Task pt;
		//pt.prepare ("MeshData//Nodes.txt", "MeshData//Triangles.txt", "Source Files//Iron//conditions.txt");
		//pt.solve_task (0, 80);
		//pt.save_solution_in_points (0, "Source Files//Iron//field_points.txt", "Result//solution.txt", false);
		//
		//Painter painter1;
		//painter = &painter1;
		//painter->set_painter (&pt, 1000, 0);
		//painter->Initialize (argc, argv);
		//painter->draw_field_and_isovalues (0, 20, true, false, L"Pictures//field.png");
		//painter->draw_mesh (false, L"Pictures//mesh.png");
		//painter = NULL; 
	}
	{
		//compressed_matrix crm;
		//MathVector res (5);
		//crm.get_from_files ();
		//crm.factorize_Cholesky ();
		//crm.fprint ();
	}
	{
		//double c0[2], cN[2], h[2]; 
		//c0[0] = c0[1] = -1.0;
		//cN[0] = cN[1] = 1.0; 
		//h[0] = h[1] = 0.1;
		//Spline * spline = new Spline ();
		//spline->prepare (c0, cN, 0.0, "Source Files//Comparison//density.txt");
		//spline->solve_task (-1, 200);
		//spline->save_evenly_2D (0, h, "Result//spline.txt", true);
	}
	{ 
	//	//Free_mesh mesh;
	//	//mesh.build_Mesh_wo_renumerate ("Source Files//Melt_SF//Mesh_free.txt", "Source Files//Melt_SF//density.txt", "Source Files//Melt_SF//order.txt", "MeshData//Nodes.txt", "MeshData//Triangles.txt");
	//	
	//	Sleep (2000);

	//	Task_Melt_SF<Triangular_Mesh> pt ("Source Files//Melt_SF//Parameters.txt");
	//	Task_pointer * tp = &pt;
	//	char * bound[] = { "Source Files//Melt_SF//Conditions_T.txt",
	//		"Source Files//Melt_SF//Conditions_P.txt",
	//		"Source Files//Melt_SF//Conditions_W.txt" };
	//	pt.prepare ("MeshData//Nodes.txt", "MeshData//Triangles.txt", bound, "Source Files//Melt_SF//time_mapping.txt", "Source Files//Melt_SF//time_stamps.txt");


	//	Painter * painter = new Painter;
	//	painter->initialize (argc, argv);
	//	painter->set_task (tp);
	//	painter->set_max_resolution (1000);
	//	painter->draw_mesh ();
	//	painter->draw_to_file (L"Pictures//Melt_SF//mesh.png");
	//	painter->reset ();
	//	pt.painter_pointer (painter);

	//	//pt.use_density_inversion ();
	//	pt.get_scales ("Source Files//Melt_SF//Scales.txt");
	//	int param[][5] = { { SOLVER_METHOD_GMRES, SOLVER_DECOMP_TYPE_LU, 30, 4, SOLVER_MKL_NO},
	//	{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 5, 4, SOLVER_MKL_NO },
	//	{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 5, 4, SOLVER_MKL_NO } };

	//	double GR, Ste, PC, phys_T0_F;
	//	GR = 570.491136260751;
	//	Ste = 0.00492537313432836;
	//	PC = 22.0 - 21.0;
	//	phys_T0_F = 23.0 - 21.0;
	//	double theta = 11.0;
	//	pt.reset_parameters (theta, GR * theta, Ste * theta, PC / theta, phys_T0_F / theta);
	//	pt.solve_task (param);

	//	FILE * file = fopen ("Result//Melt_SF//final_front_position.txt", "w");
	//	double x[] = { 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 
	//		0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 
	//		1.1, 1.2, 1.3, 1.4, 1.5, 2.0};

	//	for (int i = 0; i < 27; i++)
	//	{
	//		fprintf (file, "%.16lf %.16lf\n", x[i], pt.get_front_position (x[i]));
	//	}
	//	fclose (file);
	}
	{
		////Free_mesh mesh;
		////mesh.build_Mesh_wo_renumerate ("Source Files//MESH//Mesh_free.txt", "Source Files//MESH//density.txt", "Source Files//MESH//order.txt", "Source Files//MESH//Nodes.txt", "Source Files//MESH//Triangles.txt");
		////mesh.build_Mesh ("Source Files//Melt_SF//Mesh_free.txt", "Source Files//Melt_SF//density.txt", "Source Files//Melt_SF//order.txt", "MeshData//Nodes.txt", "MeshData//Triangles.txt");
		//  
		////Sleep (2000);
		//{
		//	Triangular_Mesh mesh;
		//	double c0[] = { 0.0, 0.0 };
		//	double cN[] = {1000.0, 1000.0};
		//	int n_axis[] = {2, 3};
		//	mesh.build_mesh (c0, cN, n_axis, "Source Files//MESH//Nodes.txt", "Source Files//MESH//Triangles.txt");

		//}
		//Task_Melt_SF<Triangular_Mesh> pt;
		//Task_pointer * tp = &pt;
		//char * bound[] = { "Source Files//Melt_SF//Conditions_T.txt",
		//	"Source Files//Melt_SF//Conditions_P.txt",
		//	"Source Files//Melt_SF//Conditions_W.txt" };
		//pt.prepare ("Source Files//MESH//Nodes.txt", "Source Files//MESH//Triangles.txt", bound, "Source Files//Melt_SF//time_mapping.txt", "Source Files//Melt_SF//time_stamps.txt");
		////pt.get_scales ();

		//Painter * painter = new Painter;
		//painter->initialize (argc, argv);
		//painter->set_task (tp);
		//painter->set_max_resolution (1000);
		//painter->draw_mesh ();
		//painter->draw_to_file (L"Pictures//mesh.png");
		//painter->reset ();
		//painter->initialize (argc, argv);
		//painter->set_task (tp);
		//pt.painter_pointer (painter);

		//pt.solve_task (0, 1, 20);
	}
	{
	//	//Free_mesh mesh_rel;
	//	//mesh_rel.build_Mesh ("Source Files//Melt_SF//res_ult//Mesh_free.txt", "Source Files//Melt_SF//res_ult//density.txt", "Source Files//Melt_SF//res_ult//order.txt", "MeshData//Nodes_rel.txt", "MeshData//Triangles_rel.txt");

	//	Task_Melt_SF<Mixed_Triangular_Mesh> pt_rel;
	//	char * bound[] = { "Source Files//Melt_SF//Conditions_T.txt",
	//		"Source Files//Melt_SF//Conditions_P.txt",
	//		"Source Files//Melt_SF//Conditions_W.txt" };

	//	pt_rel.prepare ("MeshData//Nodes_rel.txt", "MeshData//Triangles_rel.txt", bound, "Source Files//Melt_SF//time_mapping.txt", "Source Files//Melt_SF//time_stamps.txt");

	//	//Free_mesh mesh_comp;
	//	//mesh_comp.build_Mesh ("Source Files//Melt_SF//comp//Mesh_free.txt", "Source Files//Melt_SF//comp//density.txt", "Source Files//Melt_SF//comp//order.txt", "MeshData//Nodes.txt", "MeshData//Triangles.txt");

	//	Task_Melt_SF<Mixed_Triangular_Mesh> pt_comp;
	//	 
	//	pt_comp.prepare ("MeshData//Nodes.txt", "MeshData//Triangles.txt", bound, "Source Files//Melt_SF//time_mapping.txt", "Source Files//Melt_SF//time_stamps.txt");
	//
	//	char * sol_rel[] = { "Result//Melt_SF//res_ult//s0_t_0.20000.txt",
	//		"Result//Melt_SF//res_ult//s1_t_0.20000.txt",
	//		"Result//Melt_SF//res_ult//s2_t_0.20000.txt",
	//		"Result//Melt_SF//res_ult//s0_t_0.40000.txt",
	//		"Result//Melt_SF//res_ult//s1_t_0.40000.txt",
	//		"Result//Melt_SF//res_ult//s2_t_0.40000.txt",
	//		"Result//Melt_SF//res_ult//s0_t_0.60000.txt",
	//		"Result//Melt_SF//res_ult//s1_t_0.60000.txt",
	//		"Result//Melt_SF//res_ult//s2_t_0.60000.txt",
	//		"Result//Melt_SF//res_ult//s0_t_2.00000.txt",
	//		"Result//Melt_SF//res_ult//s1_t_2.00000.txt",
	//		"Result//Melt_SF//res_ult//s2_t_2.00000.txt" };
	//	char * sol_comp[] = { "Result//Melt_SF//res_3//s0_t_0.20000.txt",
	//		"Result//Melt_SF//res_3//s1_t_0.20000.txt",
	//		"Result//Melt_SF//res_3//s2_t_0.20000.txt",
	//		"Result//Melt_SF//res_3//s0_t_0.40000.txt",
	//		"Result//Melt_SF//res_3//s1_t_0.40000.txt",
	//		"Result//Melt_SF//res_3//s2_t_0.40000.txt",
	//		"Result//Melt_SF//res_3//s0_t_0.60000.txt",
	//		"Result//Melt_SF//res_3//s1_t_0.60000.txt",
	//		"Result//Melt_SF//res_3//s2_t_0.60000.txt",
	//		"Result//Melt_SF//res_3//s0_t_2.00000.txt",
	//		"Result//Melt_SF//res_3//s1_t_2.00000.txt",
	//		"Result//Melt_SF//res_3//s2_t_2.00000.txt" };

	//	FILE * file = fopen ("Result//Melt_SF//Comp_result.txt", "w");
	//	for (int j = 0; j < 4; j++)
	//	{
	//		double r;
	//		for (int i = 0; i < 3; i++)
	//		{
	//			pt_rel.read_solution (i, sol_rel[j+i]);
	//			pt_comp.read_solution (i, sol_comp[j + i]);
	//			r = compare_solutions (i, &pt_comp, &pt_rel);
	//			printf ("%e\n", r);
	//			fprintf (file, "%e\t", r);
	//		}
	//		fprintf (file, "\n");
	//		printf ("\n");
	//	}
	//	fclose (file);
	}

	{
		//// full mesh for the fem task
		////Free_mesh mesh;
		////mesh.set_alpha_beta (1e-5, 1e-10);
		////mesh.build_Mesh_wo_renumerate ("Source Files//2d_source_fem//mesh.txt", "Source Files//2d_source_fem//density.txt", "Source Files//2d_source_fem//order.txt", "Source Files//2d_source_fem//Nodes.txt", "Source Files//2d_source_fem//Triangles.txt");
		//int n_f;

		////// tr mesh
		////Triangular_Mesh tm;
		////tm.build_Mesh ("Source Files//2d_source_fem//Mesh_nf_2.txt", "Source Files//2d_source_fem//Nodes.txt", "Source Files//2d_source_fem//Triangles.txt", &n_f);
		////// fem task
		////Task_2D_point_source task_1;
		////task_1.prepare ("Source Files//2d_source_fem//Nodes.txt", "Source Files//2d_source_fem//Triangles.txt", "Source Files//2d_source_fem//conditions.txt", "Source Files//2d_source_fem//point_sources.txt");
		////task_1.get_lambda_values ("Source Files//2d_source_fem//lambda.txt");
		////int solver_param[] = {SOLVER_METHOD_GMRES, SOLVER_DECOMP_TYPE_D, 20, 4, SOLVER_MKL_NO};

		//// field selection task
		//int n_insertions = 2;
		//char * mesh_insertion[] = { "Source Files//2d_source_fem//mesh_insertion_nf_1.txt", "Source Files//2d_source_fem//mesh_insertion_nf_2.txt" };
		//char * mesh_insertion_density[] = { "Source Files//2d_source_fem//density_insertion_1.txt", "Source Files//2d_source_fem//density_insertion_2.txt" };
		//char * mesh_insertion_nodes[] = { "Source Files//2d_source_fem//nodes_insertion_1.txt", "Source Files//2d_source_fem//nodes_insertion_2.txt" };
		//char * mesh_insertion_triangles[] = { "Source Files//2d_source_fem//triangles_insertion_1.txt" , "Source Files//2d_source_fem//triangles_insertion_2.txt" };
		//for (int i = 0; i < n_insertions; i++)
		//{
		//	//Triangular_Mesh tm;
		//	//tm.build_Mesh (mesh_insertion[i], mesh_insertion_nodes[i], mesh_insertion_triangles[i], &n_f);
		//	Painter * painter = new Painter;
		//	painter->initialize (argc, argv);
		//	Cartesian_Test_Task task_empty;
		//	task_empty.prepare (mesh_insertion_nodes[i], mesh_insertion_triangles[i], "Source Files//2d_source_fem//conditions.txt");
		//	Task_pointer * task_empty_pointer = &task_empty;
		//	painter->set_task (task_empty_pointer);
		//	painter->set_max_resolution (800);
		//	painter->draw_mesh ();
		//	double c0[] = { -300.0, -300.0 };
		//	double cN[] = { 500.0, 300.0 };
		//	//double c0[] = { -1000.0, -1000.0 };
		//	//double cN[] = { 1000.0, 1000.0 };
		//	painter->select_section (c0, cN);
		//	std::wstring file_name = L"Pictures//2d_source//mesh";
		//	file_name += std::to_wstring (i);
		//	file_name += L".png";
		//	painter->draw_to_file (file_name);
		//	delete painter;
		//}
		////task_1.solve_task (solver_param);

		//Task_2D_source_w_field_selection task_2;
		//task_2.prepare (n_insertions, "Source Files//2d_source_fem//point_sources.txt", "Source Files//2d_source_fem//lambda.txt", mesh_insertion_nodes, mesh_insertion_triangles);
		//task_2.solve_task ();


		////Task_pointer * tp = &task_1;
		//Task_pointer * tp_f = &task_2;
		////task_1.eq_matrixes[0].f->FPrint ("Result//2d_source_fem//f.txt");
		//Task_2D_Mesh_Interpolation task_interp;
		//{
		//	//double c0[] = { -300.0, 150.0 };
		//	//double cN[] = { -150.0, 300.0 };
		//	double c0[] = { 300.0, -250.0 };
		//	double cN[] = { 500.0, -150.0 };
		//	int N_axis[] = { 100, 100 };
		//	task_interp.prepare (tp_f, 0, c0, cN, N_axis);
		//}
		//task_interp.solve_task (2, 0, 0);
		//Task_pointer * tp2 = &task_interp;

		//// painter init
		//Painter * painter = new Painter;
		//painter->initialize (argc, argv);
		////painter->set_task (tp);

		////// mesh for full fem task
		//// mesh
		////painter->set_max_resolution (1000);
		////painter->draw_mesh ();
		////painter->draw_to_file (L"Pictures//2d_source//mesh_fem.png");
		////painter->reset ();
		////// mesh close-up
		////{
		////	double c0[] = { -10.0, -250.0 };
		////	double cN[] = { 500.0, 10.0 };
		////	painter->select_section (c0, cN);
		////}
		////painter->set_max_resolution (800);
		////painter->draw_mesh ();
		////painter->draw_to_file (L"Pictures//2d_source//mesh_closeup.png");
		////painter->reset ();


		//// field everywhere
		////painter->set_max_resolution (1000);
		//////painter->set_min_max (0.0, 1e-5);
		////painter->draw_field (0, COLOR_SCALE_RAINBOW);
		////painter->draw_contour_lines (0, 20);
		////painter->draw_to_file (L"Pictures//2d_source//field_full.png");
		////painter->reset ();

		//// field close up 1
		////painter->set_max_resolution (800);
		////painter->draw_field (0, COLOR_SCALE_RAINBOW);
		////painter->draw_contour_lines (0, 30);
		////{
		////	//double c0[] = { 300.0, -250.0 };
		////	//double cN[] = { 500.0, -150.0 };
		////	double c0[] = { -300.0, 150.0 };
		////	double cN[] = { -150.0, 300.0 };
		////	painter->select_section (c0, cN);
		////}
		////painter->draw_to_file (L"Pictures//2d_source//field_close1.png");
		////painter->reset ();
		//// field close up 3
		////painter->set_max_resolution (1000);
		////painter->draw_field (0, COLOR_SCALE_RAINBOW);
		////painter->draw_contour_lines (0, 100);
		////{
		////	double c0[] = { 300.0, -250.0 };
		////	double cN[] = { 500.0, -150.0 };
		////	//double c0[] = { -500.0, 400.0 };
		////	//double cN[] = { -400.0, 500.0 };
		////	painter->select_section (c0, cN);
		////}
		////painter->draw_to_file (L"Pictures//2d_source//field_close_3.png");
		////painter->reset ();
		////
		//painter->set_task (tp2);
		////painter->set_max_resolution (800);
		////painter->draw_field (0, COLOR_SCALE_RAINBOW);
		////painter->draw_contour_lines (0, 30);
		////{
		////	double c0[] = { 300.0, -250.0 };
		////	double cN[] = { 500.0, -150.0 };
		////	painter->select_section (c0, cN);
		////}
		////painter->draw_to_file (L"Pictures//2d_source//field_close2.png");
		////painter->reset ();
		////// field close up 2
		////painter->set_max_resolution (800);
		////painter->draw_field (0, COLOR_SCALE_RAINBOW);
		////painter->draw_contour_lines (0, 50);
		////{
		////	double c0[] = {4.0, -70.0 };
		////	double cN[] = { 110.0, -10.0 };
		////	painter->select_section (c0, cN);
		////}
		////painter->draw_to_file (L"Pictures//2d_source//field_close2.png");
		////painter->reset ();


		//// field close up 3
		//painter->set_max_resolution (1000);
		//painter->draw_field (0, COLOR_SCALE_RAINBOW);
		//painter->draw_contour_lines (0, 100);
		//{
		//	double c0[] = { 300.0, -250.0 };
		//	double cN[] = { 500.0, -150.0 };
		//	//double c0[] = { -300.0, 150.0 };
		//	//double cN[] = { -150.0, 300.0 };
		//	painter->select_section (c0, cN);
		//}
		//painter->draw_to_file (L"Pictures//2d_source//field2.png");
		//painter->reset ();

		////// field close up 4
		////painter->set_max_resolution (800);
		////painter->draw_field (0, COLOR_SCALE_RAINBOW);
		////painter->draw_contour_lines (0, 100);
		////{
		////	double c0[] = { 15.0, -110.0 };
		////	double cN[] = { 110.0, 110.0 };
		////	painter->select_section (c0, cN);
		////}
		////painter->draw_to_file (L"Pictures//2d_source//field3.png");
		////painter->reset ();

		//// save slices
		//double bound[] = { 100.0, 300.0 };
		////task_1.save_slice (0, bound, "Source Files//2d_source_fem//slice_vert_1.txt", "Result//2d_source_fem//slice_vert_1.txt", false);
		//task_2.save_slice (0, bound, "Source Files//2d_source_fem//slice_vert_1.txt", "Result//2d_source_fem//slice_vert_1_fs.txt", false);
		//double bound2[] = {0.0, 700.0};
		////task_1.save_slice (0, bound2, "Source Files//2d_source_fem//slice_vert_1.txt", "Result//2d_source_fem//slice_vert_1_full.txt", false);
		//task_2.save_slice (0, bound2, "Source Files//2d_source_fem//slice_vert_1.txt", "Result//2d_source_fem//slice_vert_1_fs_full.txt", false);
		//double bound3[] = { -300.0, -100.0 };
		////task_1.save_slice (0, bound3, "Source Files//2d_source_fem//slice_vert_2.txt", "Result//2d_source_fem//slice_vert_2.txt", false);
		//task_2.save_slice (0, bound3, "Source Files//2d_source_fem//slice_vert_2.txt", "Result//2d_source_fem//slice_vert_2_fs.txt", false);
		//double bound4[] = { -700.0, 0.0 };
		////task_1.save_slice (0, bound4, "Source Files//2d_source_fem//slice_vert_2.txt", "Result//2d_source_fem//slice_vert_2_full.txt", false);
		//task_2.save_slice (0, bound4, "Source Files//2d_source_fem//slice_vert_2.txt", "Result//2d_source_fem//slice_vert_2_fs_full.txt", false);
	}

	// solver
	{	
		//FILE * log = fopen ("D://Projects//F-1709//F-1709//Result//solver_log.txt", "a");

		//char * folder = "D://Projects//F-1709//F-1709//Result//Melt_SF//Matrix//2";
		//compressed_matrix crm;
		//crm.get_from_files (folder);

		//MathVector *v = new MathVector (crm.Size ());
		//char start[256];
		//strcpy (start, folder);
		//strcat (start, "//sv.txt");
		//FILE * file = fopen (start, "r");
		//v->getFromFile (file);
		//fclose (file);

		//crm.convert_to_CSR ();
		//int K_REPEATS = 10;
		//for (int i = 1; i < 9; i++)
		//{
		//	double ms = 0.0;
		//	for (int k = 0; k < K_REPEATS; k++)
		//	{
		//		crm.set_starting_point (v);
		//		crm.Inner_copy ();

		//		clock_t t1, t2;
		//		t1 = clock ();

		//		crm.solve_LOS_MKL (1, i);

		//		t2 = clock ();
		//		ms += (double)(t2 - t1);

		//	}
		//	printf ("time: %lf s\n", ms / 1000 / K_REPEATS);
		//	fprintf (log, "2 LOS LU MKL time: %lf s, threads: %i \n", ms / 1000 / K_REPEATS, i);
		//}
		//crm.get_solution (v);
		//char result[256];
		//strcpy (result, folder);
		//strcat (result, "//x.txt");	
		//v->FPrint (result);
		//delete v;
		//fclose (log);
	}
	{
		////Free_mesh mesh;
		////mesh.build_Mesh ("Source Files//Melt_SF//Mesh_free.txt", "Source Files//Melt_SF//density.txt", "Source Files//Melt_SF//order.txt", "MeshData//Nodes.txt", "MeshData//Triangles.txt");

		//Painter * painter = new Painter;
		//painter->initialize (argc, argv);
		//
		//double alpha = 1e-7;
		//for (int i = 0; i < 10; i++)
		//{
		//	Task_Find_Theta tft;
		//	tft.solve_ ("Source Files//Melt_SF//RT.txt", painter, alpha);
		//	alpha = alpha / 10.0;
		//}
	}
	{
		//{
		//	NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay tmd;
		//	tmd.build_mesh ("Source Files//Triangulation//PG_data.txt", "Source Files//Triangulation//nodes.txt", "Source Files//Triangulation//elements.txt");
		//	Mesh_Prototype * mesh_pointer = &tmd;
		//	Cartesian_Test_Task pt;
		//	pt.prepare ("Source Files//Triangulation//nodes.txt", "Source Files//Triangulation//elements.txt", "Source Files//Conditions.txt");
		//	pt.solve_task (SOLVER_METHOD_LU, SOLVER_DECOMP_TYPE_D, 0);
		//	Task_pointer * tp = &pt;

		//	Painter * painter = new Painter;
		//	painter->initialize (argc, argv);
		//	painter->set_task (tp);
		//	painter->set_max_resolution (800);
		//	painter->set_point_size (5.0);
		//	painter->draw_mesh ();
		//	painter->draw_to_file (L"Source Files//Triangulation//mesh.png");
		//	painter->reset ();

		//	painter->set_task (tp);
		//	painter->set_max_resolution (800);
		//	painter->draw_contour_lines (0, 10, COLOR_SCALE_BLACK);
		//	painter->draw_field (0, COLOR_SCALE_RAINBOW);
		//	painter->draw_to_file (L"Source Files//Triangulation//pic.png");
		//	delete painter;
		//}
	}
	{
		{
			//Painter * painter = new Painter;
			//painter->initialize (argc, argv);

			////NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay tmd;
			////tmd.build_mesh ("Source Files//Melt_SF//PG_data.txt", "Source Files//Melt_SF//nodes.txt", "Source Files//Melt_SF//elements.txt");
			////Mesh_Prototype * mesh_pointer = &tmd;
			////painter->set_mesh (mesh_pointer);
			////painter->set_max_resolution (1000);
			////painter->set_point_size (0.5);
			////painter->draw_mesh ();
			////painter->draw_to_file (L"Pictures//Melt_SF//mesh.png");
			// 
			//Sleep (2000);

			//painter->reset ();
			//Task_Melt_SF<Triangular_Mesh> pt ("Source Files//Melt_SF//Parameters.txt");
			//Task_pointer * tp = &pt;
			//char * bound[] = { "Source Files//Melt_SF//Conditions_T.txt",
			//	"Source Files//Melt_SF//Conditions_P.txt",
			//	"Source Files//Melt_SF//Conditions_W.txt" };
			//pt.prepare ("Source Files//Melt_SF//nodes.txt", "Source Files//Melt_SF//elements.txt", bound, "Source Files//Melt_SF//time_mapping.txt", "Source Files//Melt_SF//time_stamps.txt");

			//pt.use_density_inversion ();
			//pt.get_scales ("Source Files//Melt_SF//Scales.txt");
			//pt.painter_pointer (painter);

			//painter->set_task (tp);
			//painter->set_max_resolution (1000);
			//painter->set_point_size (0.5);
			//painter->draw_mesh ();
			//painter->draw_to_file (L"Pictures//Melt_SF//mesh.png");
			//painter->reset ();

			//char * solutions[] = { "Result//Melt_SF//cryst_centered//s0_t_0.27030.txt",
			//	"Result//Melt_SF//cryst_centered//s1_t_0.27030.txt",
			//	"Result//Melt_SF//cryst_centered//s2_t_0.27030.txt" };
			//pt.start_from_selected_time_layer (0.27030, solutions);

			//int param[][5] = { { SOLVER_METHOD_GMRES, SOLVER_DECOMP_TYPE_D, 30, 4, SOLVER_MKL_NO},
			//{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 5, 4, SOLVER_MKL_NO },
			//{ SOLVER_METHOD_LOS, SOLVER_DECOMP_TYPE_LU, 5, 4, SOLVER_MKL_NO } };

			//pt.solve_task (param);
		}
	}
	{
		//// mesh
		//Painter * painter = new Painter;
		//painter->initialize (argc, argv);

		////NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay tmd;
		////tmd.build_mesh ("Source Files//NLST//PG_data.txt", "Source Files//NLST//nodes.txt", "Source Files//NLST//elements.txt");
		////Mesh_Prototype * mesh_pointer = &tmd;
		////painter->set_mesh (mesh_pointer);
		////painter->set_max_resolution (1000);
		////painter->set_point_size (0.5);
		////painter->draw_mesh ();
		////painter->draw_to_file (L"Pictures//NLST//mesh.png");

		//Sleep (2000);

		//painter->reset ();

		//Non_linear_simple_task nlst;
		//nlst.prepare ("Source Files//NLST//nodes.txt", "Source Files//NLST//elements.txt", "Source Files//NLST//bound.txt");
		//nlst.painter_pointer (painter);

		//Task_pointer * tp = &nlst;
		//painter->set_task (tp);
		//painter->set_max_resolution (1000);
		//painter->reset ();

		//int param[][5] = { { SOLVER_METHOD_GMRES, SOLVER_DECOMP_TYPE_D, 30, 4, SOLVER_MKL_NO } };

		//nlst.solve_task (param);
	}
	{
		//// mesh
		//Painter * painter = new Painter;
		//painter->initialize (argc, argv);

		////NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay tmd;
		////tmd.build_mesh ("Source Files//NLNSTT//PG_data.txt", "Source Files//NLNSTT//nodes.txt", "Source Files//NLNSTT//elements.txt");
		////Mesh_Prototype * mesh_pointer = &tmd;
		////painter->set_mesh (mesh_pointer);
		////painter->set_max_resolution (1000);
		////painter->set_point_size (0.5);
		////painter->draw_mesh ();
		////painter->draw_to_file (L"Pictures//NLNSTT//mesh.png");

		//Sleep (2000);

		//painter->reset ();

		//Nonlinear_nonstationary_task nlst;
		//nlst.prepare ("Source Files//NLNSTT//nodes.txt", "Source Files//NLNSTT//elements.txt", "Source Files//NLNSTT//bound.txt", "Source Files//NLNSTT//time_layers.txt");
		//nlst.painter_pointer (painter);

		//Task_pointer * tp = &nlst;
		//painter->set_task (tp);
		//painter->set_max_resolution (1000);
		//painter->reset ();

		//int param[][5] = { { SOLVER_METHOD_PARDISO, SOLVER_DECOMP_TYPE_LU, 30, 4, SOLVER_MKL_YES } };

		//double ms;
		//clock_t t1, t2;
		//t1 = clock ();

		//nlst.solve_task (param);

		//t2 = clock ();
		//ms = (double)(t2 - t1);
		//printf ("time: %lf s\n", ms / 1000);
	}

	{
		//int param[][5] = { { SOLVER_METHOD_PARDISO, SOLVER_DECOMP_TYPE_LU, 30, 4, SOLVER_MKL_YES } };

		//double ms;

		//compressed_matrix CM;
		//CM.get_from_files ("Test//Matrix");
		//MathVector v, res;
		//v.setSize (CM.Size ());
		//res.setSize (CM.Size ());
		//for (int i = 0; i < CM.Size (); i++)
		//	v.setElem (i, i * 0.5 + 1.0);

		//int par_num_threads = 1;
		//double * res_omp = NULL;
		//if (par_num_threads > 1)
		//	res_omp = new double[(par_num_threads - 1) * CM.Size ()];

		//clock_t t1, t2;
		//t1 = clock ();

		//CM.mult_A_v (v, &res, res_omp, par_num_threads);

		//t2 = clock ();
		//ms = (double)(t2 - t1);
		//printf ("time: %lf s\n", ms / 1000);
	}

	{
		//// mesh
		//Painter * painter = new Painter;
		//painter->initialize (argc, argv);

		////NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay tmd;
		////tmd.build_mesh ("Source Files//NLNSSTT//PG_data.txt", "Source Files//NLNSSTT//nodes.txt", "Source Files//NLNSSTT//elements.txt");
		////Mesh_Prototype * mesh_pointer = &tmd;
		////painter->set_mesh (mesh_pointer); 
		////painter->set_max_resolution (1000);
		////painter->set_point_size (0.5);
		////painter->draw_mesh ();
		////painter->draw_to_file (L"Pictures//NLNSSTT//mesh.png");

		//Sleep (2000);

		//painter->reset ();

		//Nonlinear_nonstationary_system_task nlsst;
		//char * bound[] = { "Source Files//NLNSSTT//cond_0.txt",
		//	"Source Files//NLNSSTT//cond_1.txt",
		//	"Source Files//NLNSSTT//cond_2.txt" };
		//nlsst.prepare ("Source Files//NLNSSTT//nodes.txt", "Source Files//NLNSSTT//elements.txt", bound, "Source Files//NLNSSTT//time_layers.txt");
		//nlsst.painter_pointer (painter);

		//Task_pointer * tp = &nlsst;
		//painter->set_task (tp);
		//painter->set_max_resolution (1000);
		//painter->reset ();

		//int param[][5] = { { SOLVER_METHOD_PARDISO, SOLVER_DECOMP_TYPE_LU, 30, 1, SOLVER_MKL_NO }, 
		//{ SOLVER_METHOD_PARDISO, SOLVER_DECOMP_TYPE_LU, 30, 1, SOLVER_MKL_NO } ,
		//{ SOLVER_METHOD_PARDISO, SOLVER_DECOMP_TYPE_LU, 30, 1, SOLVER_MKL_NO } };

		//double ms;
		//clock_t t1, t2;
		//t1 = clock ();

		//nlsst.solve_task (param);

		//t2 = clock ();
		//ms = (double)(t2 - t1);
		//printf ("time: %lf s\n", ms / 1000);
	}


	{
		//// mesh
		//Painter * painter = new Painter;
		//painter->initialize (argc, argv);

		////NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay tmd;
		////tmd.build_mesh ("Source Files//TC//PG_data.txt", "Source Files//TC//nodes.txt", "Source Files//TC//elements.txt");
		////Mesh_Prototype * mesh_pointer = &tmd;
		////painter->set_mesh (mesh_pointer); 
		////painter->set_max_resolution (1000);
		////painter->set_point_size (0.5);
		////painter->draw_mesh ();
		////painter->draw_to_file (L"Pictures//TC//mesh.png");
		// 
		//Sleep (2000);

		//painter->reset ();

		//Task_convec tc ("Source Files//TC//param.txt");
		//char * bound[] = { "Source Files//TC//cond_0.txt",
		//	"Source Files//TC//cond_1.txt",
		//	"Source Files//TC//cond_2.txt" };
		//tc.prepare ("Source Files//TC//nodes.txt", "Source Files//TC//elements.txt", bound, "Source Files//TC//time_layers.txt");
		//tc.painter_pointer (painter);

		//Task_pointer * tp = &tc;
		//painter->set_task (tp);
		//painter->set_max_resolution (900);
		//painter->reset ();

		//int param[][5] = { { SOLVER_METHOD_PARDISO, SOLVER_DECOMP_TYPE_LU, 30, 1, SOLVER_MKL_NO },
		//{ SOLVER_METHOD_PARDISO, SOLVER_DECOMP_TYPE_LU, 30, 1, SOLVER_MKL_NO } ,
		//{ SOLVER_METHOD_PARDISO, SOLVER_DECOMP_TYPE_LU, 30, 1, SOLVER_MKL_NO } };

		//char * sol[] = { "Result//TC//s0_t_0.07000.txt",
		//	"Result//TC//s1_t_0.07000.txt",
		//	"Result//TC//s2_t_0.07000.txt",
		//};

		//tc.start_from_selected_time_layer (0.07000, sol);
		//double ms;
		//clock_t t1, t2;
		//t1 = clock ();

		//tc.solve_task (param);

		//t2 = clock ();
		//ms = (double)(t2 - t1);
		//printf ("time: %lf s\n", ms / 1000);
	}

	{
		// mesh
		Painter * painter = new Painter;
		painter->initialize (argc, argv);

		//NS_Triangular_Mesh_Delaunay::Triangular_Mesh_Delaunay tmd;
		//tmd.build_mesh ("Source Files//TC//PG_data.txt", "Source Files//TC//nodes.txt", "Source Files//TC//elements.txt");
		//Mesh_Prototype * mesh_pointer = &tmd;
		//painter->set_mesh (mesh_pointer); 
		//painter->set_max_resolution (1000);
		//painter->set_point_size (0.5);
		//painter->draw_mesh ();
		//painter->draw_to_file (L"Pictures//TC//mesh.png");

		Sleep (2000);

		painter->reset ();

		Task_convec tc ("Source Files//TC//param.txt");
		char * bound[] = { "Source Files//TC//cond_0.txt",
			"Source Files//TC//cond_1.txt",
			"Source Files//TC//cond_2.txt" };
		tc.prepare ("Source Files//TC//nodes.txt", "Source Files//TC//elements.txt", bound, "Source Files//TC//time_layers.txt");
		tc.painter_pointer (painter);

		Task_pointer * tp = &tc;
		painter->set_task (tp);
		painter->set_max_resolution (900);
		painter->reset ();

		int param[][5] = { { SOLVER_METHOD_PARDISO, SOLVER_DECOMP_TYPE_LU, 30, 1, SOLVER_MKL_NO },
		{ SOLVER_METHOD_PARDISO, SOLVER_DECOMP_TYPE_LU, 30, 1, SOLVER_MKL_NO } ,
		{ SOLVER_METHOD_PARDISO, SOLVER_DECOMP_TYPE_LU, 30, 1, SOLVER_MKL_NO } };

		char * sol[3];
		for (int i = 0; i < 3; i++)
			sol[i] = new char[128];

		double time_layer;
		double st_time_layer = 0.14000;

		for (int k = 0; k < 1; k++)
		{
			printf ("%i ", k);

			time_layer = st_time_layer + k * 0.001;
			for (int i = 0; i < 3; i++)
			{
				sol[i] = new char[128];
				sprintf (sol[i], "Result//TC//water_bottom_cryst//s%i_t_%.5lf.txt", i, time_layer);
			}

			tc.start_from_selected_time_layer (time_layer, sol);

			painter->set_min_max (-8.0, 8.0);
			painter->set_max_resolution (800);
			painter->draw_contour_lines (2, 12);
			painter->draw_field (2, COLOR_SCALE_RAINBOW);
			painter->draw_to_file (L"Pictures//TC//pic2.png");
			//char file_name[128];
			//double bound[] = {0.08, 0.2};

			//sprintf (file_name, "Result//TC//water_bottom_cryst//det//cell_Vx_1_na_%.5lf.txt", time_layer);
			//double w[] = {1.0, 1.0};
			//tc.save_slice_first_derivative (2, 1, "Source Files//TC//slice_cell.txt", file_name, true);
			//tc.save_slice (0, "Source Files//TC//slice_T_vert.txt", file_name, true);
			//tc.save_front_pos (file_name);
		}
	}
}