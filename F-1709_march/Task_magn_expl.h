#pragma once
#include "Mesh_Types.h"
#include "Sparse_Matrixes.h"

#define MAX_STEPS 1000

class Task_magn_expl 
{
private:
	int n_areas;
	double ** p; // P by areas

	int n_received;
	double ** b_received; // receivers data
	double gamma_start;
	double FUNCTIONAL_THRESHOLD_PERCENTILE;
	double DES_DIF;
	double P_ZERO;

	// system
	Matrix * sA;
	MathVector * sB;
	MathVector * sP;
	Matrix * L;

	// regularization
	MathVector * Gamma;

	// receivers data
	int rec_n_z_levels; // amount of levels by z
	double * rec_z_levels; // their values
						   
	int rec_n_x_levels; // amount of x receiver levels
	double * rec_x_y_levels; // y components of the receivers
	double rec_x0, rec_xN; // boundaries 
	int rec_n_x_by_level; // amount on the level

	int rec_n_y_levels; // amount of y receiver levels
	double * rec_y_x_levels; // x components of the receivers
	double rec_y0, rec_yN; // boundaries 
	int rec_n_y_by_level; // amount on the level

	void initialize ();
	void build_system ();
	void output ();
	void regularize ();
	int adapt (); // increases gammas
	double calc_functional (); // functional with regularization (is that correct)
	void initialize_l ();
public:
	Cubic_Mesh mesh; // mesh

	Task_magn_expl (); 
	Task_magn_expl (const Task_magn_expl & task);
	~Task_magn_expl ();

	void set_parameters (double Gamma_start, double threshold, double des_dif, double p_zero);
	void prepare_field (char * mesh_file, char * file_areas, char * file_area_data); // mesh file, file with points and values, file with areas	
	void set_receivers_data (char * file_receivers_data); // set receivers
	void get_receivers_data (char * file_receivers_data); // gets data from them

	void save_field (char * file_result, bool calc_p); // saves field in requested points that is calculated through formula
	double receive (int i_receiver, bool calc); // calculates the field in specific point
	double solve ();
	double get_P (int k_element, int k_component);
	double get_P_true (int k_element, int k_component);
}; 