#include "Prism.h"

Prism::Prism ()
{
	n_def_nodes = 6;
	n_base_nodes = 6;
	defining_nodes = NULL;
	base_nodes = NULL;
	Alpha = NULL;
	D = NULL;
	n_functions = 6;
	determinant = 0.0;
	z0 = 0.0;
	zN = 0.0;

	dim = 3;
}

Prism::Prism (const Prism & prism)
{
	if (n_def_nodes != prism.n_def_nodes) // if triangle's amount of nodes doesn't equal to amount of nodes of current triangle
	{
		if (defining_nodes != NULL) // if memory was even allocated
			delete[] defining_nodes; // free it

		n_def_nodes = prism.n_def_nodes; // set new amount of nodes
		defining_nodes = new int[n_def_nodes]; // allocate the memory

		if (base_nodes != NULL) // if memory was even allocated
			delete[] base_nodes; // free it

		n_base_nodes = prism.n_base_nodes; // set new amount of nodes
		base_nodes = new int[n_base_nodes]; // allocate the memory
	}

	if (defining_nodes == NULL)
		defining_nodes = new int[n_def_nodes]; // allocate the memory

	for (int i = 0; i < n_def_nodes; i++) // copy nodes
		defining_nodes[i] = prism.defining_nodes[i];

	if (base_nodes == NULL)
		base_nodes = new int[n_base_nodes]; // allocate the memory

	for (int i = 0; i < n_base_nodes; i++) // copy nodes
		base_nodes[i] = prism.base_nodes[i];

	n_functions = prism.n_functions;
	area = prism.area;
	z0 = prism.z0;
	zN = prism.zN;
	determinant = prism.determinant;
	dim = prism.dim;

	// set Alpha if it exists
	if (prism.Alpha != NULL)
		Alpha->Full_Copy (*prism.Alpha);
	else
		Alpha = NULL;
	// set G if it exists
	if (prism.G != NULL)
		G->Full_Copy (*prism.G);
	else
		G = NULL;
	// set M if it exists
	if (prism.M != NULL)
		M->Full_Copy (*prism.M);
	else
		M = NULL;
	// set D if it exists
	if (prism.D != NULL)
		D->Full_Copy (*prism.D);
	else
		D = NULL;
}

Prism::~Prism ()
{
	if (Alpha != NULL)
		delete Alpha;
	if (D != NULL)
		delete D;
}

double Prism::get_geometrical_area ()
{
	// volume ?
	return fabs (determinant) * (zN - z0) / 2.0;
}

bool Prism::point_inside (const Mesh_Prototype & mesh, double * coordinates)
{
	// check by z
	if (coordinates[2] < zN + ZERO_prism && coordinates[2] > z0 - ZERO_prism)
	{
		// check by the triangle
		double tr_coord[3][3];
		// get triangle's coordinates
		for (int i = 0; i < 3; i++)
		{
			mesh.get_node_coordinates (base_nodes[i], tr_coord[i]);
		}

		// get triangle's geometrical area
		double tr_S = fabs (determinant) / 2.0;

		// calculate geometrical areas of triangles made with point (coordinates)
		double square = 0.0;
		for (int i = 0; i < 3; i++)
		{
			square += fabs ((tr_coord[(i + 1) % 3][0] - tr_coord[i][0]) * (coordinates[1] - tr_coord[i][1]) -
				(coordinates[0] - tr_coord[i][0]) * (tr_coord[(i + 1) % 3][1] - tr_coord[i][1]));
		}
		square /= 2.0;

		// if half of their sum is equal to element's geometrical area
		if (abs (square - tr_S) < ZERO_element)
		{
			// point is in the element
			return true;
		}
	}
	return false;
}

void Prism::inside_prepare (const Mesh_Prototype & mesh)
{
	Alpha = new Matrix (3, 3);
	D = new Matrix (3, 3);
	double coordinates[3];

	for (int i = 0; i < 3; i++)
	{
		mesh.get_node_coordinates (base_nodes[i], coordinates);
		D->setElem (0, i, 1.0);
		D->setElem (1, i, coordinates[0]);
		D->setElem (2, i, coordinates[1]);
	}
	z0 = coordinates[2];
	mesh.get_node_coordinates (base_nodes[3], coordinates);
	zN = coordinates[2];

	determinant = D->get_determinant_3 ();

	D->inverse_matrix_size_3 (Alpha);
}

double Prism::get_basis_function_value (int k_func, double * coordinates)
{
	double r = 0.0;

	// get function value on triangle
	int tr_func = k_func % 3;
	r += Alpha->Elem (tr_func, 0);
	r += Alpha->Elem (tr_func, 1) * coordinates[0];
	r += Alpha->Elem (tr_func, 2) * coordinates[1];

	// get Z-function value
	int z_func = k_func / 3;
	if (z_func == 0) 
	{
		r *= (zN - coordinates[2]) / (zN - z0);
	}
	else
	{
		r *= (coordinates[2] - z0) / (zN - z0);
	}
	return r;
}

double Prism::get_basis_function_derivative (int k_func, int k_var, double * coordinates)
{
	double r = 0.0;

	int z_func = k_func / 3;
	int tr_func = k_func % 3;

	if (k_var == 2)
	{	
		// get function value on triangle
		r += Alpha->Elem (tr_func, 0);
		r += Alpha->Elem (tr_func, 1) * coordinates[0];
		r += Alpha->Elem (tr_func, 2) * coordinates[1];
		r *= pow (-1.0 , z_func + 1);
		r /= (zN - z0);
	}
	else
	{
		r += Alpha->Elem (tr_func, k_var + 1);
		// get Z-function value
		if (z_func == 0)
		{
			r *= (zN - coordinates[2]) / (zN - z0);
		}
		else
		{
			r *= (coordinates[2] - z0) / (zN - z0);
		}
	}
	return r;
}

void Prism::prepare_GM ()
{
	Matrix * Gxy = new Matrix (3, 3);
	Matrix * Mxy = new Matrix (3, 3);

	Matrix * Gz = new Matrix (2, 2);
	Matrix * Mz = new Matrix (2, 2);

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Gxy->setElem (i, j, (Alpha->Elem (i, 1) * Alpha->Elem (j, 1) + Alpha->Elem (i, 2) * Alpha->Elem (j, 2)));
			Mxy->setElem (i, j, 1.0);
			if (i == j)
			{
				Mxy->setElem (i, j, 2.0);
			}
		}
	}
	Gxy->MultiplyByNumber (fabs (determinant) / 2.0);
	Mxy->MultiplyByNumber (fabs (determinant) / 24.0);

	Gz->IdentityMatrix ();
	Gz->setElem (0, 1, -1.0);
	Gz->setElem (1, 0, -1.0);
	Gz->MultiplyByNumber (1.0 / (zN - z0));

	Mz->IdentityMatrix ();
	Mz->MultiplyByNumber (2.0);
	Mz->setElem (0, 1, 1.0);
	Mz->setElem (1, 0, 1.0);
	Mz->MultiplyByNumber ((zN - z0) / 6.0);

	G = new Matrix (n_functions, n_functions);
	M = new Matrix (n_functions, n_functions);
	int mi, ni, mj, nj;
	for (int i = 0; i < 6; i++)
	{
		mi = i % 3;
		ni = i / 3;
		for (int j = 0; j < 6; j++)
		{
			mj = j % 3;
			nj = j / 3;

			G->setElem (i, j, Gxy->Elem (mi, mj) * Mz->Elem (ni, nj) + Mxy->Elem (mi, mj) * Gz->Elem (ni, nj));
			M->setElem (i, j, Mxy->Elem (mi, mj) * Mz->Elem (ni, nj));
		}
	}

	delete Gxy;
	delete Mxy;
	delete Gz;
	delete Mz;
}

int Prism::get_amount_non_zero_functions ()
{
	return n_functions;
}

void Prism::get_D (int var_der, Matrix * D)
{
	double r;
	double val;
	int mj, vi, vj, mi;

	switch (var_der)
	{
	case 0:
	case 1:	
		r = fabs (determinant);
		r *= (zN - z0);
		r /= 36.0;
		for (int i = 0, i_end = D->Size0 (); i < i_end; i++)
		{
			vi = i / 3;
			for (int j = 0, j_end = D->Size1 (); j < j_end; j++)
			{
				mj = j % 3;
				vj = j / 3;
				val = Alpha->Elem (mj, var_der + 1);
				if (vi == vj)
					val *= 2.0;
				D->setElem (i, j, r * val);
			}
		}
		break;
	case 2:
		r = fabs (determinant);
		r /= 48.0;
		for (int i = 0, i_end = D->Size0 (); i < i_end; i++)
		{
			mi = i % 3;
			for (int j = 0, j_end = D->Size1 (); j < j_end; j++)
			{
				mj = j % 3;
				vj = j / 3;
				val = pow (-1, vj + 1);
				if (mi == mj)
					val *= 2.0;
				D->setElem (i, j, r * val);
			}
		}
		break;
	}
}
