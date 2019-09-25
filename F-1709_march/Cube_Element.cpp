#include "Cube_Element.h"

Cube::Cube ()
{
	n_def_nodes = 8;
	n_base_nodes = 8;
	defining_nodes = NULL;
	base_nodes = NULL;
	n_functions = 8;

	dim = 3;
	
	for (int i = 0; i < dim; i++)
	{
		c0[i] = 0;
		cN[i] = 0;
	}
}

Cube::Cube (const Cube & cube)
{
	if (n_def_nodes != cube.n_def_nodes) // if triangle's amount of nodes doesn't equal to amount of nodes of current triangle
	{
		if (defining_nodes != NULL) // if memory was even allocated
			delete[] defining_nodes; // free it

		n_def_nodes = cube.n_def_nodes; // set new amount of nodes
		defining_nodes = new int[n_def_nodes]; // allocate the memory

		if (base_nodes != NULL) // if memory was even allocated
			delete[] base_nodes; // free it

		n_base_nodes = cube.n_base_nodes; // set new amount of nodes
		base_nodes = new int[n_base_nodes]; // allocate the memory
	}

	if (defining_nodes == NULL)
		defining_nodes = new int[n_def_nodes]; // allocate the memory

	for (int i = 0; i < n_def_nodes; i++) // copy nodes
		defining_nodes[i] = cube.defining_nodes[i];

	if (base_nodes == NULL)
		base_nodes = new int[n_base_nodes]; // allocate the memory

	for (int i = 0; i < n_base_nodes; i++) // copy nodes
		base_nodes[i] = cube.base_nodes[i];

	n_functions = cube.n_functions;
	area = cube.area;
	dim = cube.dim;
	for (int i = 0; i < dim; i++)
	{
		c0[i] = cube.c0[i];
		cN[i] = cube.cN[i];
	}

	// set G if it exists
	if (cube.G != NULL)
		G->Full_Copy (*cube.G);
	else
		G = NULL;
	// set M if it exists
	if (cube.M != NULL)
		M->Full_Copy (*cube.M);
	else
		M = NULL;
}

Cube::~Cube ()
{
}

double Cube::get_geometrical_area ()
{
	// volume ?
	double r = 1.0;
	for (int i = 0; i < dim; i++)
		r *= (cN[i] - c0[i]);
	return r;
}

bool Cube::point_inside (const Mesh_Prototype & mesh, double * coordinates)
{
	for (int i = 0; i < dim; i++)
	{
		if (!(coordinates[i] < cN[i] + ZERO_cube && coordinates[i] > c0[i] - ZERO_cube))
		{
			return false;
		}
	}
	return true;
}

void Cube::inside_prepare (const Mesh_Prototype & mesh)
{
	double coordinates[3];
	mesh.get_node_coordinates (base_nodes[0], coordinates);
	for (int i = 0; i < 3; i++)
	{
		c0[i] = coordinates[i];
	}
	mesh.get_node_coordinates (base_nodes[7], coordinates);
	for (int i = 0; i < 3; i++)
	{
		cN[i] = coordinates[i];
	}
}

double Cube::get_basis_function_value (int k_func, double * coordinates)
{
	int f[] = { k_func % 2, (k_func / 2) % 2, (k_func / 4)};

	double V[2];
	double r = 1.0;
	for (int i = 0; i < 3; i++)
	{
		V[0] = (cN[i] - coordinates[i]) / (cN[i] - c0[i]);
		V[1] = (coordinates[i] - c0[i]) / (cN[i] - c0[i]);
		r *= V[f[i]];
	}

	return r;
}

double Cube::get_basis_function_derivative (int k_func, int k_var, double * coordinates)
{
	int f[] = { k_func % 2, (k_func / 2) % 2, (k_func / 4) };

	double V[2];
	double r = 1.0;

	for (int i = 0; i < 3; i++)
	{
		V[0] = (cN[i] - coordinates[i]) / (cN[i] - c0[i]);
		V[1] = (coordinates[i] - c0[i]) / (cN[i] - c0[i]);
		if (i != k_var)
		{
			r *= V[f[i]];
		}
		else
		{
			r *= pow (-1, k_var + 1) / (cN[i] - c0[i]);
		}
	}
	
	return r;
}

void Cube::prepare_GM ()
{
	Matrix * Gx = new Matrix (2, 2);
	Matrix * Mx = new Matrix (2, 2);
	Gx->IdentityMatrix ();
	Gx->setElem (0, 1, -1.0);
	Gx->setElem (1, 0, -1.0);
	Gx->MultiplyByNumber (1.0 / (cN[0] - c0[0]));

	Mx->IdentityMatrix ();
	Mx->MultiplyByNumber (2.0);
	Mx->setElem (0, 1, 1.0);
	Mx->setElem (1, 0, 1.0);
	Mx->MultiplyByNumber (1.0 / 6.0);
	Mx->MultiplyByNumber (cN[0] - c0[0]);

	Matrix * Gy = new Matrix (2, 2);
	Matrix * My = new Matrix (2, 2);
	Gy->IdentityMatrix ();
	Gy->setElem (0, 1, -1.0);
	Gy->setElem (1, 0, -1.0);
	Gy->MultiplyByNumber (1.0 / (cN[1] - c0[1]));

	My->IdentityMatrix ();
	My->MultiplyByNumber (2.0);
	My->setElem (0, 1, 1.0);
	My->setElem (1, 0, 1.0);
	My->MultiplyByNumber (1.0 / 6.0);
	My->MultiplyByNumber (cN[1] - c0[1]);

	Matrix * Gz = new Matrix (2, 2);
	Matrix * Mz = new Matrix (2, 2);
	Gz->IdentityMatrix ();
	Gz->setElem (0, 1, -1.0);
	Gz->setElem (1, 0, -1.0);
	Gz->MultiplyByNumber (1.0 / (cN[2] - c0[2]));

	Mz->IdentityMatrix ();
	Mz->MultiplyByNumber (2.0);
	Mz->setElem (0, 1, 1.0);
	Mz->setElem (1, 0, 1.0);
	Mz->MultiplyByNumber (1.0 / 6.0);
	Mz->MultiplyByNumber (cN[2] - c0[2]);

	int fi[3], fj[3];
	G = new Matrix (n_functions, n_functions);
	M = new Matrix (n_functions, n_functions);

	for (int i = 0; i < 8; i++)
	{
		fi[0] = i % 2;
		fi[1] = (i / 2) % 2;
		fi[2] = (i / 4);

		for (int j = 0; j < 8; j++)
		{
			fj[0] = j % 2;
			fj[1] = (j / 2) % 2;
			fj[2] = (j / 4);

			G->setElem (i, j, Gx->Elem (fi[0], fj[0]) * My->Elem (fi[1], fj[1]) * Mz->Elem (fi[2], fj[2]) +
				Mx->Elem (fi[0], fj[0]) * Gy->Elem (fi[1], fj[1]) * Mz->Elem (fi[2], fj[2]) +
				Mx->Elem (fi[0], fj[0]) * My->Elem (fi[1], fj[1]) * Gz->Elem (fi[2], fj[2]));
			M->setElem (i, j, Mx->Elem (fi[0], fj[0]) * My->Elem (fi[1], fj[1]) * Mz->Elem (fi[2], fj[2]));
		}
	}

	delete Gx;
	delete Mx;
	delete Gy;
	delete My;
	delete Gz;
	delete Mz;
}

int Cube::get_amount_non_zero_functions ()
{
	return 8;
}

void Cube::get_boundaries (double * C0, double * CN)
{
	for (int i = 0; i < 3; i++)
	{
		C0[i] = c0[i];
		CN[i] = cN[i];
	}
}

int Cube::amount_of_integration_points ()
{
	return 14;
}

void Cube::integration_points (double ** points, double * weigths, double * jac)
{
	// Gauss-4 
	int n_int_nodes = 14;
	// set w(i)
	MathVector * w = new MathVector (n_int_nodes);
	{
		w->setElem (0, 320.0 / 361.0);
		w->setElem (1, 320.0 / 361.0);
		w->setElem (2, 320.0 / 361.0);
		w->setElem (3, 320.0 / 361.0);
		w->setElem (4, 320.0 / 361.0);
		w->setElem (5, 320.0 / 361.0);
		w->setElem (6, 121.0 / 361.0);
		w->setElem (7, 121.0 / 361.0);
		w->setElem (8, 121.0 / 361.0);
		w->setElem (9, 121.0 / 361.0);
		w->setElem (10, 121.0 / 361.0);
		w->setElem (11, 121.0 / 361.0);
		w->setElem (12, 121.0 / 361.0);
		w->setElem (13, 121.0 / 361.0);
	}

	// set master point's cordinates	
	MathVector * xi = new MathVector (n_int_nodes);
	MathVector * etta = new MathVector (n_int_nodes);
	MathVector * zeta = new MathVector (n_int_nodes);
	{
		xi->setElem (0, -pow (19.0 / 30.0, 0.5));
		xi->setElem (1, pow (19.0 / 30.0, 0.5));
		xi->setElem (2, 0.0);
		xi->setElem (3, 0.0);
		xi->setElem (4, 0.0);
		xi->setElem (5, 0.0);
		xi->setElem (6, -pow (19.0 / 33.0, 0.5));
		xi->setElem (7, -pow (19.0 / 33.0, 0.5));
		xi->setElem (8, -pow (19.0 / 33.0, 0.5));
		xi->setElem (9, -pow (19.0 / 33.0, 0.5));
		xi->setElem (10, pow (19.0 / 33.0, 0.5));
		xi->setElem (11, pow (19.0 / 33.0, 0.5));
		xi->setElem (12, pow (19.0 / 33.0, 0.5));
		xi->setElem (13, pow (19.0 / 33.0, 0.5));

		etta->setElem (0, 0.0);
		etta->setElem (1, 0.0);
		etta->setElem (2, -pow (19.0 / 30.0, 0.5));
		etta->setElem (3, pow (19.0 / 30.0, 0.5));
		etta->setElem (4, 0.0);
		etta->setElem (5, 0.0);
		etta->setElem (6, -pow (19.0 / 33.0, 0.5));
		etta->setElem (7, -pow (19.0 / 33.0, 0.5));
		etta->setElem (8, pow (19.0 / 33.0, 0.5));
		etta->setElem (9, pow (19.0 / 33.0, 0.5));
		etta->setElem (10, -pow (19.0 / 33.0, 0.5));
		etta->setElem (11, -pow (19.0 / 33.0, 0.5));
		etta->setElem (12, pow (19.0 / 33.0, 0.5));
		etta->setElem (13, pow (19.0 / 33.0, 0.5));

		zeta->setElem (0, 0.0);
		zeta->setElem (1, 0.0);
		zeta->setElem (2, 0.0);
		zeta->setElem (3, 0.0);
		zeta->setElem (4, -pow (19.0 / 30.0, 0.5));
		zeta->setElem (5, pow (19.0 / 30.0, 0.5));
		zeta->setElem (6, -pow (19.0 / 33.0, 0.5));
		zeta->setElem (7, pow (19.0 / 33.0, 0.5));
		zeta->setElem (8, -pow (19.0 / 33.0, 0.5));
		zeta->setElem (9, pow (19.0 / 33.0, 0.5));
		zeta->setElem (10, -pow (19.0 / 33.0, 0.5));
		zeta->setElem (11, pow (19.0 / 33.0, 0.5));
		zeta->setElem (12, -pow (19.0 / 33.0, 0.5));
		zeta->setElem (13, pow (19.0 / 33.0, 0.5));
	}

	for (int k = 0; k < n_int_nodes; k++)
	{
		points[k][0] = (cN[0] - c0[0]) * (xi->getElem (k) + 1.0) / 2.0 + c0[0];
		points[k][1] = (cN[1] - c0[1]) * (etta->getElem (k) + 1.0) / 2.0 + c0[1];
		points[k][2] = (cN[2] - c0[2]) * (zeta->getElem (k) + 1.0) / 2.0 + c0[2];
		weigths[k] = w->getElem (k);
	}
	// jacobian
	*jac = (cN[0] - c0[0]) * (cN[1] - c0[1]) * (cN[2] - c0[2]) / 8.0;

	delete w;
	delete xi;
	delete etta;
	delete zeta;
}