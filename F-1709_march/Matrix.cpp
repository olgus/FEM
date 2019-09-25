#include "Matrix.h"

Matrix::Matrix ()
{
	size0 = size1 = 0;
	elem = NULL;
}

Matrix::Matrix (int n)
{
	size0 = size1 = n;
	elem = new double * [size0];
	for (int i = 0; i < size0; i++)
	{
		elem[i] = new double[size1];
	}
	Zero ();
}

Matrix::Matrix (const Matrix & matrix)
{
	if (size0 != matrix.size0 || size1 != matrix.size1) // if matrixes have different amount of elements, reset it
	{
		for (int i = 0; i < size0; i++)
		{
			delete[] elem[i];
		}
		delete[] elem;

		// resize
		size0 = matrix.size0;
		size1 = matrix.size1;

		// allocate the memory
		elem = new double * [size0];
		for (int i = 0; i < size0; i++)
		{
			elem[i] = new double[size1];
		}
	}

	// copy the elements
	for (int i = 0; i < size0; i++)
	{
		for (int j = 0; j < size1; j++)
		{
			elem[i][j] = matrix.elem[i][j];
		}
	}
}

Matrix::Matrix (int n, int m)
{
	size0 = n;
	size1 = m;
	elem = new double * [size0];
	for (int i = 0; i < size0; i++)
	{
		elem[i] = new double[size1];
	}
	Zero ();
}

Matrix::~Matrix ()
{
	if (elem != NULL)
	{
		for (int i = 0; i < size0; i++)
		{
			delete[] elem[i];
		}
		delete[] elem;
	}
	elem = NULL;
}

void Matrix::Zero ()
{
	for (int i = 0; i < size0; i++)
	{
		for (int j = 0; j < size1; j++)
		{
			elem[i][j] = 0.0L;
		}
	}
}

void Matrix::solve (MathVector * x, const MathVector & f)
{
	Matrix * A = new Matrix (*this);
	MathVector * b = new MathVector (f);

	int i, j, k, kmax;
	double max, s;
	double * t;
	
	// uppertriangle
	for (i = 0; i < size0; i++)
	{
		for (k = i + 1, max = elem[i][i], kmax = i; k < size0; k++)
		{
			if (fabs (max) < fabs (elem[k][i]))
			{
				max = elem[k][i];
				kmax = k;
			}
		}

		if (kmax != i)
		{
			t = elem[i];
			elem[i] = elem[kmax];
			elem[kmax] = t;

			s = b->getElem (i);
			b->setElem (i, b->getElem (kmax));
			b->setElem (kmax, s);
		}

		s = elem[i][i];

		for (j = i; j < size0; j++)
		{
			elem[i][j] /= s;
		}

		b->setElem (i, b->getElem (i) / s);

		for (k = i + 1; k < size0; k++)
		{
			s = elem[k][i];
			for (j = i; j < size0; j++)
			{
				elem[k][j] = elem[k][j] - elem[i][j] * s;
			}
			b->setElem (k, b->getElem (k) - b->getElem (i) * s);
		}
	}
	// retrace

	for (i = x->getSize() - 1; i >= 0; i--)
	{
		x->setElem (i, b->getElem (i) / elem[i][i]);
		for (k = i - 1; k >= 0; k--)
		{
			b->setElem (k, b->getElem (k) - x->getElem (i) * elem[k][i]);
		}
	}

	this->Copy (*A);
	delete A;
	delete b;
}

void Matrix::IdentityMatrix ()
{
	for (int i = 0; i < size0; i++)
	{
		for (int j = 0; j < size1; j++)
		{
			elem[i][j] = 0.0L;
			if (i == j) 
				elem[i][j] = 1.0L;
		}
	}
}

void Matrix::Print ()
{
	for (int i = 0; i < size0; i++)
	{
		for (int j = 0; j < size1; j++)
			printf ("%.7lf  ", elem[i][j]);
		printf ("\n");
	}
	printf ("\n");
}

void Matrix::FPrint (FILE * file)
{
	for (int i = 0; i < size0; i++)
	{
		for (int j = 0; j < size1; j++)
			fprintf (file, "%.16lf  ", elem[i][j]);
		fprintf (file, "\n");
	}
	fprintf (file, "\n");
}

void Matrix::FPrint (char * file_name)
{
	FILE * file = fopen (file_name, "w");
	for (int i = 0; i < size0; i++)
	{
		for (int j = 0; j < size1; j++)
			fprintf (file, "%.16lf  ", elem[i][j]);
		fprintf (file, "\n");
	}
	fclose (file);
}

void Matrix::MultiplyByNumber (double c)
{
	for (int i = 0; i < size0; i++)
		for (int j = 0; j < size1; j++)
			elem[i][j] *= c;
}

void Matrix::Add (const Matrix & m)	//прибавить матрицу к матрице
{
	for (int i = 0; i < size0; i++)
		for (int j = 0; j < size1; j++)
			elem[i][j] += m.elem[i][j];
}

void Matrix::Insert_Row (int nR, const MathVector & row)
{
	for (int j = 0; j < size1; j++)
	{
		elem[nR][j] = row.getElem (j);
	}
}

void Matrix::Get_Row (int nR, MathVector * r)
{
	int vN = r->getSize ();
	for (int j = 0; j < size1 && j < vN; j++)
	{
		r->setElem (j, elem[nR][j]);
	}
}

void Matrix::Cut_Last_Row ()
{
	int old_size = size0;
	int new_size = size0 - 1;

	// allocate the memory to store new_size * size1 elements of the matrix 
	double ** e;
	e = new double * [new_size];
	for (int i = 0; i < new_size; i++)
	{
		e[i] = new double[size1];
		for (int j = 0; j < size1; j++)
		{
			e[i][j] = elem[i][j]; // store elements
		}
	}

	// free this memory
	for (int i = 0; i < old_size; i++)
	{
		delete[] elem[i];
	}
	delete[] elem;

	size0 = new_size; // reset the size0
	elem = new double * [new_size];
	for (int i = 0; i < new_size; i++)
	{
		elem[i] = new double[size1];
		for (int j = 0; j < size1; j++)
		{
			elem[i][j] = e[i][j]; // store elements back
		}
	}

	for (int i = 0; i < new_size; i++)
	{
		delete[] e[i];
	}
	delete[] e;

}

bool Matrix::MultiplyByMatrix (Matrix const & sour, Matrix * dest)
{
	int i, j, k;
	double s;
	if (size1 != sour.size0) 
		return false;

	for (i = 0; i < size0; i++)
	{
		for (j = 0; j < sour.size1; j++)
		{
			s = 0;
			for (k = 0; k < sour.size0; k++)
			{
				s += elem[i][k] * sour.elem[k][j];
			}
			dest->elem[i][j] = s;
		}
	}

	return true;
}

void Matrix::Rotation_Matrix (double c, double s, int start_block)
{
	IdentityMatrix ();

	elem[start_block][start_block] = c;
	elem[start_block + 1][start_block + 1] = c;

	elem[start_block + 1][start_block] = -s;
	elem[start_block][start_block + 1] = s;
}

void Matrix::Copy (Matrix const & source)
{
	for (int i = 0; i < size0 && i < source.size0; i++)
	{
		for (int j = 0; j < size1 && j < source.size1; j++)
		{
			elem[i][j] = source.elem[i][j];
		}
	}
}

void Matrix::Full_Copy (Matrix const & source)
{
	if (size0 != source.size0 || size1 != source.size1) // if matrixes have different amount of elements, reset it
	{
		for (int i = 0; i < size0; i++)
		{
			delete[] elem[i];
		}
		delete[] elem;

		// resize
		size0 = source.size0;
		size1 = source.size1;

		// allocate the memory
		elem = new double *[size0];
		for (int i = 0; i < size0; i++)
		{
			elem[i] = new double[size1];
		}
	}

	// copy the elements
	for (int i = 0; i < size0; i++)
	{
		for (int j = 0; j < size1; j++)
		{
			elem[i][j] = source.elem[i][j];
		}
	}
}

void Matrix::Insert_Column (int nC, MathVector * column)
{
	for (int i = 0; i < size0; i++)
	{
		elem[i][nC] = column->getElem (i);
	}
}

void Matrix::calc_URV (const MathVector & v, MathVector * vect)
{
	int i, k;
	double s;

	MathVector * vcopy = new MathVector (size0);
	vcopy->Copy (v);

	for (i = size0 - 1; i >= 0; i--)
	{
		vect->setElem (i, vcopy->getElem (i) / elem[i][i]);

		for (k = 0; k < i; k++)
		{
			s = elem[k][i] * vect->getElem (i);
			vcopy->addToElem (k, -s);
		}
	}
	delete vcopy;
}

void Matrix::Get_Column (int nC, MathVector * r)
{
	for (int i = 0; i < size0; i++)
	{
		r->setElem(i, elem[i][nC]);
	}
}

void Matrix::MultiplyMatrixByVector (const MathVector & v, MathVector * r)
{
	double s;
	for (int i = 0; i < size0; i++)
	{
		s = 0;
		for (int j = 0; j < size1 && j < v.getSize(); j++)
		{
			s += elem[i][j] * v.getElem (j);
		}
		r->setElem (i, s);
	}
}

int Matrix::Size0 () const
{
	return size0;
}

int Matrix::Size1 () const
{
	return size1;
}

double Matrix::Elem (int i, int j)
{
	return elem[i][j];
}

void Matrix::setElem (int i, int j, double v)
{
	elem[i][j] = v;
}

void Matrix::addToElem (int i, int j, double v)
{
	elem[i][j] += v;
}

bool Matrix::inverse_matrix_size_3 (Matrix * result)
{
	if (size0 != 3 || size0 != 3)
		return false;

	double determinant;

	result->Zero ();
	result->setElem (0, 0, (elem[1][1] * elem[2][2] - elem[1][2] * elem[2][1]));
	result->setElem (0, 1, -(elem[0][1] * elem[2][2] - elem[0][2] * elem[2][1]));
	result->setElem (0, 2, (elem[0][1] * elem[1][2] - elem[0][2] * elem[1][1]));

	result->setElem (1, 0, -(elem[1][0] * elem[2][2] - elem[1][2] * elem[2][0]));
	result->setElem (1, 1, (elem[0][0] * elem[2][2] - elem[0][2] * elem[2][0]));
	result->setElem (1, 2, -(elem[0][0] * elem[1][2] - elem[0][2] * elem[1][0]));

	result->setElem (2, 0, (elem[1][0] * elem[2][1] - elem[1][1] * elem[2][0]));
	result->setElem (2, 1, -(elem[0][0] * elem[2][1] - elem[0][1] * elem[2][0]));
	result->setElem (2, 2, (elem[0][0] * elem[1][1] - elem[0][1] * elem[1][0]));

	determinant = elem[0][0] * result->Elem (0, 0) + elem[0][1] * result->Elem (1, 0) + elem[0][2] * result->Elem (2, 0);
	
	if (fabs (determinant) > ZERO_DETERMINANT)
		result->MultiplyByNumber (1.0 / determinant);
	else
		return false;
	
	return true;
}

double Matrix::get_determinant_3 ()
{
	if (size0 != 3 || size0 != size1)
		return 0.0;

	double a, b, c;
	a = elem[1][1] * elem[2][2] - elem[1][2] * elem[2][1];
	b = -(elem[1][0] * elem[2][2] - elem[1][2] * elem[2][0]);
	c = (elem[1][0] * elem[2][1] - elem[1][1] * elem[2][0]);

	return elem[0][0] * a + elem[0][1] * b + elem[0][2] * c;
}

void Matrix::Transpose ()
{
	// NOT TESTED
	double c;

	for (int i = 1; i < Size0 (); i++)
	{
		for (int j = 0; j < i; j++)
		{
			c = elem[i][j];
			elem[i][j] = elem[j][i];
			elem[j][i] = c;
		}
	}
}
