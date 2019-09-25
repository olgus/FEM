#include "MathVector.h"

MathVector::MathVector ()
{
	size = 0;
	elem = NULL;
}

MathVector::MathVector (int n)
{
	size = n;
	elem = new double[size];
	Zero ();
}

MathVector::MathVector (const MathVector & mvector)
{
	if (size == 0) // if vector has not been initialized
	{
		size = mvector.size; // set size
		if (size > 0)
			elem = new double[size]; // allocate the memory
	}
	if (size != mvector.size) // if vector has been initialized, but with different size
	{
		delete[] elem; // free allocated memory
		size = mvector.size; // reset size
		if (size > 0)
			elem = new double[size]; // allocate the memory
	}
	for (int i = 0; i < size; i++) // copy elements
	{
		elem[i] = mvector.elem[i];
	}
}

void MathVector::Zero ()
{
	for (int i = 0; i < size; i++)
	{
		elem[i] = 0.0L;
	}
}

void MathVector::Initialize (double * e) 
{
	for (int i = 0; i < size; i++)
	{
		elem[i] = e[i];
	}
}

void MathVector::Initialize_const (double e)
{
	for (int i = 0; i < size; i++)
	{
		elem[i] = e;
	}
}

MathVector::~MathVector ()
{
	if (size > 0 && elem != NULL)
		delete[] elem;
	//elem = NULL;
}

void MathVector::Print ()
{
	for (int i = 0; i < size; i++)
	{
		printf ("%lf\t", elem[i]);
	}
	printf ("\n");
}

void MathVector::MultiplyByNumber (double b)
{
	for (int i = 0; i < size; i++)
	{
		elem[i] *= b;
	}
}

void MathVector::Substract (const MathVector & c)	
{
	for (int i = 0; i < size && i < c.size; i++)
	{
		elem[i] -= c.elem[i];
	}
}

double MathVector::Norm ()
{
	return sqrt (Scalar_Product (*this));
}

double MathVector::Norm (int par_num_threads)
{
	return sqrt (Scalar_Product (*this, par_num_threads));
}

double MathVector::Norm_MKL ()
{
	return sqrt (Scalar_Product_MKL (*this));
}

void MathVector::Copy (const MathVector & c)	
{
	// no resizing?
	if (size != c.size)
	{
		if (elem != NULL)
			delete[] elem;
		size = c.size;
		elem = new double[size];
	}
	for (int i = 0; i < size/* && i < c.size*/; i++)
	{
		elem[i] = c.elem[i];
	}
}

void MathVector::FPrint (FILE * f)
{
	for (int i = 0; i < size; i++)
	{
		fprintf (f, "%.13lf\n", elem[i]);
	}
	//fprintf (f, "\n\n");
}

void MathVector::FPrint (char * file_name)
{
	FILE * file = fopen (file_name, "w");
	for (int i = 0; i < size; i++)
	{
		fprintf (file, "%.16lf\n", elem[i]);
	}
	fclose (file);

}

void MathVector::Add (const MathVector & c)	
{
	for (int i = 0; i < size && i < c.size; i++)
	{
		elem[i] += c.elem[i];
	}
}

void MathVector::Linear_Combination (double a, const MathVector & c, double b)
{
	for (int i = 0; i < size && i < c.size; i++)
	{
		elem[i] = a * elem[i] + b * c.elem[i];
	}
}

void MathVector::Linear_Combination (double a, const MathVector & c, double b, int par_num_threads)
{
	int sh_size = size < c.size ? size : c.size;

	if (sh_size > PARALLEL_MIN_SIZE && par_num_threads > 1)
	{
		P_BLAS::Linear_Combination (sh_size, a, elem, b, c.elem, par_num_threads);
	}
	else
	{
		Linear_Combination (a, c, b);
	}
}

void MathVector::Linear_Combination (const MathVector & d, double a, const MathVector & c, double b, int par_num_threads)
{
	int sh_size = size < c.size ? size : c.size;

	if (sh_size > PARALLEL_MIN_SIZE && par_num_threads > 1)
	{
		P_BLAS::Linear_Combination (sh_size, a, d.elem, b, c.elem, elem, par_num_threads);
	}
	else
	{
		Linear_Combination (d, a, c, b);
	}
}

void MathVector::Linear_Combination (const MathVector & d, double a, const MathVector & c, double b)
{
	for (int i = 0; i < size && i < c.size; i++)
	{
		elem[i] = a * d.elem[i] + b * c.elem[i];
	}
}

void MathVector::Linear_Combination_MKL (const MathVector & c, double a, double b)
{
	int incr = 1;
	cblas_daxpby (size, b, c.elem, incr, a, elem, incr);
}

double MathVector::Scalar_Product (const MathVector & c)
{
	double r = 0;
	for (int i = 0; i < size && i < c.size; i++) // yep.
	{
		r += c.elem[i] * elem[i];
	}
	return r;
}

double MathVector::Scalar_Product (const MathVector & c, int par_num_threads)
{
	int sh_size = size < c.size ? size : c.size;
	double r = 0.0;
	if (sh_size > PARALLEL_MIN_SIZE && par_num_threads > 1)
	{
		return P_BLAS::Scalar_Product (sh_size, elem, c.elem, par_num_threads);
	}
	else
	{
		return Scalar_Product (c);
	}
}

double MathVector::Scalar_Product_MKL (const MathVector & c)
{
	int incr = 1;
	return ddot (&size, elem, &incr, c.elem, &incr);
}

void MathVector::Cut_Last ()
{
	double * e;
	int old_size = size;
	int new_size = size - 1; 
	e = new double[new_size];

	for (int j = 0; j < new_size; j++) // copy first size - 1 elements to temp vector
	{
		e[j] = elem[j];
	}
	delete[] elem; // delete current vector's elements

	size = new_size; // reset the size
	elem = new double[new_size]; // allocate memory for elements with new sizeentionality
	for (int j = 0; j < new_size; j++)
	{
		elem[j] = e[j]; // copy previously stored elements
	}
	delete[] e; // delete allocated memory
}

double MathVector::sum_sq_dif (const MathVector & c)
{
	double r = 0.0;
	for (int i = 0; i < size; i++)
	{
		r += pow (elem[i] - c.elem[i], 2.0);
	}
	return sqrt (r);
}

double MathVector::getElem (int k) const
{
	return elem[k];
}

void MathVector::setElem (int k, double e)
{
	elem[k] = e;
}

int MathVector::getSize () const
{
	return size;
}

void MathVector::addToElem (int k, double e)
{
	elem[k] += e;
}

void MathVector::multiplyElem (int k, double e)
{
	elem[k] *= e;
}

void MathVector::getFromFile (FILE * file)
{
	for (int i = 0; i < size; i++)
	{
		fscanf (file, "%lf", &elem[i]);
	}
}

void MathVector::setSize (int n)
{
	if (elem != NULL)
		delete[] elem;

	size = n;
	elem = new double[size];
	Zero ();
}

MathVector & MathVector::operator=(const MathVector & v)
{
	if (elem != NULL)
		delete[] elem;

	size = v.size;
	elem = new double[size];
	for (int i = 0; i < size; i++)
	{
		elem[i] = v.elem[i];
	}
	return *this;
}

