#include <stdio.h>
#include "Point.h"

void Point_Prototype::set_point (double * Coordinates)
{
	if (dim != 0 && coordinates == NULL) // allocate memory if possible and it has not been allocated yet
		coordinates = new double[dim]; 

	for (int i = 0; i < dim; i++) // set coordinates
	{
		coordinates[i] = Coordinates[i];
	}
}

void Point_Prototype::get_point (double * Coordinates) const
{
	for (int i = 0; i < dim; i++)
	{
		Coordinates[i] = coordinates[i];
	}
}

Point_Prototype::Point_Prototype ()
{
	dim = 0;
	coordinates = NULL;
}

Point_Prototype::~Point_Prototype ()
{
	if (coordinates != NULL)
		delete[] coordinates;
	coordinates = NULL;
}

int Point_Prototype::Dim () const
{
	return dim;
}

Point_Prototype& Point_Prototype::operator=(const Point_Prototype & point)
{
	if (dim != point.dim && coordinates != NULL)
	{
		delete[] coordinates;
	}

	dim = point.Dim ();
	coordinates = new double[dim];
	for (int i = 0; i < dim; i++)
	{
		coordinates[i] = point.coordinates[i];
	}

	return *this;
}

void Point_Prototype::set_dim (int Dim)
{
	dim = Dim; // set dimentionality
}

Point_2D::Point_2D ()
{
	dim = 2;
}

Point_2D::Point_2D (double x, double y)
{
	dim = 2;
	if (coordinates == NULL) // if memory has not been allocated yet
	{
		coordinates = new double[dim]; // do it
	}
	coordinates[0] = x;
	coordinates[1] = y;
}

Point_2D::Point_2D (const Point_2D & point)
{
	dim = 2;
	if (coordinates != NULL)
	{
		delete[] coordinates;
	}
	coordinates = new double[dim]; // do it

	if (point.coordinates != NULL)
	{
		for (int i = 0; i < dim; i++)
		{
			coordinates[i] = point.coordinates[i];
		}
	}
}

Point_2D::~Point_2D ()
{
}

double Point_2D::X ()
{
	return coordinates[0];
}

double Point_2D::Y ()
{
	return coordinates[1];
}

void Point_2D::set_point (double x, double y)
{
	if (coordinates == NULL) // if memory has not been allocated yet
		coordinates = new double[dim];  // do it
	coordinates[0] = x;
	coordinates[1] = y;
}

void Point_2D::set_point (Point_2D point)
{
	if (coordinates == NULL) // if memory has not been allocated yet
		coordinates = new double[dim]; // do it
	if (point.coordinates != NULL)
	{
		for (int i = 0; i < dim; i++)
			coordinates[i] = point.coordinates[i];
	}
}

Point_2D_Polar::Point_2D_Polar ()
{
	dim = 2;
}

Point_2D_Polar::Point_2D_Polar (double x, double y)
{
	dim = 2;
	if (coordinates == NULL) // if memory has not been allocated yet
	{
		coordinates = new double[dim]; // do it
	}
	coordinates[0] = x;
	coordinates[1] = y;
}

Point_2D_Polar::Point_2D_Polar (const Point_2D_Polar & point)
{
	dim = 2;
	if (coordinates != NULL)
	{
		delete[] coordinates;
	}
	coordinates = new double[dim]; // do it

	if (point.coordinates != NULL)
	{
		for (int i = 0; i < dim; i++)
		{
			coordinates[i] = point.coordinates[i];
		}
	}
}

Point_2D_Polar::~Point_2D_Polar ()
{
}

double Point_2D_Polar::R ()
{
	return coordinates[0];
}

double Point_2D_Polar::P ()
{
	return coordinates[1];
}

void Point_2D_Polar::set_point (double x, double y)
{
	if (coordinates == NULL) // if memory has not been allocated yet
		coordinates = new double[dim];  // do it
	coordinates[0] = x;
	coordinates[1] = y;
}

void Point_2D_Polar::set_point (Point_2D_Polar point)
{
	if (coordinates == NULL) // if memory has not been allocated yet
		coordinates = new double[dim]; // do it
	for (int i = 0; i < dim; i++)
		coordinates[i] = point.coordinates[i];
}

Point_1D::Point_1D ()
{
	dim = 1;
}

Point_1D::Point_1D (double x)
{
	dim = 1;
	if (coordinates == NULL) // if memory has not been allocated yet
	{
		coordinates = new double[dim]; // do it
	}
	coordinates[0] = x;
}

Point_1D::Point_1D (const Point_1D & point)
{
	dim = 1;
	if (coordinates != NULL)
	{
		delete[] coordinates;
	}
	coordinates = new double[dim]; // do it

	if (point.coordinates != NULL)
	{
		for (int i = 0; i < dim; i++)
		{
			coordinates[i] = point.coordinates[i];
		}
	}
}

Point_1D::~Point_1D ()
{
}

void Point_1D::set_point (double x)
{
	if (coordinates == NULL) // if memory has not been allocated yet
		coordinates = new double[dim];  // do it
	coordinates[0] = x;
}

void Point_1D::set_point (Point_1D point)
{
	if (coordinates != NULL)
	{
		delete[] coordinates;
	}
	coordinates = new double[dim]; // do it

	if (point.coordinates != NULL)
	{
		for (int i = 0; i < dim; i++)
		{
			coordinates[i] = point.coordinates[i];
		}
	}
}

double Point_1D::X ()
{
	return coordinates[0];
}

Point_3D::Point_3D ()
{
	dim = 3;
	coordinates = NULL;
}

Point_3D::Point_3D (double x, double y, double z)
{
	dim = 3;
	if (coordinates == NULL) // if memory has not been allocated yet
	{
		coordinates = new double[dim]; // do it
	}
	coordinates[0] = x;
	coordinates[1] = y;
	coordinates[2] = z;
}

Point_3D::Point_3D (const Point_3D & point)
{
	dim = 3;
	if (coordinates != NULL)
	{
		delete[] coordinates;
	}
	coordinates = new double[dim]; // do it

	if (point.coordinates != NULL)
	{
		for (int i = 0; i < dim; i++)
		{
			coordinates[i] = point.coordinates[i];
		}
	}
}

Point_3D::~Point_3D ()
{
}

void Point_3D::set_point (double x, double y, double z)
{
	if (coordinates == NULL) // if memory has not been allocated yet
		coordinates = new double[dim];  // do it
	coordinates[0] = x;
	coordinates[1] = y;
	coordinates[2] = z;
}

void Point_3D::set_point (Point_3D point)
{
	if (coordinates == NULL) // if memory has not been allocated yet
		coordinates = new double[dim]; // do it
	if (point.coordinates != NULL)
	{
		for (int i = 0; i < dim; i++)
			coordinates[i] = point.coordinates[i];
	}
}

double Point_3D::X ()
{
	return coordinates[0];
}

double Point_3D::Y ()
{
	return coordinates[1];
}

double Point_3D::Z ()
{
	return coordinates[2];
}
