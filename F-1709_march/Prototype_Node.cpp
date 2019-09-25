#include "Prototype_Node.h"

template class Node <Point_Prototype>;
template class Node <Point_2D>;
template class Node <Point_2D_Polar>;
template class Node <Point_1D>;
template class Node <Point_3D>;

template <class Point_Type>
Node<Point_Type>::Node ()
{
	coordinates.set_dim (0);
}

template<class Point_Type>
Node<Point_Type>::Node (const Node & node)
{
	coordinates = node.coordinates;
}

template <class Point_Type>
Node <Point_Type>::~Node ()
{
}

template <class Point_Type>
void Node <Point_Type>::get_coordinates (double * Coordinates) const
{
	int n = coordinates.Dim ();
	double * c = new double[n]; 
	coordinates.get_point (c);
	for (int i = 0; i < coordinates.Dim(); i++)
		Coordinates[i] = c[i];
	delete[] c;
}

template <class Point_Type>
void Node <Point_Type>::set_coordinates (double * Coordinates)
{
	double * c = new double[coordinates.Dim ()];
	for (int i = 0; i < coordinates.Dim (); i++)
		c[i] = Coordinates[i];
	coordinates.set_point (c);
	delete[] c;
}

template<class Point_Type>
int Node<Point_Type>::get_dim () const
{
	return coordinates.Dim ();
}

template <class Point_Type>
void Node <Point_Type>::set_node (double * Coordinates)
{
	double * c = new double[coordinates.Dim ()];
	for (int i = 0; i < coordinates.Dim (); i++)
		c[i] = Coordinates[i];
	coordinates.set_point (c);
	delete[] c;
}

template<class Point_Type>
double Node<Point_Type>::distance (const Node <Point_Type> & node)
{
	double r = 0.0;
	double * c0 = new double[get_dim ()];
	get_coordinates (c0);
	double * c1 = new double[node.get_dim ()];
	node.get_coordinates (c1);

	for (int i = 0; i < coordinates.Dim (); i++)
	{
		r += pow (c0[i] - c1[i], 2.0);
	}

	delete[] c0;
	delete[] c1;

	return sqrt (r);
}

template<class Point_Type>
double Node<Point_Type>::distance (double * point)
{
	double * c0 = new double[get_dim ()];
	get_coordinates (c0);
	double r = 0.0;
	for (int i = 0; i < coordinates.Dim (); i++)
	{
		r += pow (point[i] - c0[i], 2.0);
	}
	delete[] c0;

	return sqrt (r);
}

template<class Point_Type>
double Node<Point_Type>::distance_origin ()
{
	double * c0 = new double[get_dim ()];
	get_coordinates (c0);
	double r = 0.0;
	for (int i = 0; i < coordinates.Dim (); i++)
	{
		r += pow (c0[i], 2.0);
	}
	delete[] c0;

	return sqrt (r);
}
