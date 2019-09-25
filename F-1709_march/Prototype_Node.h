#pragma once
#include <stdio.h>
#include <math.h.>
#include "Point.h"

#define ZERO_node 1e-10

// node prototype
template <class Point_Type> 
class Node
{
private:
protected:
	Point_Type coordinates; // coordinates
public:
	Node (); // constructor
	Node (const Node & node); // copy-constructor
	virtual ~Node (); // destructor
	virtual void get_coordinates (double * Coordinates) const /*override*/; // get node's coordinates
	virtual void set_coordinates (double * Coordinates); // set node's coordinates
	virtual int get_dim () const; // get node's dimentionality

	virtual void set_node (double * Coordinates); // set the node with a number and coordinates
	virtual double distance (const Node <Point_Type> & node); // get the distance between nodes
	virtual double distance (double * point); // get the distance between a node and a point
	virtual double distance_origin (); // get the distance between node and origin (0, 0)
};

// comparison betweeen nodes
template<class Point_Type>
inline bool operator==(const Node<Point_Type>& node1, const Node<Point_Type>& node2)
{
	if (node1.get_dim () != node2.get_dim ())
		return false;

	int n = node1.get_dim ();
	double * c0 = new double[n];
	node1.get_coordinates (c0);
	double * c1 = new double[n];
	node2.get_coordinates (c1);

	for (int i = 0; i < node1.get_dim (); i++)
	{
		if (fabs(c0[i] - c1[i]) > ZERO_node)
		{
			delete[] c0;
			delete[] c1;
			return false;
		}
	}

	delete[] c0;
	delete[] c1;
	return true;
}

