#include "Node_1D.h"

Node_1D::Node_1D ()
{
	coordinates.set_dim (1);
}

Node_1D::Node_1D (const Node_1D & node)
{
	if (coordinates.Dim () == 0) // if node hasn't been set before
		coordinates.set_dim (1); // set dimentionality
	coordinates.set_point (node.coordinates); // set node's point 
}

Node_1D::~Node_1D ()
{
}

void Node_1D::set_Point (Point_1D p)
{
	coordinates.set_point (p); // set node's point 
}

Point_1D * Node_1D::get_Point ()
{
	return new Point_1D (coordinates);
}
