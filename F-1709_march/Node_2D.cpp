#include "Node_2D.h"

Node_2D::Node_2D ()
{
	coordinates.set_dim (2);
}

Node_2D::Node_2D (const Node_2D & node)
{
	if (coordinates.Dim () == 0) // if node hasn't been set before
		coordinates.set_dim (2); // set dimentionality
	coordinates.set_point (node.coordinates); // set node's point 
}

Node_2D::~Node_2D ()
{
}

void Node_2D::set_Point (Point_2D p)
{
	coordinates.set_point (p); // set node's point 
}

Point_2D * Node_2D::get_Point ()
{
	return new Point_2D (coordinates);
}
