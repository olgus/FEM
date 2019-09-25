#include "Node_3D.h"

Node_3D::Node_3D ()
{
	coordinates.set_dim (3);
}

Node_3D::Node_3D (const Node_3D & node)
{
	if (coordinates.Dim () == 0) // if node hasn't been set before
		coordinates.set_dim (3); // set dimentionality
	coordinates.set_point (node.coordinates); // set node's point 
}

Node_3D::~Node_3D ()
{
}

void Node_3D::set_Point (Point_3D p)
{
	coordinates.set_point (p); // set node's point 
}

Point_3D * Node_3D::get_Point ()
{
	return new Point_3D (coordinates);
}
