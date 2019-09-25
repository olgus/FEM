#include "Node_2D_Polar.h"

Node_2D_Polar::Node_2D_Polar ()
{
	coordinates.set_dim (2);
}

Node_2D_Polar::Node_2D_Polar (const Node_2D_Polar & node)
{
	if (coordinates.Dim () == 0) // if node hasn't been set before
		coordinates.set_dim (2); // set dimentionality
	coordinates.set_point (node.coordinates); // set node's point 
}

Node_2D_Polar::~Node_2D_Polar ()
{
}

void Node_2D_Polar::set_Point (Point_2D_Polar p)
{
	coordinates.set_point (p); // set node's point 
}

Point_2D_Polar * Node_2D_Polar::get_Point ()
{
	return new Point_2D_Polar (coordinates);
}
