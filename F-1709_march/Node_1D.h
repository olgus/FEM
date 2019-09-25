#pragma once
#include "Prototype_Node.h"

class Node_1D : public Node <Point_1D>
{
private:
public:
	Node_1D ();
	Node_1D (const Node_1D & node);
	~Node_1D ();

	void set_Point (Point_1D p);
	Point_1D * get_Point ();
};