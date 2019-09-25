#pragma once
#include "Prototype_Node.h"
#include "Point.h"

// 2D node in cartesian coordinates
class Node_2D : public Node <Point_2D>
{
private:
public:
	Node_2D (); // constructor
	Node_2D (const Node_2D & node);  // constructor-copy
	~Node_2D (); // destructor

	void set_Point (Point_2D p);
	Point_2D * get_Point ();
};
