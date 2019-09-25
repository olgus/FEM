#pragma once
#include "Prototype_Node.h"
#include "Point.h"

// 2D node in cartesian coordinates
class Node_2D_Polar : public Node <Point_2D_Polar>
{
private:
public:
	Node_2D_Polar (); // constructor
	Node_2D_Polar (const Node_2D_Polar & node);  // constructor-copy
	~Node_2D_Polar (); // destructor

	void set_Point (Point_2D_Polar p);
	Point_2D_Polar * get_Point ();
};
