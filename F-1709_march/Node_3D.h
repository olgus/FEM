#pragma once
#include "Prototype_Node.h"
#include "Point.h"

class Node_3D : public Node <Point_3D>
{
private:
protected:
public:
	Node_3D ();
	Node_3D (const Node_3D & node);
	~Node_3D ();

	void set_Point (Point_3D p);
	Point_3D * get_Point ();
};