#pragma once

// point prototype
class Point_Prototype
{
protected:
	double * coordinates; // coordinates
	int dim; // dimentionality
public:
	Point_Prototype (); // constructor

	~Point_Prototype (); // destructor

	virtual void set_point (double * Coordinates); // set coordinates for the point
	void set_dim (int Dim); // set dimentionality

	virtual void get_point (double * Coordinates) const; // get coordinates for the point
	int Dim () const; // return dimentionality

	Point_Prototype& operator= (const Point_Prototype & point);
};

// point in cartesian coordinates
class Point_2D : public Point_Prototype
{
private:
public:
	using Point_Prototype::set_point;

	Point_2D (); // constructor
	Point_2D (double x, double y); // constructor with coordinates
	Point_2D (const Point_2D & point); // constructor-copy

	~Point_2D (); // destructor

	void set_point (double x, double y); // set point through x, y
	void set_point (Point_2D point); // set point with another point

	double X (); // get x
	double Y (); // get y 
};

// point in cartesian coordinates
class Point_2D_Polar : public Point_Prototype
{
private:
public:
	using Point_Prototype::set_point;

	Point_2D_Polar (); // constructor
	Point_2D_Polar (double x, double y); // constructor with coordinates
	Point_2D_Polar (const Point_2D_Polar & point); // constructor-copy

	~Point_2D_Polar (); // destructor

	void set_point (double x, double y); // set point through x, y
	void set_point (Point_2D_Polar point); // set point with another point

	double R (); // get x
	double P (); // get y 
};

class Point_1D : public Point_Prototype
{
private:
public:
	using Point_Prototype::set_point;

	Point_1D (); // constructor
	Point_1D (double x); // constructor with coordinates
	Point_1D (const Point_1D & point); // constructor-copy
	~Point_1D (); // destructor

	void set_point (double x); // set point through x, y
	void set_point (Point_1D point); // set point with another point

	double X (); // get x
};

class Point_3D : public Point_Prototype
{
private:
public:
	using Point_Prototype::set_point;

	Point_3D (); // constructor
	Point_3D (double x, double y, double z); // constructor with coordinates
	Point_3D (const Point_3D & point); // constructor-copy
	~Point_3D (); // destructor

	void set_point (double x, double y, double z); // set point through x, y
	void set_point (Point_3D point); // set point with another point

	double X (); // get x
	double Y (); // get y
	double Z (); // get z
};