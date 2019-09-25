#pragma once
#include <vector>

class Spline_Pointer
{
private:
protected:
public:
	// save data_points, extra - alpha, beta for smoothing splines
	virtual void prepare (std::vector <std::pair<double, double>> points_data, double * extra) {};
	// solve
	virtual void solve_task () {};
	// get solution in point
	virtual bool get_solution_in_point (double * c, double * value) { return false; };
	// get boundaries
	virtual void get_boundaries (double * c0, double * cN) {};

	Spline_Pointer () {};
	Spline_Pointer (const Spline_Pointer & sp) {};
	virtual ~Spline_Pointer () {};
};