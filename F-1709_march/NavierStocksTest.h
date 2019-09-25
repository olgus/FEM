#pragma once
#include "Time_Dependent_Task.h"

class NavierStocksTest : public Time_Dependent_Task <Triangular_Mesh>
{
private:
protected:
	virtual void set_n_systems () override; // sets the amount of equations of the task
	virtual void build_system (int k_system) override; // build system

	 // reset first cond and f_func
	virtual double function_f (int k_system, double * coordinates, int area) override; // returns value of function f in coordinates
	virtual double function_FCondition (int k_system, double * coordinates, int area, int boundary) override; // returns value for first condition in coordinates
public:
	NavierStocksTest ();
	NavierStocksTest (const NavierStocksTest & nst);
	~NavierStocksTest ();
};