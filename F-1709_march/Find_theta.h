#pragma once
#include "Melt_SF.h"
#include <time.h>

#define FIND_THETA_MAX_ITER 50
#define FIND_THETA_PRECIS 2e-7
#define BETA_PRECIS 0.00005
#define ZERO_PRECIS 2e-10

struct receiver_data
{
	double x;
	double y;
};

class Task_Find_Theta
{
private:
public:
	void solve (char * file_front_name, Painter * painter);
	void solve_ (char * file_front_name, Painter * painter);
	void solve_ (char * file_front_name, Painter * painter, double alpha);
};