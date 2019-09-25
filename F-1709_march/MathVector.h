#pragma once
#include <math.h>
#include <stdio.h>

#include <omp.h>
#include <mkl.h>
#include "Parallel_realisations.h"

#define PI 3.14159265358979
#define PARALLEL_MIN_SIZE 1000

class MathVector // vector that stores 1D-doubles
{	
private:
	int size; // amount of elements
public:
	double * elem; // data
	MathVector (); // contructor 
	MathVector (int n); // contructor that sets the amount of elements and allocates the memory, filling it with 0-s
	MathVector (const MathVector & mvector); // contructor-copy 
	~MathVector (); // destructor
	void MultiplyByNumber (double b); // multiplies each member of the vector by a number
	void Substract (const MathVector & c); // substracts the vector c from current
	void Add (const MathVector & c); // add the vector c to current
	void Linear_Combination (double a, const MathVector & c, double b); // add the vector b * c to current with a modifier a
	void Linear_Combination (double a, const MathVector & c, double b, int par_num_threads); // add the vector b * c to current with a modifiera
	void Linear_Combination (const MathVector & d, double a, const MathVector & c, double b, int par_num_threads); // add the vector b * c to current with a modifiera
	void Linear_Combination (const MathVector & d, double a, const MathVector & c, double b); // add the vector b * c to current with a modifiera
	void Linear_Combination_MKL (const MathVector & c, double a, double b); // add the vector b * c to current with a modifiera

	void Initialize (double * e); // initialize vector with elements from an array
	void Initialize_const (double e); // initialize vector with elements equals e
	void Print (); // prints the vector
	void FPrint (FILE * f); // prints the vector into the file in question
	void FPrint (char * file_name); // prints the vector into the file in question
	double Norm (); // return norm of the vector
	double Norm (int par_num_threads); // return norm of the vector
	double Norm_MKL (); // return norm of the vector
	void Copy (const MathVector & c); // copies the vector
	void Zero (); // sets the vector to zeroes, DOESN'T CHANGE THE DIM
	double Scalar_Product (const MathVector & c); // return scalar product of this and c vectors
	double Scalar_Product (const MathVector & c, int par_num_threads); // return scalar product of this and c vectors, parallel
	double Scalar_Product_MKL (const MathVector & c); // return scalar product of this and c vectors, parallel
	void Cut_Last (); // resizes the vector, deleting the last element
	double sum_sq_dif (const MathVector & c); // returns sqrt ((this - c)^2)

	double getElem (int k) const; // get k-element from vector
	void setElem (int k, double e); // set k-element's value to e
	void addToElem (int k, double e); // add e to k-element
	void multiplyElem (int k, double e); // multiplies k-element by e
	int getSize () const; // returns size of the vector
	void getFromFile (FILE * file); // sets elements of the vector from file
	void setSize (int n); // sets the amount of elements and allocates the memory, filling it with 0-s
	
	MathVector& operator= (const MathVector & v); // assignment operator
};

