#pragma once
#include <math.h>
#include <stdio.h>
#include "MathVector.h"

#define ZERO_DETERMINANT 1e-10 // check if determinant is close to the zero for inversion

class Matrix // 2D matrix. Do i need other dimentionalities? Doubt it
{
private:
	// dim = size0, dim2 = size1
	int size0, size1; // sizes of the matrix
	double ** elem; // elements of the matrix
public:	
	Matrix (); // contructor
	Matrix (int n); // constructor of the square matrix. Allocates the memory
	Matrix (const Matrix & matrix); // costructor-copy
	Matrix (int n, int m); // contructor of rectangular matrix. Allocates the memory
	~Matrix (); // destructor
	void Zero (); // sets 0-matrix
	void solve (MathVector * x, const MathVector & f);
	void IdentityMatrix (); // sets current matrix as diagonal with 1-s. Even rectangular, yep
	void Print (); // prints the matrix
	void FPrint (FILE * file); // prints the matrix into specified file
	void FPrint (char * file_name); // prints the matrix into specified file
	void MultiplyByNumber (double c); // multiplies each elements of the matrix by c
	void Add (const Matrix & m); // adds another matrix to this

	// the whole row and column thing is confusing tho
	// i-row goes by j to size1 and j-column is by i to size0
	// turns out it's not THAT confusing
	void Insert_Row (int nR, const MathVector & row); // inserts row into nR-row of the matrix 
	void Get_Row (int nR, MathVector * r); // gets nR-row from a matrix, puts it into r
	void Cut_Last_Row (); // cuts the last row of the matrix
	bool MultiplyByMatrix (Matrix const & sour, Matrix * dest); // multiply this by sour into dest
	void Rotation_Matrix (double c, double s, int start_block); // 1-diag matrix with a "cssc" block in start_block place
	void Copy (Matrix const & source); // copy source matrix
	void Full_Copy (Matrix const & source); // copy source matrix
	void Insert_Column (int nC, MathVector * column); // insert column in nC column of this
		
	void calc_URV (const MathVector & v, MathVector * vect); // solve U * r = v, where r and v are vectors
	void Get_Column (int nC, MathVector * r); // get nC-column of this and put it into r
	void MultiplyMatrixByVector (const MathVector & v, MathVector * r); // multiply this by v into r

	int Size0 () const; // returns size0
	int Size1 () const; // returns size1
	double Elem (int i, int j); // return ij-element

	void setElem (int i, int j, double v); // sets ij-element to v
	void addToElem (int i, int j, double v); // sets ij-element to v
	bool inverse_matrix_size_3 (Matrix * result); // inverse matrix, if its size = 3
	double get_determinant_3 (); // return matrix's determinant, if its size=3

	void Transpose ();
};
