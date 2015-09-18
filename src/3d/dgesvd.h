#ifndef DGESVD_H
#define DGESVD_H
/*
  This file has my implementation of the LAPACK routine dgesdd for
  C++.  This program solves for the singular value decomposition of a
  rectangular matrix A.  The function call is of the form

    void dgesdd(double **A, int m, int n, double *S, double *U, double *VT)

    A: the m by n matrix that we are decomposing
    m: the number of rows in A
    n: the number of columns in A (generally, n<m)
    S: a min(m,n) element array to hold the singular values of A
    U: a [m, min(m,n)] element rectangular array to hold the right
       singular vectors of A.  These vectors will be the columns of U,
       so that U[i][j] is the ith element of vector j.
    VT: a [min(m,n), n] element rectangular array to hold the left
        singular vectors of A.  These vectors will be the rows of VT
	(it is a transpose of the vector matrix), so that VT[i][j] is
	the jth element of vector i.

  Note that S, U, and VT must be initialized before calling this
  routine, or there will be an error.  Here is a quick sample piece of
  code to perform this initialization; in many cases, it can be lifted
  right from here into your program.
  
    S = new double[minmn];
    U = new double*[m]; for (int i=0; i<m; i++) U[i] = new double[minmn];
    VT = new double*[minmn]; for (int i=0; i<minmn; i++) VT[i] = new double[n];

  Scot Shaw
  24 January 2000 */

#include <math.h>
#include <iostream>
using namespace std;

// double cond(double **A, int m, int n);
	
void interfacedgesvd(double **A, int m, int n, double *S, double **U, double **VT);

double *dgesvd_ctof(double **in, int rows, int cols);

void dgesvd_ftoc(double *in, double **out, int rows, int cols);

/*
extern "C" void dgesvd_(char *jobu, char *jobvt, int *m, int *n,
			double *a, int *lda, double *s, double *u,
			int *ldu, double *vt, int *ldvt, double *work,
			int *lwork, int *info);
*/

#endif

