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

// #include <math.h>
#include <iostream>
#include "mkl_lapacke.h"
#include "dgesvd.h"

using namespace std;

void interfacedgesvd(double **A, int m, int n, double *S, double **U, double **VT)
{
  char jobu, jobvt;
  lapack_int lda, ldu, ldvt, lwork, info;
  double *a, *u, *vt, *work;

  int minmn, maxmn;

  jobu = 'N'; /* Specifies options for computing U.
		 A: all M columns of U are returned in array U;
		 S: the first min(m,n) columns of U (the left
		    singular vectors) are returned in the array U;
		 O: the first min(m,n) columns of U (the left
		    singular vectors) are overwritten on the array A;
		 N: no columns of U (no left singular vectors) are
		    computed. */

  jobvt = 'A'; /* Specifies options for computing VT.
		  A: all N rows of V**T are returned in the array
		     VT;
		  S: the first min(m,n) rows of V**T (the right
		     singular vectors) are returned in the array VT;
		  O: the first min(m,n) rows of V**T (the right
		     singular vectors) are overwritten on the array A;
		  N: no rows of V**T (no right singular vectors) are
		     computed. */

  lda = m; // The leading dimension of the matrix a.
  a = dgesvd_ctof(A, lda, n); /* Convert the matrix A from double pointer
			  C form to single pointer Fortran form. */

  ldu = m;

//  for(int ii = 0; ii < m; ii++){
//	  for(int jj = 0; jj < n; jj++)
//		  cout <<"  " << A[ii][jj];
//	  cout <<'\n';
//  }

  /* Since A is not a square matrix, we have to make some decisions
     based on which dimension is shorter. */

  if (m>=n) { minmn = n; maxmn = m; } else { minmn = m; maxmn = n; }

  ldu = m; // Left singular vector matrix
//  u = new double[ldu*minmn];

  ldvt = n; // Right singular vector matrix
  vt = new double[ldvt*n];

  lwork = 5*maxmn; // Set up the work array, larger than needed.
  work = new double[lwork];

  double superb[min(m,n)-1];

//  cout <<" calling dgesvd_\n";
  // lapack Fortran 95 interface
  // gesvd(&jobu, &jobvt, &m, &n, a, &lda, S, u, &ldu, vt, &ldvt, work, &lwork, &info);

  //lapack_int LAPACKE_dgesvd( int matrix_order, char jobu, char jobvt, lapack_int m, lapack_int n, double* a, lapack_int lda, double* s, double* u, lapack_int ldu, double* vt, lapack_int ldvt, double* superb );
LAPACKE_dgesvd(LAPACK_COL_MAJOR, jobu, jobvt, m, n, a, lda, S, u,
	       ldu, vt, ldvt, superb);
// info = LAPACKE_dgesvd( LAPACK_COL_MAJOR, 'A', 'A', m, n, a, lda,
//                         s, u, ldu, vt, ldvt, superb );



// dgesvd(LAPACK_COL_MAJOR, &jobu, &jobvt, &m, &n, a, &lda, S, u,
//          &ldu, vt, &ldvt, work, &lwork, &info);
//  cout <<" returned from dgesvd_\n";




//  dgesvd_ftoc(u, U, ldu, minmn);
  dgesvd_ftoc(vt, VT, ldvt, n);
 
//  for(int ii = 0; ii < n; ii++){
//	  for(int jj = 0; jj < n; jj++)
//		  cout <<"  " << VT[ii][jj];
//	  cout <<'\n';
//  }
  
  delete [] a;
//  delete [] u;
  delete [] vt;
  delete[] work;
}


double* dgesvd_ctof(double **in, int rows, int cols)
{
  double *out;
  int i, j;

  out = new double[rows*cols];
  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*rows] = in[i][j];
  return(out);
}


void dgesvd_ftoc(double *in, double **out, int rows, int cols)
{
  int i, j;

  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*rows];
}
