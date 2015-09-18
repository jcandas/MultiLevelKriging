// Fast Inverse Code
// Test

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

#include <fftw3.h>
#include "matrixio.h"
#include "linearSolver.h"
#include <gsl/gsl_sf_bessel.h>

#ifndef M_PI
#define 	M_PI   3.14159265358979323846
#endif

extern "C" void FMM(int,double*,double*,int,double*,double*,double,double);

bool getBlurMapWithFastInverse( Context* ctx, fftw_complex* blurmap, fftw_complex* rbfCenterValues){
  int m = ctx->numCenters - 2;

	DblNumVec srcden(m);
	double norm = ctx->gridSize * ctx->gridSize * ctx->gridSize;


	for(int i = 0; i < m; i++)
		srcden(i) = rbfCenterValues[i][0];
	DblNumVec targval(m);

	ctx->fmm->eval(srcden,targval);
	for(int i = 0; i < m; i++)
		blurmap[i][0] = targval(i) * norm; // Checking

	return true;
}

