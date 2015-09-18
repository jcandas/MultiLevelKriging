#ifndef CCV_FASTINVERSE_H
#define CCV_FASTINVERSE_H
#include "linearSolver.h"
// It should be ok to use the same pointer for blurmap and rbfCenterValues
bool getBlurMapWithFastInverse( Context* ctx, fftw_complex* blurmap, fftw_complex* rbfCenterValues);
#endif


