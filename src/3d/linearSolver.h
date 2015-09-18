#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H
#include "mpi.h"
#include "petscksp.h"
#include "fwtd3.h"
#include "fmm3d.hpp"
#include "diag_blk.h"
#include "IO.h"

// extern "C" void dgetrf_(int*,int*,double*,int*,int*,int*);
// extern "C" void dgetrs_(char*,int*,int*,double*,int*,int*,double*,int*,int*);


double getTime();

typedef struct {
	int level;
	int numCenters;
	int numFreq;
	fftw_complex* RBF3D;
	int gridSize;
	double alpha;
	int interpFuncExtent;
	int rbfWidth;
	double* xk;
	double* yk;
	double* zk;
	FWTd3* fwtd3;
	KSP* solver;
	Mat* mD;
	Kernel3d* knl;
	FMM3d* fmm;
	FMM3d* fmm2;
	DblNumMat* centers;
	char* options_file;
	double* diagD;
	double* mask;
	int* msup;
	bool useFastDiag;
	bool usePMult;
	Diag_Blk* db;
	bool printResidual;
	Vec* rhs; 
	double* actualResidual;
} Context;


double maternderivativecopy(double nu, double rho, double d);
double maternfourthderivativecopy(double nu, double rho, double d);
double initializeadaptivecopy(double nu, double rho, double d, double Delta, double tol);
double adaptivehermitecopy(double nu, double rho, double d, double Delta, double tol);

void jacobi(int n, int lda, double* a, double* b, double* x, double& tol, int& max_iter);
void gaussSeidel(int n, int lda, double* a, double* b, double* x, double& tol, int& max_iter);
void lu_init(int n, int lda, double* a, int* ipiv);
void lu(int n, int lda, double* a, int* ipiv, double* b, double* x);
PetscErrorCode gmres_init(int argc, char* argv[], int n, int numCenters, int gridSize, double alpha, int d, int interpFuncExtent, double tol, int max_iter,int van,int pvan, double*rbfCenterPositions, fftw_complex* RBF3D, KSP& solver, Mat& mA, Mat& mD, int minlevel, double* b);
PetscErrorCode gmres(int n, KSP solver, Mat mA, double* b, double* x, double& tol, int& max_iter,double*rbfCenterPositions,int numCenters, double *rbfTargetCenterPositions, int numTargetCenters);
PetscErrorCode gmres_finalize(KSP& solver, Mat& mA, Mat& mD);
//Multiply with the A matrix
PetscErrorCode AMult(Mat mA, Vec x, Vec b);
//Multiply with a diagonal block of A matrix
PetscErrorCode DMult(Mat mD, Vec x, Vec b);
//Get the diagonal entries of D;
PetscErrorCode DGetDiag(Mat D, Vec diag);
//Multiply with the inverse of the SSOR preconditioner
PetscErrorCode PMult(PC pc, Vec x, Vec b);
//Multiplly with the sparse diagonal block
PetscErrorCode sparseMult(Mat mD, Vec x, Vec b);

// Evaluate \Phi_i' A \Phi_j for all \Phi_i in the list
void FMMevalHB(HB* hb, int* si, int siLength, int sj, double* eval, int* temp, Neighbor3* nei, double* dist, double* rout, double* site);

// Monitor the actual residual for each iteration
PetscErrorCode KSPMonitorResidual(KSP ksp, PetscInt n, PetscReal rnorm, void *ctx);
#endif
