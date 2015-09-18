# MultiLevelKriging
Mutli Level Kriging Project


This code implements the method described in:

"Multi-Level Restricted Maximum Likelihood Covariance Estimation and
Kriging for Large Non-Gridded Spatial Datasets".  Julio E. Castrillon-Candas
,Marc G. Genton, and Rio Yokota

and

"A Discrete Adapted Hierarchical Basis Solver For Radial
Basis Function Interpolation", Julio E. Castrillon-Candas, 
Jun Li, Victor Eijkhout.

Note that we use some files from Volume Rover by Anthony Thane from
the University of Texas at Austin.


INSTALLATION:
-------------
Ubuntu 14.04, bash.

1.  Before compiling the code, you will need to install the following
software using sudo apt-get install:
sudo apt-get install gfortran
sudo apt-get install gcc
sudo apt-get install g++
sudo apt-get install libmpich2-dev
sudo apt-get install liblapack-dev
sudo apt-get install libcr-dev
sudo apt-get update && sudo apt-get install build-essential
sudo apt-get install libtbb-dev
sudo apt-get install libfftw3-dev
sudo apt-get install libboost-all-dev
sudo apt-get install libgsl0-dev

Install and Configure Intel MKL.

For example, for 64 bit machine use the following:

source /opt/intel/mkl/bin/mklvars.sh intel64 $*

2. Download and install PETSc 3.3

Try the following configuration:

./configure --configModules=PETSc.Configure
--optionsModule=PETSc.compilerOptions --with-shared-libraries=1
--with-clanguage=cxx --with-debugging=0 --COPTFLAGS="-O3
-mtune=native" --CXXOPTFLAGS="-O3 -mtune=native" --FOPTFLAGS="-O3
-mtune=native"

3. Compile the dependent packages:
CXSparse
SuiteSparse_config-4.2.1

4. 
a) Go to the Makefile in \krigingsparse
b) Change PETSC_DIR to point to PETSc
c) Set VPATH := ./source2d for 2D or VPATH := ./source3d for 3D.
d) Execute "make" under lonestar\krigingsparse

5. Once the code is compiled, you can test it with the following command line in
/krigingsparse:

export PETSC_OPTIONS='-ksp_type cg -ksp_view'

This sets petsc to use the preconditioned conjugate gradient for solving the  
Kriging problem.


3D Test:
--------

For n = 32,000 centers (Look at input32_coords.50.txt for the input
format).

Solving the Kriging Problem in 3D:
----------------------------------

./Kriging 2 2 2 2 ./data/centers/input32_coords.50.txt
./data/rhs/rand32000.dat soln.dat ./data/options/options_mat 4 4
9.88e-4 1 ./data/target/Target.txt 0 0

See below for the explanation for each term.

Computing the Sparse matrix with the distance criterion tau=1 in 2D, t
= 5 and i = 4.  Create a directory named: ~/Projects/Data/Kriging and
then run:

./Kriging 2 2 2 2 ./data/centers/input32_coords_2D.50.txt
./data/rhs/rand32000.dat soln.dat ./data/options/options_mat 5 4 1e-4
33 ./data/target/Target.txt 1 4

A Kw.dat file will be saved in the current directory where ./Kriging
is located.

6. MKL configuration at run time.

It is recommended to disconnect Hyperthreading from your machine. Also
the number of threads used by MKL can be controlled with the following
environment variables:

export MKL_NUM_THREADS=N
export OMP_NUM_THREADS=N

where N is the number of threads. Note, that using the making number
of threads is not necessarly the most efficient. Play around with N
until you obtain the best performance.

Notes:
------

Usage : ./Kriging alpha interpFuncExtent numFreq gridSize
inputFile1(coords) inputFile2(coefs) outputFileName rbfFile (ftilde+1)
(f+1) GMRES_tolerance [bitmap] Targetfile distance-criterion-cte
minlevel

   bitmap 1 = use fast diagonal computation (diagonal preconditioner)
   bitmap 2 = use PMult (SSOR preconditioner)
   bitmap 4 = write diagonal blocks
   bitmap 8 = read in the recomputed diagonal blocks
   bitmap 16 = compute actual residual for each iteration
   bitmap 32 = compute sparse Kw matrix and write it to file 

a) alpha interpFuncExtent numFreq gridSize: dummy variables left over
   for some legacy code. Just put 2 2 2 2 and forget about it.

b) inputFile1(coords): Coordinates of observation points. (Look at
input32_coords.50.txt for 3D and input32_coords_2D.50.txt for 2D).

c) inputFile2(coefs) : Observation values with the same order as a)

d) outputFileName : Solution to Kriging (\hat{\gamma}).

e) rbfFile: option file for KIFMM Method, Accuracy can be changed here
   (see options_mat file). Also, covariance function is defined
   here. kt set to 906 for Matern kernel, with alpha and rho.  Ignore
   numsrc and numtrg in this file, leave them at 1000.

f) ptilde: Total Degree of vanishing moments of Hierarchical Basis

g) p : Total Degree polynomial order

h) Tol : Tol of CG solver

i) Mode: 1 - Diagonal, for solving Kriging problem, 
        17 - Diagonal, for solving Kriging problem, output residual error for each iteration.
	33 - Construct and output sparse Kw matrix in /Projects/Data/Kriging

j) Targetfile: New locations where the kriging estimate is computed

k) (d) distance-criterion-cte: if Mode 1 is set then set d = 0
			       if Mode 33 is set then set d = 0 or d = 1, or d = 2, ...

m) minlevel (i): If mode 33 is set then the output matrix is computed from minlevel to maxHBlevel 
			


Major Modules:
--------------


- linearsolver.cpp
  Uses PETSc for GMRES itertaions

  Major functions: 
  PetscErrorCode gmres_init(int argc, char* argv[],
  int n, double tol, int max_iter, KSP& solver, Mat& mA, Mat& mD);

  PetscErrorCode gmres(int n, KSP solver, Mat mA, double* b, double*
  x, double& tol, int& max_iter);

  PetscErrorCode gmres_finalize(KSP& solver, Mat& mA, Mat& mD);
  PetscErrorCode AMult(Mat mA, Vec x, Vec b); //Use fast summation for
  A*x

  PetscErrorCode DMult(Mat mD, Vec x, Vec b); //Use fast summation for
  multiplying only with diagonal blocks of W'AW

  PetscErrorCode PMult(void* ctx, Vec x, Vec b); //Use blocked
  triagular solve to apply the preconditioner

- fwtd3.cpp
  Fast wavelet transform in 3D

  Major functions: 
  double multiLift(int cube,int van); //Lifting with
  SVD
  
  void trafo(double* d,int direction, int level); //Wavelet transform
  on vector


  related files includes:
  dgesvd.cpp  //computing SVD
  hb.cpp  //Hierarchical Basis
  matrixio.cpp // I/O
  neighbor3.cpp  // domain decomposition
  vande3.cpp  //moment matrix

- FastInverse.cpp
  Fast summation for scattered RBFs evaluated at scattered locations

  Major functions:
  
  bool getBlurMapWithFastInverse( fftw_complex* blurmap, fftw_complex*
  rbfCenterValues, double* rbfCenterPositions, fftw_complex* rbf1D,
  int numCenters, int numFreq, int gridSize, double alpha, int
  interpFuncExtent, int rbfWidth, SmoothingFunction*
  smoothingFunction, double* xk, double* yk, double* zk);

  related files includes:
  IO.cpp  // I/O
  

