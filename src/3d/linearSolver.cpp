#include "linearSolver.h"
#include "petscksp.h"
#include "FastInverse.h"
#include "fwtd3.h"
#include "matrixio.h"
#include "IO.h"
#include "dgesvd.h"
#include "vande3.h"
#include "mkl.h"

#include <gsl/gsl_sf_bessel.h>
#include <boost/math/special_functions/bessel.hpp>

#define A(i,j) (a[lda*(i)+(j)])
#define MAX_INT 1000
#define NDORDER 1

using namespace std;

// Interpolant
static int interpolantflag = 123;
static DFTaskPtr task;  // MKL Data Fitting task descriptor
static double K0;
static double rho;

extern bool getNormalizedPositions( double* rbfCenterPositions, double* xk, double* yk, double* zk, 
            int rbfWidth, int gridSize, int numCenters, int interpFuncExtent, int numFreq );

// Matern derivative
double maternderivativecopy(double nu, double rho, double d)
{
     double coef = (1 / (pow(2,nu-1) * tgamma(nu))) * (sqrt(2 * nu) / rho)  ;
     double z = sqrt(2 * nu) * d / rho; 
     double dkcenter =  (boost::math::cyl_bessel_k(nu-1, z) + boost::math::cyl_bessel_k(nu+1, z) );     
     double dkmodcenter = boost::math::cyl_bessel_k(nu + 1, z) - ( nu * boost::math::cyl_bessel_k(nu, z) / z );
     return coef * (nu * pow(z,nu-1) * boost::math::cyl_bessel_k(nu, z) - pow(z,nu) * dkmodcenter );
}

// Matern fourth derivative (used to compute the adaptive grid error)
double maternfourthderivativecopy(double nu, double rho, double d)
{
  double coef = (1 / (pow(2,nu-1) * tgamma(nu))); 
  double t2 = sqrt(nu);
  double t3 = 1.0/rho;
  double t4 = d*d;
  double t5 = rho*rho;
  double t6 = t5*t5;
  double t7 = nu*nu;
  double t8 = nu*(1.0/2.0);
  double t9 = sqrt(2.0);
  double t10 = d*t2*t3*t9;
  double t11 = d*t2*t3;
  double t12 = pow(t11,nu);
  double t0 = pow(2.0,t8+1.0)*coef*1.0/(d*d*d*d)*nu*1.0/(rho*rho*rho*rho)*t12*boost::math::cyl_bessel_k(nu,t10)*(t6*-6.0+nu*t6*2.2E1+t4*t5*3.0-t6*t7*2.4E1+nu*(t4*t4)*2.0-nu*t4*t5*1.2E1+nu*t6*t7*8.0+t4*t5*t7*1.2E1)-pow(2.0,t8+3.0/2.0)*coef*1.0/(d*d*d)*1.0/(rho*rho*rho)*t2*t12*(nu*2.0-1.0)*boost::math::cyl_bessel_k(nu+1.0,t10)*(t5*3.0+nu*t4*2.0-nu*t5*5.0+t5*t7*2.0);
  return coef * t0;
}

// Initialize adaptive matern
double initializeadaptivecopy(double nu, double rho, double d, double Delta, double tol)
{
  // Error at initial increment point
  double fourth_der_matern = fabs(maternfourthderivativecopy(nu, rho, d + Delta));
  double error = pow(Delta,4) *  fourth_der_matern / 384;
  while (error > tol){
      Delta = Delta / 2;
      // Error at new point
      error = fabs(pow(Delta,4) * fmax(fourth_der_matern, fabs(maternfourthderivativecopy(nu, rho, d + Delta))) / 384);
    }
  return Delta; 
}
// Adaptive matern
double adaptivehermitecopy(double nu, double rho, double d, double Delta, double tol)
{
  // Error at initial increment point
  double fourth_der_matern_1 = fabs(maternfourthderivativecopy(nu, rho, d));
  double fourth_der_matern_2 = fabs(maternfourthderivativecopy(nu, rho, d + Delta));
  double error = pow(Delta,4) *  fmax(fourth_der_matern_1, fourth_der_matern_2) / 384;
  double fourth_der_matern_new = 0;
  while (error > tol){
      Delta = Delta / 2;
      // Error at new point
      fourth_der_matern_new = fabs(maternfourthderivativecopy(nu, rho, d + Delta));
      error = fabs(pow(Delta,4) * fmax(fourth_der_matern_new,fourth_der_matern_1) / 384);
    }
  return Delta; 
}




void jacobi(int n, int lda, double* a, double* b, double* x, double &tol, int &max_iter){
	if (n <= 0) return; //empty system
	int c = 0; 
	double e = 10;  //norm squared of change in x
	double ne = 0;  //norm of rhs
	double t;
	double* y = new double[n];
	for(int i = 0; i < n; i++)
		ne += b[i] * b[i];
	while (c < max_iter && e > tol*tol){
		e = 0;
		for(int i = 0; i < n; i++){
			y[i] = b[i];
			// dealing with the tridiagonal system
			for(int j = i-1; j <= i+1; j+=2)
				if(j >= 0 && j < n) y[i] -= A(i,j)*x[j];
			t = y[i] - A(i,i)*x[i];
			y[i] /= A(i,i);
			e += t*t;
		}
		e /= (n+1)*ne;
		for(int i = 0; i < n; i++)
			x[i] = y[i];
		c++;
	}
	if (c < max_iter) max_iter = c;
	tol = sqrt(e);
	delete[] y;
}
void gaussSeidel(int n, int lda, double* a, double* b, double* x, double &tol, int &max_iter){
	if (n <= 0) return; //empty system
	int c = 0; 
	double e = 10;  //norm squared of change in x
	double ne = 0;  //norm of rhs
	double t;
	double* y = new double[n];
	for(int i = 0; i < n; i++)
		ne += b[i] * b[i];
	while (c < max_iter && e > tol*tol){
		e = 0;
		for(int i = 0; i < n; i++){
			y[i] = b[i];
			// dealing with the tridiagonal system
			for(int j = i-1; j <= i+1; j+=2)
				if(j >= 0 && j < n) y[i] -= A(i,j)*x[j];
			t = y[i] - A(i,i)*x[i];
			y[i] /= A(i,i);
			x[i] = y[i];
			e += t*t;
		}
		e /= (n+1)*ne;
		c++;
	}
	if (c < max_iter) max_iter = c;
	tol = sqrt(e);
	delete[] y;
}

void lu_init(int n, int lda, double* a, int* ipiv){
	int info;
	if (n <= 0) return; //empty system
	dgetrf_(&n,&n,a,&lda,ipiv,&info);
}

void lu(int n, int lda, double* a, int* ipiv, double* b, double* x){
	int nrhs = 1;
	int ldb = n;
	int info;
	if (n <= 0) return; //empty system
	for(int i = 0; i < n; i++)
		x[i] = b[i];
	dgetrs_("N",&n,&nrhs,a,&lda,ipiv,x,&ldb,&info);
}

static char help[] = "Solves a linear system with KSP.\n\n";

int optionsCreate(const char* optfile, map<string,string>& options)
{
  options.clear();
  ifstream fin(optfile);
  if(fin.good()==false) {
    printf("Options File %s \n", optfile);
         cerr<<"wrong option file"<<endl;        exit(1);
  }
  string name;  fin>>name;
  while(fin.good()) {
         char cont[100];         fin.getline(cont, 99);
         options[name] = string(cont);
         fin>>name;
  }
  fin.close();
  return 0;
}

PetscErrorCode gmres_init(int argc, char* argv[], int n, int numCenters, int gridSize, double alpha, int d, 
			  int interpFuncExtent, double tol, int max_iter,int van, int pvan,double*rbfCenterPositions, 
			  fftw_complex* RBF3D, KSP& solver, Mat& mA, Mat& mD, int minlevel, double* b){
	PetscErrorCode ierr;	
	PC pc,pc2;
	int rank;
	int bitmap = 3;

	printf("Enerting gmres_init ...\n");
	Context* ctx = new Context[1];
	ctx->numCenters = numCenters;
	ctx->gridSize = gridSize;
	ctx->numFreq = gridSize;
        ctx->RBF3D = RBF3D;
        int rbfSize = gridSize ;
	ctx->rbfWidth = gridSize /2;
	ctx->alpha = alpha;
	ctx->interpFuncExtent = interpFuncExtent;
        ctx->xk = new double[numCenters];
        ctx->yk = new double[numCenters];
        ctx->zk = new double[numCenters];
	double** coord = new double*[ctx->numCenters];
	for (int i = 0; i < ctx->numCenters; i++)
		coord[i] = &rbfCenterPositions[i*3];
	ctx->fwtd3 = new FWTd3(coord,n,van,pvan,d,minlevel);
	printf("gmres_init: FWTd3 initialized.\n");///
	ctx->options_file = argv[8];
	ctx->useFastDiag = true;
	ctx->usePMult = true;
	if (argc > 12)	{
		bitmap = atoi(argv[12]);
		ctx->useFastDiag = (bitmap & 1); //first bit
		ctx->usePMult = (bitmap & 2); //second bit
		ctx->printResidual = (bitmap & 16); //fifth bit
	}
	ctx->diagD = new double[numCenters];
	ctx->mask = new double[numCenters];
	ctx->msup = new int[numCenters];
	ctx->db = new Diag_Blk[1];
	ctx->db->max_level = 20; //set maximum level
	ctx->db->sparse_level = 3; 
	ctx->db->sparse_blk = new cs*[ctx->db->max_level
				- ctx->db->sparse_level + 1];
	ctx->db->blk_size = new int[ctx->db->max_level];
	for(int i = 0; i <= ctx->fwtd3->hb->wavMaxLevel; i++)
		ctx->db->blk_size[i] = ctx->fwtd3->hb->wavlevels[i][1]
			- ctx->fwtd3->hb->wavlevels[i][0];
	ctx->db->dense_blk = new double*[ctx->db->max_level];
	ctx->db->dense_lu = new double*[ctx->db->max_level];
	ctx->db->ipiv = new int*[ctx->db->max_level];
	printf("gmres_init: ctx->db initializaed.\n");///

	map<string,string> opts;
	optionsCreate(ctx->options_file, opts);

	vector<double> tmp(2);         

  	map<string,string>::iterator mi;
	mi = opts.find("-alpha"); assert(mi!=opts.end());

	{ istringstream ss((*mi).second);  ss>>tmp[0]; }
	mi = opts.find("-kt"); assert(mi!=opts.end());
	int kt; // = 901; //r^alpha
	{ istringstream ss((*mi).second);  ss>>kt; }
	printf("Initializing ctx->knl (kt = %d, alpha = %f)...\n",kt,tmp[0]);

	mi = opts.find("-rho"); assert(mi!=opts.end());
	{ istringstream ss((*mi).second);  ss>>tmp[1]; }
	rho = tmp[1];

	/* Initialize Matern kernel */
	double nu = tmp[0];
	rho   = tmp[1];
	double tol_mat = 5e-9; // Tolerance error for matern kernel
	double gridmaxincrement = 0.1;
	int    maxsize = 1000;
	double coef = 1 / (pow(2,nu-1) * tgamma(nu));	
	double maxrange = 10;
	double a2     = sqrt(2 * nu);	
	double mid = 1e-8;
	double zr = sqrt(2 * nu) * mid / rho;
	K0 = coef * boost::math::cyl_bessel_k(nu, zr) * pow(zr, nu); // MATERN
	double Z0 = maternderivativecopy(nu, rho, mid);
	int counter = 0;

	double xint[maxsize], yint[maxsize], zint[maxsize];
	double dsource = mid;
	xint[0] = mid;
	yint[0] = K0;
	zint[0] = Z0;
	// Initialize grid size
	double Delta = initializeadaptivecopy(nu, rho, dsource, gridmaxincrement, tol_mat);
   
	// Loop to put the rest of the values and create the adaptive grid such that 
	// error < tol

	while (dsource<maxrange){
	  counter = counter + 1;
	  dsource = dsource + Delta;
     
	  // MATERN
	  xint[counter] = dsource;
	  yint[counter] = coef * boost::math::cyl_bessel_k(nu, dsource * sqrt(2 * nu) / rho) * pow(sqrt(2 * nu) * dsource / rho, nu);
	  zint[counter] = maternderivativecopy(nu, rho, dsource);
	     
	  // Next adaptive time step
	  Delta = adaptivehermitecopy(nu, rho, dsource, gridmaxincrement, tol_mat);
	}
	// interpolant	
	int N = counter;
	int NSCOEFF = (N-1)*DF_PP_CUBIC;
	int status;

	// Setup Hermite spline
	MKL_INT sorder;                     // spline order
	MKL_INT stype;                      // spline type
	MKL_INT nx;                         // number of break points
	MKL_INT xhint;                      // additional info about break points
	MKL_INT ny;                         // number of functions
	MKL_INT yhint;                      // additional info about function
	MKL_INT nscoeff;                    // number of spline coefficients
	MKL_INT scoeffhint;                 // additional info about spline
	// coefficients
	MKL_INT bc_type;                    // boundary conditions type
	MKL_INT nbc;                        // number of boundary conditions
	MKL_INT ic_type;                    // internal conditions type
	MKL_INT nic;                        // number of internal conditions	   
	int NIC = N - 2;
	double ic[NIC];                      // array of internal conditions
	double bc[2];                        // array of boundary conditions
	double *datahint;                    // additional info about structure
	                                       // of arrays x and y
	double scoeff[ 2 * NSCOEFF];        // array of spline coefficients
	   
	// Initializing parameters for Data Fitting task
	sorder = DF_PP_CUBIC;
	stype  = DF_PP_HERMITE;
	
	// Parameters describing interpolation interval
	nx          = N;
	xhint       = DF_NON_UNIFORM_PARTITION;
	   
	// Parameters describing function
	ny          = 1;
	yhint       = DF_NO_HINT;
	   
	// Parameters describing spline coefficients storage
	nscoeff     = NSCOEFF;
	scoeffhint  = DF_NO_HINT;
	bc_type     = DF_BC_1ST_LEFT_DER | DF_BC_1ST_RIGHT_DER;
	   
	// Parameters describing internal conditions
	ic_type     = DF_IC_1ST_DER;
	nic         = NIC;
	scoeffhint  = DF_NO_HINT;
	bc_type     = DF_BC_1ST_LEFT_DER | DF_BC_1ST_RIGHT_DER;
	ic_type     = DF_IC_1ST_DER;
	xhint       = DF_NON_UNIFORM_PARTITION;

	// Internal derivative
	for(int i = 0; i < NIC; i++){
	  ic[i] = zint[i + 1];
	}

	// Boundary conditions data
	bc[0] = zint[0];
	bc[1] = zint[nx - 1];

	// Create Data Fitting task
	status = dfdNewTask1D( &task, nx, xint, xhint, ny, yint, yhint );
	   
	printf("Error code = %d \n", status);

	// Edit task parameters for Hermite cubic spline construction
	status = dfdEditPPSpline1D( task, sorder, stype, bc_type, bc,
				    ic_type, ic, scoeff, scoeffhint );
	// CheckDfError(errcode);
	printf("Error code = %d \n", status);


	// Construct Hermite cubic spline using STD method
	status = dfdConstruct1D( task, DF_PP_SPLINE, DF_METHOD_STD );
	printf("Error code = %d \n", status);
	printf("K0 = %f coef = %f Set Interpolant %d \n",K0,coef,interpolantflag);
	interpolantflag = 999;
	printf("nu = %f rho = %f Set Interpolant %d \n",nu,rho,interpolantflag);
	printf("Size of adaptive grid %d for error < %e \n", counter, tol_mat);
	

	ctx->fwtd3->hb->options_file = ctx->options_file;
	ctx->fwtd3->hb->kernelType = kt;
	ctx->fwtd3->hb->kernelParameter = tmp[0];

	// Evaluate
	ctx->knl = new Kernel3d(kt,tmp);

	// Assuming that the numCenters has 2 extra for framing
	printf("Allocating memory for ctx->centers (%d by %d)...\n",3,ctx->numCenters-2);
	ctx->centers = new DblNumMat(3,ctx->numCenters-2);
	//fixed here a bug of indexing error
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < ctx->numCenters-2; j++)
			(*(ctx->centers))(i,j) =  rbfCenterPositions[i+3*j];
	printf("done.\n");
	ctx->fmm = new FMM3d("fmm3d_");
  	ctx->fmm->srcpos() = ctx->centers;  
        ctx->fmm->srcnor() = ctx->centers;
  	ctx->fmm->trgpos() = ctx->centers;
	ctx->fmm->ctr() = Point3(0,0,0); // CENTER OF THE TOPLEVEL BOX
	ctx->fmm->rootlvl() = 0;         // 2^(-rootlvl) is the RADIUS OF THE TOPLEVEL BOX
	ctx->fmm->knl() = *(ctx->knl);
  	ctx->fmm->setup(opts);

	/// testing starts here
	/*
	printf("TESTING: n = %d\n",n);
	double** vv = new double*[n]; 
	for(int i = 0; i < n; i++)
		vv[i] = new double[n];
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++)
			vv[i][j] = 0;
		vv[i][i] = 1;
		ctx->fwtd3->trafo(vv[i],FWT_BACKWARD,-1);
	}
	writeMatrix("vv.dat",n,n,vv);
	printf("TESTING: vv written to file\n");
	for(int i = 0; i < n; i++)
		delete[] vv[i];
	delete[] vv;
	*/
	/// testing ends here
	
	//Compute the coarsest level scaling functions
	int wn = ctx->fwtd3->hb->waveletSize;
	printf("Computing L: n = %d, wn = %d\n",n,wn);
	double** vv1 = new double*[n-wn]; 
	for(int i = 0; i < n-wn; i++)
	  vv1[i] = new double[n];
	for(int i = wn; i < n; i++){
	  for(int j = 0; j < n; j++)
	    vv1[i-wn][j] = 0;
	  vv1[i-wn][i] = 1;
	  ctx->fwtd3->trafo(vv1[i-wn],FWT_BACKWARD,-1);
	}
	char fname[256];
	sprintf(fname,"%sL.dat",argv[7]);
	writeMatrix(fname,n-wn,n,vv1);
	// TESTING: vv1 written to file
	for(int i = 0; i < n-wn; i++)
	  delete[] vv1[i];
	delete[] vv1;

        if( !getNormalizedPositions( rbfCenterPositions, ctx->xk, ctx->yk, ctx->zk, ctx->rbfWidth, 
				     gridSize, numCenters, ctx->interpFuncExtent, ctx->numFreq ) )
        {
                printf("Could not normalize\n");
                fftw_free( RBF3D );
                delete [] rbfCenterPositions;
                return -1;
        }

	printf("Starting Petsc...\n");///
	PetscFunctionBegin;
	PetscInitialize(&argc,&argv,(char *)0,help);

	printf("initializing mA ...\n");///
	//mA is the operator of WAW'
	ierr = MatCreateShell(PETSC_COMM_SELF,wn,wn,PETSC_DECIDE,PETSC_DECIDE,
		(void*)ctx,&mA); CHKERRQ(ierr);
	ierr = MatShellSetOperation(mA, MATOP_MULT, (void(*)(void)) AMult); CHKERRQ(ierr);
	ierr = MatSetFromOptions(mA); CHKERRQ(ierr);
	ierr = MatAssemblyBegin(mA,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(mA,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	double* ee;
	Vec Vee;

	/*
	/// initializing TEST of TKT' starts here
	printf("Computing test TKT': wn = %d\n",wn);
	double** kk = new double*[wn];
	for(int i = 0; i < wn; i++)
		kk[i] = new double[wn];
	ierr = VecCreate(MPI_COMM_SELF,&Vee); CHKERRQ(ierr);
	ierr = VecSetSizes(Vee,PETSC_DECIDE,wn); CHKERRQ(ierr);
	ierr = VecSetFromOptions(Vee); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(Vee); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(Vee); CHKERRQ(ierr);
	for(int i = 0; i < wn; i++){
		ierr = VecSet(Vee,0.0); CHKERRQ(ierr);
		ierr = VecGetArray(Vee, &ee); CHKERRQ(ierr);
		ee[i] = 1;
		ierr = VecRestoreArray(Vee, &ee); CHKERRQ(ierr);
		ierr = AMult(mA,Vee,Vee);
		ierr = VecGetArray(Vee, &ee); CHKERRQ(ierr);
		for(int j = 0; j < wn; j++)
			kk[j][i] = ee[j];
		ierr = VecRestoreArray(Vee, &ee); CHKERRQ(ierr);
	}
	VecDestroy(&Vee);
	writeMatrix("kk.dat",wn,wn,kk);
	printf("TESTING: vv written to file\n");
	for(int i = 0; i < wn; i++)
		delete[] kk[i];
	delete[] kk;
	*/

	/// initializing TEST ends here

	printf("initializing mD ...\n");///
	//mD is the operator of (W_i'AW_i)^{-1} with levels in ct x
	ierr = MatCreateShell(PETSC_COMM_SELF,wn,wn,PETSC_DECIDE,PETSC_DECIDE,
		(void*)ctx,&mD); CHKERRQ(ierr);
	if(ctx->useFastDiag){
		ierr = MatShellSetOperation(mD, MATOP_MULT, (void(*)(void)) sparseMult); CHKERRQ(ierr);
	}else{
		ierr = MatShellSetOperation(mD, MATOP_MULT, (void(*)(void)) DMult); CHKERRQ(ierr);
	}
	ierr = MatSetFromOptions(mD); CHKERRQ(ierr);
	
	ierr = MatAssemblyBegin(mD,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(mD,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	if(0){
	double tint0, tint1;
	/// initializing diagonal of TKT' starts here
	//  do this for the coarsest level
	printf("Computing diagonal: wn = %d\n",wn);
	ierr = VecCreate(MPI_COMM_SELF,&Vee); CHKERRQ(ierr); 
	ierr = VecSetSizes(Vee,PETSC_DECIDE,wn); CHKERRQ(ierr);
	ierr = VecSetFromOptions(Vee); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(Vee); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(Vee); CHKERRQ(ierr);
	
	tint0 = dsecnd();
	for(int level = 0; level >=0; level--){
	  printf("level: %d\n",level);
		ctx->level = level;
		for(int i = ctx->fwtd3->hb->wavlevels[level][0]; 
			i < ctx->fwtd3->hb->wavlevels[level][1]; i++){
		//printf("Computing wavelet node %d\n",i);///
			ierr = VecSet(Vee,0.0); CHKERRQ(ierr);
			ierr = VecGetArray(Vee, &ee); CHKERRQ(ierr);
			ee[i] = 1;
			ierr = VecRestoreArray(Vee, &ee); CHKERRQ(ierr);
			ierr = DMult(mD,Vee,Vee);
			ierr = VecGetArray(Vee, &ee); CHKERRQ(ierr);
			ctx->diagD[i] = fabs(ee[i]);
			ierr = VecRestoreArray(Vee, &ee); CHKERRQ(ierr);
		}
	}
	tint1 = dsecnd();
	printf("Used %lf seconds \n",tint1 - tint0);
	VecDestroy(&Vee);

	//do this for the finer levels
	double* mask = ctx->mask;
	int* msup = ctx->msup;
	int mlen;
	DblNumMat* centers;
	DblNumVec* srcden;
	DblNumVec* targval;
	char numchar[100];
	double diagval;
	double* wx = new double[n];
	for(int level = ctx->fwtd3->hb->wavMaxLevel; level > 0; level--){
	  printf("level: %d\n",level);
  	  for(int i = ctx->fwtd3->hb->wavlevels[level][0]; 
		  i < ctx->fwtd3->hb->wavlevels[level][1]; i++){
                // printf("Computing wavelet node %d\n",i);
		for(int j = 0; j < n; j++)
			wx[j] = 0.0;
		wx[i] = 1.0;
		ctx->fwtd3->trafo(wx, FWT_BACKWARD, level);	
		mlen = 0;
	        for(int j = 0; j < n; j++)
	                if(wx[j] != 0){
        	                mask[mlen] = wx[j];
                	        msup[mlen] = j;
                        	mlen++;
	                }
 
		centers = new DblNumMat(3,mlen);
		for(int ii = 0; ii < 3; ii++)
			for(int jj = 0; jj < mlen; jj++)
				(*(centers))(ii,jj) = 
					(*(ctx->centers))(ii,msup[jj]);

		ctx->fmm2 = new FMM3d("fmm3d_");
		ctx->fmm2->srcpos() = centers;
		ctx->fmm2->srcnor() = centers;
		ctx->fmm2->trgpos() = centers;
		ctx->fmm2->ctr() = Point3(0,0,0);
		ctx->fmm2->rootlvl() = 0;
		ctx->fmm2->knl() = *(ctx->knl);
		sprintf(numchar,"%d",mlen);
		opts["-numsrc"] = string(numchar);
		opts["-numtrg"] = string(numchar);

		/* Initializing FMM */
		iC(ctx->fmm2->setup(opts));

		srcden = new DblNumVec(mlen);
		for(int j = 0; j < mlen; j++)
			(*srcden)(j) = mask[j];
		targval = new DblNumVec(mlen);
		ctx->fmm2->eval(*srcden, *targval);

		diagval = 0;
		for(int j = 0; j < mlen; j++)
			diagval += mask[j]*(*targval)(j);
		ctx->diagD[i] = fabs(diagval);

		delete ctx->fmm2;
		delete centers;
		delete targval;
		delete srcden;
	  }
	}
	delete[] wx;

	/* initializing diagonal ends here */
    }
	printf("initializing ksp ...\n");///
	ierr = KSPCreate(PETSC_COMM_SELF,&solver); CHKERRQ(ierr);
	ierr = KSPSetOperators(solver,mA,mA,DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
	ierr = KSPGetPC(solver,&pc); CHKERRQ(ierr);
	ierr = KSPSetType(solver,KSPGMRES); CHKERRQ(ierr);
	if(ctx->usePMult){
		ierr = PCSetType(pc,PCSHELL); CHKERRQ(ierr);
		ierr = PCShellSetApply(pc,PMult); CHKERRQ(ierr);
		ierr = PCShellSetContext(pc,(void*)ctx); CHKERRQ(ierr);
	}else{
		if(ctx->useFastDiag){
			ierr = PCSetType(pc,PCJACOBI); CHKERRQ(ierr);
			ierr = MatShellSetOperation(mA, MATOP_GET_DIAGONAL, (void(*)(void)) DGetDiag); CHKERRQ(ierr);
		}else 
			ierr = PCSetType(pc,PCNONE); CHKERRQ(ierr);
	}

	ierr = KSPSetTolerances(solver,tol,PETSC_DEFAULT,PETSC_DEFAULT,max_iter); CHKERRQ(ierr);
	ierr = KSPGMRESSetRestart(solver,200); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
	KSP* solver2 = new KSP;
	ierr = KSPCreate(PETSC_COMM_SELF,solver2); CHKERRQ(ierr);
	ierr = KSPSetOperators(*solver2,mD,mD,DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
	ierr = KSPGetPC(*solver2,&pc2); CHKERRQ(ierr);
	ierr = KSPSetType(*solver2,KSPGMRES); CHKERRQ(ierr);
	printf("Using GMRES for inner iterations.\n");
		ierr = PCSetType(pc2,PCNONE); CHKERRQ(ierr);

	ierr = KSPSetFromOptions(*solver2); CHKERRQ(ierr);

	// Setting tolarance here smaller than tol
	ierr = KSPSetTolerances(*solver2,tol,PETSC_DEFAULT,PETSC_DEFAULT,max_iter); CHKERRQ(ierr);
	ierr = KSPGMRESSetRestart(*solver2,100); CHKERRQ(ierr);

	ctx->solver = solver2;
	ctx->mD = &mD;


	// Output file
	double* vb = new double[n];
	for(int i = 0; i < n; i++)
		vb[i] = b[i];
	((Context*)ctx)->fwtd3->trafo(vb, FWT_FORWARD, -1);//-1 = FWT_ALL_LEVELS
        writeMatrix("rhs_hb_w.dat", 1, n, &vb);///
	delete[] vb;

	if(ctx->useFastDiag){

	  // Build multilevel covariance (RBF) matrix
	  if(bitmap & 32)
	    ctx->fwtd3->multilevel_trafo_diag_blk(ctx->db);

	  else {

	    if(!(bitmap &8) || bitmap & 4)
	      ctx->fwtd3->trafo_diag_blk(ctx->db);

	    if(bitmap & 4)
	      write_diag_blk(argc, argv, ctx->fwtd3->hb, ctx->db);
	    else if(bitmap & 8)
	      read_diag_blk(argc, argv,ctx->db);

	    double t;
	    printf("Evaluating diagonal for level ");///
	    if(!ctx->usePMult) //using diagonal preconditioner only
	      for(int li = 0; li <= ctx->fwtd3->hb->wavMaxLevel; li++){
		printf("%d ",li);///
		// printf("(size %d) \n",ctx->db->blk_size[li]);///
		for(int i = 0; i < ctx->db->blk_size[li]; i++){
		  if((li < ctx->db->sparse_level && ctx->db->blk_size[li] < 1e4)||(van>pvan && li==ctx->fwtd3->hb->wavMaxLevel) ){
		    //printf("Accessing ctx->db->dense_blk[%d][%d]\n",li,i*ctx->db->blk_size[li]+i);///
		    t = ctx->db->dense_blk[li][i * ctx->db->blk_size[li] + i];
		  }else{
		    //printf("Accessing ctx->db->sparse_blk[%d]\n",li-ctx->db->sparse_level);///
		    t = ctx->fwtd3->hb->cs_getEntry(ctx->db->sparse_blk[li - ctx->db->sparse_level],i,i);
		  }
		  // printf("Accessing ctx->diagD[%d]\n",ctx->fwtd3->hb->wavlevels[li][0]+i);///
		  ctx->diagD[ctx->fwtd3->hb->wavlevels[li][0] + i] = t;
		  // printf("%d %e \n", i, t);
		}
	      }
	    printf("\n");///
	  }
	}
	if(ctx->printResidual)
	  ctx->actualResidual = new double[max_iter + 1];
	printf("initialization done.\n");///
	PetscFunctionReturn(0);
}

//Solve W*A*W'x = W*d for x, given d
PetscErrorCode gmres(int n, KSP solver, Mat mA, double* b, double* x, double& tol, int& max_iter,double *rbfCenterPositions, 
		     int numCenters, double *rbfTargetCenterPositions, int numTargetCenters){
	PetscErrorCode ierr;	
	PetscInt its;
	KSPConvergedReason reason;
	Vec rhs,sol,rhs1;
	int rank;
	PetscScalar *xx;
	PetscScalar neg_one = -1.0;
	PetscReal norm;
	double gridnorm;
	void* ctx;
	ierr = MatShellGetContext(mA, &ctx); CHKERRQ(ierr);
	int van = ((Context*)ctx)->fwtd3->vanishingMoment;
	int waveletSize = ((Context*)ctx)->fwtd3->hb->waveletSize;

	if (n <= 0)  //empty system
		PetscFunctionReturn(0);
	ierr = PetscPrintf(PETSC_COMM_SELF,"Initializing rhs...\n"); ///
	double* vb = new double[n];
	for(int i = 0; i < n; i++)
		vb[i] = b[i];

	((Context*)ctx)->fwtd3->trafo(vb, FWT_FORWARD, -1);//-1 = FWT_ALL_LEVELS
        writeMatrix("rhs_w.dat", 1, n, &vb);///

	ierr = VecCreate(MPI_COMM_SELF,&rhs); CHKERRQ(ierr); 
	ierr = VecSetSizes(rhs,PETSC_DECIDE,waveletSize); CHKERRQ(ierr);
	ierr = VecSetFromOptions(rhs); CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_SELF,&rank); CHKERRQ(ierr);
	
	if(rank == 0)
		for(int i = 0; i < waveletSize; i++){
			ierr = VecSetValues(rhs,1,&i,&vb[i],INSERT_VALUES); CHKERRQ(ierr);
		}
	ierr = PetscPrintf(PETSC_COMM_SELF,"Values set for rhs.\n"); ///
	ierr = VecAssemblyBegin(rhs); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(rhs); CHKERRQ(ierr);

	if(((Context*)ctx)->printResidual)
	{
		((Context*)ctx)->rhs = &rhs;
		ierr = VecNorm(rhs,NORM_2, &norm); CHKERRQ(ierr);///
		printf("gmres_init: rhs norm is %g\n",norm);///
	}
	
	double* vx = new double[n];
	for(int i = 0; i < n; i++)
		vx[i] = x[i];
	//Not necessary if vx is initialized as zeros
	((Context*)ctx)->fwtd3->trafo(vx, FWT_FORWARD, -1);//-1 = FWT_ALL_LEVELS

	ierr = PetscPrintf(PETSC_COMM_SELF,"Initializing sol...\n"); ///
	ierr = VecCreate(MPI_COMM_SELF,&sol); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_SELF,"Setting values for sol 1...\n"); ///
	ierr = VecSetSizes(sol,PETSC_DECIDE,waveletSize); CHKERRQ(ierr);
	ierr = VecSetFromOptions(sol); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_SELF,"Setting values for sol 2...\n"); ///
	if(rank == 0)
		for(int i = 0; i < waveletSize; i++){
			ierr = VecSetValues(sol,1,&i,&vx[i],INSERT_VALUES); CHKERRQ(ierr);
		}
	ierr = PetscPrintf(PETSC_COMM_SELF,"Values set for sol.\n"); ///
	ierr = VecAssemblyBegin(sol); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(sol); CHKERRQ(ierr);

	if(((Context*)ctx)->printResidual)
		ierr = KSPMonitorSet(solver, KSPMonitorResidual, ctx, PETSC_NULL);

	ierr = PetscPrintf(PETSC_COMM_SELF,"Calling KSP solve...\n");
	ierr = KSPSolve(solver,rhs,sol); CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_SELF,"Returned from KSP solve.\n");
	ierr = KSPGetConvergedReason(solver,&reason); CHKERRQ(ierr);

	if(reason == KSP_DIVERGED_INDEFINITE_PC){
		ierr == PetscPrintf(PETSC_COMM_SELF,"Divergence because of indefinite preconditioner\n");
	}else if (reason < 0) {
		PetscPrintf(PETSC_COMM_SELF,"Other kind of divergence (%d)\n",reason);
	}else{
		ierr = KSPGetIterationNumber(solver,&its); CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_SELF,"Convergence in %D iterations.\n",its);
		max_iter = its;
	}
	ierr = VecGetArray(sol, &xx); CHKERRQ(ierr);
	for(int i = 0; i < waveletSize; i++)
		x[i] = xx[i];
	for(int i = waveletSize; i < n; i++)
		x[i] = 0;
	writeMatrix("SolutionHBOnly.dat",1,n,&x);
	((Context*)ctx)->fwtd3->trafo(x, FWT_BACKWARD, -1);//-1 = FWT_ALL_LEVELS

	ierr = VecRestoreArray(sol, &xx); CHKERRQ(ierr);
	ierr = VecDuplicate(rhs,&rhs1); CHKERRQ(ierr);
	ierr = MatMult(mA,sol,rhs1); CHKERRQ(ierr);
	ierr = VecAXPY(rhs1,neg_one,rhs); CHKERRQ(ierr);
	ierr = VecNorm(rhs1,NORM_2, &norm); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_SELF,"  Norm of rhs1: %A...\n",norm);
	tol = norm;
	ierr = VecNorm(rhs,NORM_2, &norm); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_SELF,"  Norm of rhs: %A...\n",norm);
	tol /= norm;

	// c = L(d - Kx)
	fftw_complex* cx = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * n );
	for(int i = 0; i < n; i++){
		cx[i][0] = x[i];
		cx[i][1] = 0;
		vb[i] = b[i];
	}
	getBlurMapWithFastInverse( (Context*)ctx, cx, cx);
	int gridSize = ((Context*)ctx)->gridSize;
	gridnorm = gridSize * gridSize * gridSize;
	for(int i = 0; i < n; i++){
		vb[i] -= cx[i][0]/gridnorm;
		vx[i] = vb[i];
	}


	printf("Residual number %d /n", n);
	// Save residual to file: d - Kx, apply the vandemond matrix inverse
	// in matlab
	writeMatrix("Residual.dat", 1, n, &vb);  
        
        ((Context*)ctx)->fwtd3->trafo(vb, FWT_FORWARD, -1);//-1 = FWT_ALL_LEVELS
	for(int i = 0; i < n - waveletSize; i++){
		x[n+i] = vb[waveletSize+i];
	}
        // Check the overall residual
	for(int i = 0; i < waveletSize; i++)
		vb[i] = 0;
	((Context*)ctx)->fwtd3->trafo(vb, FWT_BACKWARD, -1);//-1 = FWT_ALL_LEVELSi
	tol = 0;
	for(int i = 0; i < n; i++){
		vx[i] -= vb[i];
		tol += vx[i]*vx[i];
	}
	tol = sqrt(tol)/norm;
	        
	// Output Fast FMM at Target Centers
	map<string,string> opts;	
	optionsCreate( ((Context*)ctx)->options_file, opts);
	DblNumMat* sourcecenters = new DblNumMat(3,n);

	// fixed here a bug of indexing error
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < n; j++)
		  (*(sourcecenters))(i,j) =  rbfCenterPositions[i+3*j];
	
	
	FMM3d* finalfmm = new FMM3d("fmm3d_");	
	
	finalfmm->srcpos() = sourcecenters;  
        finalfmm->srcnor() = sourcecenters;
	finalfmm->ctr() = Point3(0,0,0); // CENTER OF THE TOPLEVEL BOX
	finalfmm->rootlvl() = 0;         // 2^(-rootlvl) is the RADIUS OF THE TOPLEVEL BOX
	finalfmm->knl() = *((Context*)ctx)->knl;
        
	// Compute Kernel at target points
	double* outputx = new double[numTargetCenters];
	DblNumMat* Targetcenters;
        Targetcenters = new DblNumMat(3,numTargetCenters);
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < numTargetCenters; j++)
		  (*Targetcenters)(i,j) =  rbfTargetCenterPositions[i+3*j];
		  // (*(Targetcenters))(i,j) = 0;
        ((Context*)ctx)->fmm->trgpos() = Targetcenters;

	finalfmm->trgpos() = Targetcenters; // Update 
	finalfmm->setup(opts);

	DblNumVec srcden(n);
	DblNumVec targval(numTargetCenters);
	
	fftw_complex* targetx = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * numTargetCenters );

	for(int i = 0; i < n; i++)
	   srcden(i) = x[i];
	finalfmm->eval(srcden,targval);
	for(int i = 0; i < numTargetCenters; i++)
	  outputx[i] = targval(i);

	writeMatrix("Target.dat", 1, numTargetCenters, &outputx);
	writeMatrix("KrigingEstimate.dat",1,numTargetCenters,&outputx);

	// Delete release memory
	delete[] outputx;	
	delete finalfmm;
	delete sourcecenters;
	delete[] vb;
	delete[] vx;
	
	fftw_free(cx);
	VecDestroy(&rhs);
	VecDestroy(&rhs1);
	VecDestroy(&sol);
	PetscFunctionReturn(0);
}

PetscErrorCode gmres_finalize(KSP& solver, Mat& mA, Mat& mD){
	PetscErrorCode ierr;
	void* t;
	Context* ctx;
	ierr = MatShellGetContext(mA, &t); CHKERRQ(ierr);
	ctx = (Context*)t;
	ierr = MatDestroy(&mA); CHKERRQ(ierr);
	ierr = MatDestroy(&mD); CHKERRQ(ierr);
	ierr = KSPDestroy(&solver); CHKERRQ(ierr);
	ierr = KSPDestroy(ctx->solver); CHKERRQ(ierr);
	
	if(ctx->printResidual)
		delete[] ctx->actualResidual;
	delete[] ctx->xk;
	delete[] ctx->yk;
	delete[] ctx->zk;
	delete ctx->fwtd3;
	delete[] ctx->diagD;
	delete ctx->knl;
	delete ctx->fmm;
	delete ctx->centers;
	delete[] ctx->mask;
	delete[] ctx->msup;
	delete ctx;
	
	ierr = PetscFinalize(); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

//Multiplied by TKT' using fast summation
//Computing only the wavelet part (scaling function part ignored)
PetscErrorCode AMult(Mat mA, Vec x, Vec b){
	//printf("Entering AMult ...\n");///
	PetscErrorCode ierr;
	void* ctx;
	ierr = MatShellGetContext(mA, &ctx); CHKERRQ(ierr);

	int n = ((Context*)ctx)->fwtd3->nei->N;
	int waveletSize = ((Context*)ctx)->fwtd3->hb->waveletSize;
	double* xk = ((Context*)ctx)->xk;
	double* yk = ((Context*)ctx)->yk;
	double* zk = ((Context*)ctx)->zk;
	int numCenters = ((Context*)ctx)->numCenters;
	int numFreq = ((Context*)ctx)->numFreq;
	int gridSize = ((Context*)ctx)->gridSize;
	double alpha = ((Context*)ctx)->alpha;
	int interpFuncExtent = ((Context*)ctx)->interpFuncExtent;
	int rbfWidth = ((Context*)ctx)->rbfWidth;
	fftw_complex* RBF3D = ((Context*)ctx)->RBF3D;
	Vec* rhs = ((Context*)ctx)->rhs;

	double* vx = new double[n];
	double* vb = new double[n];
	double* xx;
	double* bb;
	fftw_complex* cx = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * numCenters );
	double* rbfCenterPositions = ((Context*)ctx)->fwtd3->nei->coord[0];
	ierr = VecGetArray(x, &xx); CHKERRQ(ierr);
	ierr = VecGetArray(b, &bb); CHKERRQ(ierr);
	for(int i = 0; i < waveletSize; i++)
		vx[i] = xx[i];
	for(int i = waveletSize; i < n; i++)
		vx[i] = 0;
	((Context*)ctx)->fwtd3->trafo(vx, FWT_BACKWARD, -1);
	for(int i = 0; i < numCenters; i++){
		cx[i][0] = vx[i];
		cx[i][1] = 0;
	}
	if( !getBlurMapWithFastInverse( (Context*)ctx, cx, cx))
                {
                        printf("Could not compute the blurmap\n");
                        fftw_free( RBF3D );
                        return -1;
                }
	double norm = gridSize * gridSize * gridSize;
	for(int i = 0; i < n; i++)
		vb[i] = cx[i][0]/norm;

	((Context*)ctx)->fwtd3->trafo(vb, FWT_FORWARD, -1);//-1 = FWT_ALL_LEVELS

	//printf("transform done.\n");///
	for(int i = 0; i < waveletSize; i++)
		bb[i] = vb[i];
	//printf("Freeing up cx...\n");///
	fftw_free(cx);
	//printf("Restoring xx...\n");///
	ierr = VecRestoreArray(x, &xx); CHKERRQ(ierr);
	//printf("Restoring bb...\n");///
	ierr = VecRestoreArray(b, &bb); CHKERRQ(ierr);

	//printf("Time : %14.3f\n",getTime());
	//printf("Freeing up vx...\n");///
	delete[] vx;
	//printf("Freeing up vb...\n");///
	delete[] vb;
	//printf("AMult done.\n");///
	PetscFunctionReturn(0);
}

//Multiply with the diagonal block at ctx->level
//data at levels other than ctx->level are zeroed out
PetscErrorCode DMult(Mat mD, Vec x, Vec b){
//	printf("Entering DMult ...\n");///
	PetscErrorCode ierr;
	void* t;
	Context* ctx;
	double* vx;
	double* vb;
	int n;
	int nx;

	MatShellGetContext(mD,&t);
	ctx = (Context*)t;
	ierr = VecGetArray(x, &vx); CHKERRQ(ierr);
	ierr = VecGetArray(b, &vb); CHKERRQ(ierr);
	ierr = VecGetSize(x, &nx); CHKERRQ(ierr);
	n = ctx->fwtd3->nei->N; //size of the computed centers
	double* rbfCenterPositions = ctx->fwtd3->nei->coord[0];	
	double gridSize = ctx->gridSize;
	fftw_complex* cx = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
	double norm = gridSize * gridSize * gridSize;

	double* wx = new double[n];
	for(int i = 0; i < n; i++)
		wx[i] = 0.0;
	for(int i = ctx->fwtd3->hb->wavlevels[ctx->level][0]; 
		i < ctx->fwtd3->hb->wavlevels[ctx->level][1]; i++)
		wx[i] = vx[i];
	ctx->fwtd3->trafo(wx, FWT_BACKWARD, ctx->level);	
	for(int i = 0; i < n; i++){
		cx[i][0] = wx[i];
		cx[i][1] = 0;
	}
	getBlurMapWithFastInverse( ctx, cx, cx);
	for(int i = 0; i < n; i++)
		wx[i] = cx[i][0]/norm;
	ctx->fwtd3->trafo(wx,FWT_FORWARD, ctx->level);
	// still to be optimized by reducing the range of operation
	// also possibly move cx into ctx
	for(int i = ctx->fwtd3->hb->wavlevels[ctx->level][0];
		i < ctx->fwtd3->hb->wavlevels[ctx->level][1]; i++)
		vb[i] = wx[i];

	fftw_free(cx);

	ierr = VecRestoreArray(x, &vx); CHKERRQ(ierr);
	ierr = VecRestoreArray(b, &vb); CHKERRQ(ierr);
	delete[] wx;
//	printf("DMult done.\n");///
	PetscFunctionReturn(0);
}

PetscErrorCode DGetDiag(Mat D, Vec diag){
	void* t;
	Context* ctx;
	double* dd;
	PetscErrorCode ierr;

	ierr = MatShellGetContext(D, &t); CHKERRQ(ierr);
	ctx = (Context*)t;
	ierr = VecGetArray(diag, &dd); CHKERRQ(ierr);
	//printf("Setting up diagnal preconditioning ...\n");///
	for(int i = 0;
		i < ctx->fwtd3->hb->waveletSize; i++)
		dd[i] = ctx->diagD[i];
	ierr = VecRestoreArray(diag, &dd); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

//Multiply by inv(M) where M = (D+L)*inv(D)*(D+L)'
//i.e. inv((D+L)')*D*inv(D+L)
PetscErrorCode PMult(PC pc, Vec x, Vec b){
	printf("Entering PMult ...\n");///
//Check the n's!!
	char fname[256];///
	PetscErrorCode ierr;
	KSPConvergedReason reason;
	PetscInt its;
	Vec sol,rhs;
	int max_iter;

	void* t;
	PCShellGetContext(pc, &t);
	Context* ctx = (Context*)t;
	double *vx,*vb;
	int n = ctx->fwtd3->nei->N;
	double gridSize = ctx->gridSize;
	double norm = gridSize * gridSize * gridSize;
	int rank;
	int nx; //actual size of the vector, which should be ctx->fwtd3->hb->waveletSize

	bool coarseLifted = false;
	if(ctx->fwtd3->vanishingMoment > ctx->fwtd3->pvanishingMoment)
		coarseLifted = true;

	ierr = VecGetArray(x, &vx); CHKERRQ(ierr);
	ierr = VecGetArray(b, &vb); CHKERRQ(ierr);
	ierr = VecGetSize(x, &nx); CHKERRQ(ierr);
	double* xx = new double[nx];//copy of x
	double* wb = new double[n];//a temporary storage for wavelet based variable
	double *vsol,*vrhs;//
	for(int i = 0; i < nx; i++){
		xx[i] = vx[i];
	}

	ierr = VecCreate(MPI_COMM_SELF,&rhs); CHKERRQ(ierr); 
	ierr = VecSetSizes(rhs,PETSC_DECIDE,nx); CHKERRQ(ierr);
	ierr = VecSetFromOptions(rhs); CHKERRQ(ierr);

	ierr = VecCreate(MPI_COMM_SELF,&sol); CHKERRQ(ierr); 
	ierr = VecSetSizes(sol,PETSC_DECIDE,nx); CHKERRQ(ierr);
	ierr = VecSetFromOptions(sol); CHKERRQ(ierr);
	ierr = VecSet(sol, 0); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(sol); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(sol); CHKERRQ(ierr);

	ierr = MPI_Comm_rank(PETSC_COMM_SELF,&rank); CHKERRQ(ierr);
	int maxLevel = ctx->fwtd3->hb->wavMaxLevel;
	fftw_complex* cx = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
	double* rbfCenterPositions = ctx->fwtd3->nei->coord[0];

	printf("Applying inv(D+L)\n");///
	// Apply inv(D+L)
	for(int level = maxLevel; level >= 0; level--){ //not including the coarsest level
	    if(ctx->useFastDiag && 
		((level < ctx->db->sparse_level && 
			ctx->db->blk_size[level] < 1e4)||
		(coarseLifted && level == maxLevel))){
		int blk_n = ctx->db->blk_size[level];
		double* a = ctx->db->dense_lu[level];
		int* ipiv = ctx->db->ipiv[level];
		double* b = &xx[ctx->fwtd3->hb->wavlevels[level][0]];
		double* x = &vb[ctx->fwtd3->hb->wavlevels[level][0]]; 
		lu(blk_n, blk_n, a, ipiv, b, x);
	    }else{
		ierr = VecSet(rhs, 0); CHKERRQ(ierr);
		if(rank == 0){
			for(int i = ctx->fwtd3->hb->wavlevels[level][0];
				i < ctx->fwtd3->hb->wavlevels[level][1]; i++)
				ierr = VecSetValues(rhs,1,&i,&xx[i],INSERT_VALUES); 
			CHKERRQ(ierr);
		}
		ctx->level = level;
		ierr = VecAssemblyBegin(rhs); CHKERRQ(ierr);
		ierr = VecAssemblyEnd(rhs); CHKERRQ(ierr);

		// testing the idea of multigrid
			ierr = KSPSetTolerances(*(ctx->solver),1e-6,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);
			ierr = KSPSolve(*(ctx->solver), rhs, sol);
			ierr = KSPGetConvergedReason(*(ctx->solver),&reason); CHKERRQ(ierr);
			if(reason == KSP_DIVERGED_INDEFINITE_PC){
				ierr == PetscPrintf(PETSC_COMM_SELF,"Divergence because of indefinite preconditioner\n");
			}else if (reason < 0) {
				PetscPrintf(PETSC_COMM_SELF,"Other kind of divergence (ErrorCode: %d)\n",reason);
			}else{
				ierr = KSPGetIterationNumber(*(ctx->solver),&its); CHKERRQ(ierr);
				PetscPrintf(PETSC_COMM_SELF,"Convergence in %D iterations.\n",its);
				max_iter = its;
			}
		// read out the solution at this level
		if(rank == 0){
			ierr = VecGetArray(sol, &vsol); CHKERRQ(ierr);
			for(int i = ctx->fwtd3->hb->wavlevels[level][0];
				i < ctx->fwtd3->hb->wavlevels[level][1]; i++)
				vb[i] = vsol[i];
			ierr = VecRestoreArray(sol, &vsol); CHKERRQ(ierr);
		}
	    }
		// Note that levels are arranged from max to min
		for(int ll =level-1; ll >= 0; ll--){
			for(int i = 0; i < nx; i++)
				wb[i] = vb[i];
			for(int i = nx; i < n; i++)
				wb[i] = 0.0;
			ctx->fwtd3->trafo(wb, FWT_BACKWARD, level);
			for(int i = 0; i < n; i++){
				cx[i][0] = wb[i];
				cx[i][1] = 0;
			}
			getBlurMapWithFastInverse( ctx, cx, cx);
			for(int i = 0; i < n; i++)
				wb[i] = cx[i][0]/norm;
			ctx->fwtd3->trafo(wb,FWT_FORWARD, ll);
			for(int i = ctx->fwtd3->hb->wavlevels[ll][0];
				i < ctx->fwtd3->hb->wavlevels[ll][1]; i++)
				xx[i] -= wb[i];
		}
	}	
	printf("Applying D\n");///
	//Apply D
	for(int level = maxLevel; level >= 0; level--){
	    if(ctx->useFastDiag && 
		((level < ctx->db->sparse_level && 
			ctx->db->blk_size[level] < 1e4)||
		(coarseLifted && level == maxLevel))){
		int bn = ctx->db->blk_size[level];
		int i0 = ctx->fwtd3->hb->wavlevels[level][0];
		double* A = ctx->db->dense_blk[level];
		for(int ii = 0; ii < bn; ii++){
			wb[i0 + ii] = 0;
			for(int jj = 0; jj < bn; jj++)
				wb[i0 + ii] += A[ii + jj*bn] * vb[i0 + jj];
		}
		for(int ii = 0; ii < bn; ii++)
			vb[i0 + ii] = wb[i0 + ii];
	    }else{
		ierr = VecSet(sol, 0); CHKERRQ(ierr);
		if(rank == 0){
			for(int i = ctx->fwtd3->hb->wavlevels[level][0];
				i < ctx->fwtd3->hb->wavlevels[level][1]; i++)
				ierr = VecSetValues(sol,1,&i,&vb[i],INSERT_VALUES); 
			CHKERRQ(ierr);
		}
		ctx->level = level;
		ierr = VecAssemblyBegin(sol); CHKERRQ(ierr);
		ierr = VecAssemblyEnd(sol); CHKERRQ(ierr);

		ierr = MatMult(*(ctx->mD), sol, rhs);
		if(rank == 0){
			ierr = VecGetArray(rhs, &vrhs); CHKERRQ(ierr);
			for(int i = ctx->fwtd3->hb->wavlevels[level][0];
				i < ctx->fwtd3->hb->wavlevels[level][1]; i++)
				vb[i] = vrhs[i];
			ierr = VecRestoreArray(rhs, &vrhs); CHKERRQ(ierr);
		}
	    }
	}
	printf("Applying inv(D+L)'\n");///
	// Apply inv((D+L)')
	for(int i = 0; i < nx; i++)
		xx[i] = vb[i];
	for(int level = 0; level <= maxLevel; level++){
	    if(ctx->useFastDiag && 
		((level < ctx->db->sparse_level && 
			ctx->db->blk_size[level] < 1e4)||
		(coarseLifted && level == maxLevel))){
		int blk_n = ctx->db->blk_size[level];
		double* a = ctx->db->dense_lu[level];
		int* ipiv = ctx->db->ipiv[level];
		double* b = &xx[ctx->fwtd3->hb->wavlevels[level][0]];
		double* x = &vb[ctx->fwtd3->hb->wavlevels[level][0]]; 
		lu(blk_n, blk_n, a, ipiv, b, x);
	     }else{
		ierr = VecSet(rhs, 0); CHKERRQ(ierr);
		if(rank == 0){
			for(int i = ctx->fwtd3->hb->wavlevels[level][0];
				i < ctx->fwtd3->hb->wavlevels[level][1]; i++)
				ierr = VecSetValues(rhs,1,&i,&xx[i],INSERT_VALUES); 
			CHKERRQ(ierr);
		}
		ctx->level = level;
		ierr = VecAssemblyBegin(rhs); CHKERRQ(ierr);
		ierr = VecAssemblyEnd(rhs); CHKERRQ(ierr);

		///testing the idea of multigrid
		ierr = KSPSetTolerances(*(ctx->solver),1e-6,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);
		ierr = KSPSolve(*(ctx->solver), rhs, sol);
		ierr = KSPGetConvergedReason(*(ctx->solver),&reason); CHKERRQ(ierr);
		if(reason == KSP_DIVERGED_INDEFINITE_PC){
			ierr == PetscPrintf(PETSC_COMM_SELF,"Divergence because of indefinite preconditioner\n");
		}else if (reason < 0) {
			PetscPrintf(PETSC_COMM_SELF,"Other kind of divergence\n");
		}else{
			ierr = KSPGetIterationNumber(*(ctx->solver),&its); CHKERRQ(ierr);
			PetscPrintf(PETSC_COMM_SELF,"Convergence in %D iterations.\n",its);
			max_iter = its;
		}
		// read out the solution at this level
		if(rank == 0){
			ierr = VecGetArray(sol, &vsol); CHKERRQ(ierr);
			for(int i = ctx->fwtd3->hb->wavlevels[level][0];
				i < ctx->fwtd3->hb->wavlevels[level][1]; i++)
				vb[i] = vsol[i];
			ierr = VecRestoreArray(sol, &vsol); CHKERRQ(ierr);
		}
	     }
		// Note that levels are arranged from max to min
		for(int ll =level+1; ll <= maxLevel; ll++){
			for(int i = 0; i < nx; i++)
				wb[i] = vb[i];
			for(int i = nx; i < n; i++)
				wb[i] = 0.0;
			ctx->fwtd3->trafo(wb, FWT_BACKWARD, level);
			for(int i = 0; i < n; i++){
				cx[i][0] = wb[i];
				cx[i][1] = 0;
			}
			getBlurMapWithFastInverse( ctx, cx, cx);
			for(int i = 0; i < n; i++)
				wb[i] = cx[i][0]/norm;
			ctx->fwtd3->trafo(wb,FWT_FORWARD, ll);
			for(int i = ctx->fwtd3->hb->wavlevels[ll][0];
				i < ctx->fwtd3->hb->wavlevels[ll][1]; i++)
				xx[i] -= wb[i];
		}
	}

	ierr = VecRestoreArray(x, &vx); CHKERRQ(ierr);
	ierr = VecRestoreArray(b, &vb); CHKERRQ(ierr);

	fftw_free(cx);
	delete[] xx;
	delete[] wb;
	VecDestroy(&rhs);
	VecDestroy(&sol);
	printf("PMult done.\n");///
	PetscFunctionReturn(0);
}


PetscErrorCode sparseMult(Mat mD, Vec x, Vec b){
//	printf("sparseMult: Calling\n");
	PetscErrorCode ierr;
	double* vx;
	double* vb;
	int nx;
	void* t;
	Context* ctx;
	MatShellGetContext(mD,&t);
	ctx = (Context*)t;
	ierr = VecGetArray(x, &vx); CHKERRQ(ierr);
	ierr = VecGetArray(b, &vb); CHKERRQ(ierr);
	ierr = VecGetSize(x, &nx); CHKERRQ(ierr);
	int blk_n = ctx->db->blk_size[ctx->level];
	double* xx = &vx[ctx->fwtd3->hb->wavlevels[ctx->level][0]];
	double* yy = &vb[ctx->fwtd3->hb->wavlevels[ctx->level][0]];
	for(int ii = 0; ii < blk_n; ii++) yy[ii] = 0;
	cs* A = ctx->db->sparse_blk[ctx->level - ctx->db->sparse_level];
	cs_gaxpy(A,xx,yy);
	ierr = VecRestoreArray(x, &vx); CHKERRQ(ierr);
	ierr = VecRestoreArray(b, &vb); CHKERRQ(ierr);
//	printf("sparseMult: done.\n");
	PetscFunctionReturn(0);
}

// Evaluate \Phi_i' A \Phi_j for all \Phi_i in the list
void FMMevalHB(HB* hb, int* si, int siLength, int sj, double* eval, int* temp, Neighbor3* nei, double* dist, 
	       double* rout, double* site){
//temp: a temporary memory of length 2*N for sorting the support
//temp[0:N-1] union of support of si's
//temp[N : 2N-1] pointers from node list to support index

        // For Matern kernel
	int sourcethreshold = 900;
        int targetthreshold = 8000;
	
	double t0,t1;
	int N = hb->N;
	int templength = 0;
        // printf("FMMevalHB: Merging support... for %d nodes\n",siLength);///
	for(int i = 0; i < siLength; i++){
		hb->mergeSupport(temp,templength,hb->msupport[si[i]],hb->mlength[si[i]]);
	}
	int* tempindex = &temp[N];
	for(int i = 0; i < templength; i++)
		tempindex[temp[i]] = i;
	int mlen = hb->mlength[sj];
	int* msupp = hb->msupport[sj];
	double* mask = hb->mask[sj];
	
	if(mlen > sourcethreshold && templength > targetthreshold) {
	//Initialize FMM
	// printf("FMMevalHB: Initializing FMM...\n");///
	t0 = dsecnd();
	DblNumMat* src_centers = new DblNumMat(3,mlen);
	DblNumMat* trg_centers = new DblNumMat(3,templength);
	for(int ii = 0; ii < 3; ii++){
		for(int jj = 0; jj < mlen; jj++)
			(*(src_centers))(ii,jj) = nei->coord[msupp[jj]][ii];
		for(int jj = 0; jj < templength; jj++)
			(*(trg_centers))(ii,jj) = nei->coord[temp[jj]][ii];
	}
	map<string,string> opts;
	char numchar[100];
	optionsCreate(hb->options_file, opts);
	
	// vector<double> tmp(1);  // Original
	vector<double> tmp(2);  // Changed here
         
  	map<string,string>::iterator mi;
	mi = opts.find("-alpha"); assert(mi!=opts.end());
	
	{ istringstream ss((*mi).second);  ss>>tmp[0]; }
	mi = opts.find("-kt"); assert(mi!=opts.end());
	int kt; // = 901; //r^alpha
	{ istringstream ss((*mi).second);  ss>>kt; }
        
	mi = opts.find("-rho"); assert(mi!=opts.end());
	{ istringstream ss((*mi).second);  ss>>tmp[1]; }
	
	Kernel3d* knl = new Kernel3d(kt,tmp);	
	FMM3d* fmm2 = new FMM3d("fmm3d_");
	fmm2->srcpos() = src_centers;
	fmm2->srcnor() = src_centers;
	fmm2->trgpos() = trg_centers;
	fmm2->ctr() = Point3(0,0,0);
	fmm2->rootlvl() = 0;
	fmm2->knl() = *(knl);
	sprintf(numchar,"%d",mlen);
	opts["-numsrc"] = string(numchar);
	sprintf(numchar,"%d",templength);
	opts["-numtrg"] = string(numchar);
	// printf("Initializing FMM for %s(%d) points.\n",numchar,mlen);///
	iC(fmm2->setup(opts));
	DblNumVec* srcden = new DblNumVec(mlen);
	for(int j = 0; j < mlen; j++)
	  {
	    (*srcden)(j) = mask[j];
	  }

	DblNumVec* targval = new DblNumVec(templength);	
	fmm2->eval(*srcden, *targval);

	// Evaluate the inner product
	for(int i = 0; i < siLength; i++){
		mlen = hb->mlength[si[i]];
		msupp = hb->msupport[si[i]];
		mask = hb->mask[si[i]];
		eval[i] = 0;
		for(int j = 0; j < mlen; j++){
			eval[i] += mask[j]*(*targval)(tempindex[msupp[j]]);
		}
	}
	t1 = dsecnd();
	//	printf("FMMevalHB: Done.\n");///
	delete fmm2;
	delete src_centers;
	delete trg_centers;
	delete targval;
	delete srcden;
	delete knl;
	}else{
	// jmask holds the result of A phi_j 
	double* jmask = new double[templength];
	double a2  = hb->kernelParameter * hb->kernelParameter;	
	double mid = 1e-8;
	// double K0;
	
	double t0start, t1end;
	double tdist = 0;
	double tamd = 0;
	double tout = 0;
	double tdist2 = 0;
        double dummy;
	double *x;
	double *y;
	double distvar;
        double rdist;
	int r;

	// MKL Spline Setup
	MKL_INT nsite;                      // number of interpolation sites
	MKL_INT sitehint;                   // additional info about interpolation
	MKL_INT ndorder;                    // size of array describing derivative
	MKL_INT rhint;                      // interpolation results storage format
	MKL_INT *cell_idx;                  // indices of cells containing
	MKL_INT dorder[] = {1};             // spline values will be computed
	// No cell indices are computed
	cell_idx   = 0;
	
	// Define NSITE
	nsite = templength * mlen;
	
	double *datahint;                    // additional info about structure
	int counter = 0;
	int status;

	// Parameters describing interpolation sites 
	sitehint   = DF_NON_UNIFORM_PARTITION;
	
	// Additional info about structure of arrays x and y 
	datahint   = 0;
	
	// Parameter describing interpolation results storage 
	rhint      = DF_MATRIX_STORAGE_COLS;
	
	// Parameter describing size of array for derivative orders
	ndorder    = NDORDER;
	
	// Compute kernels
	for(int i = 0; i < templength; i++){
		jmask[i] = 0;		
		if(hb->kernelType == 901) //r
		  for(int j = 0; j < mlen; j++)
			jmask[i] += sqrt(dist2(nei->coord[temp[i]],
				nei->coord[msupp[j]])) * mask[j];
		else if(hb->kernelType == 902) { //sqrt(r^2 + a^2)
		  x = nei->coord[temp[i]];
		  for(int j = 0; j < mlen; j++) {
		    y = nei->coord[msupp[j]];
		  dist[counter] = (x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2]) + a2;
		    counter++;
		  }
		}
		else if(hb->kernelType == 903) { // (r^2 + a^2)^(-1/2)
		  x = nei->coord[temp[i]];
		  for(int j = 0; j < mlen; j++) {
		    y = nei->coord[msupp[j]];
		dist[counter] = (x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2]) + a2;
		    counter++;
		  }
		}
		else if(hb->kernelType == 905){// exp(-r) Need to implement
		  x = nei->coord[temp[i]];
		}
		else if(hb->kernelType == 906){ // Matern approximation with Hermite spline
		  // Compute distance vector without sqrt
		  x = nei->coord[temp[i]];
		  for(int j = 0; j < mlen; j++){
		    y = nei->coord[msupp[j]];
		    dist[counter] = (x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2]);
		    counter++;
		  }
		}
		else if(hb->kernelType == 907) { //biharmonic r
		  x = nei->coord[temp[i]];
		  for(int j = 0; j < mlen; j++) {
		    y = nei->coord[msupp[j]];
		    dist[counter] = (x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2]);
		    counter++;
		  }
		}
		else{
			printf("Undefined kernel type.\n");
			exit(-1);
		}
	}

	// Finish sqrt(a2 + r^2) kernel
	if(hb->kernelType == 902){
	  vdSqrt(nsite,dist,site); 
	  counter = 0;
	  for(int i = 0; i < templength; i++){
	    for(int j = 0; j < mlen; j++){
	      jmask[i] += mask[j] * site[counter];
	      counter++;
	    }
	  }
	}

	// Finish invsqrt(a2 + r^2) kernel
	if(hb->kernelType == 903){
	  vdInvSqrt(nsite,dist,site); 
	  counter = 0;
	  for(int i = 0; i < templength; i++){
	    for(int j = 0; j < mlen; j++){
	      jmask[i] += mask[j] * site[counter];
	      counter++;
	    }
	  }
	}
	
	// Finish sqrt(a2) kernel
	if(hb->kernelType == 907){
	  vdSqrt(nsite,dist,site); 
	  counter = 0;
	  for(int i = 0; i < templength; i++){
	    for(int j = 0; j < mlen; j++){
	      jmask[i] += mask[j] * site[counter];
	      counter++;
	    }	    
	  }
	}

	// Finish MATERN kernel
	if(hb->kernelType == 906){
	  vdSqrt(nsite,dist,site);
	  status = dfdInterpolate1D(task, DF_INTERP, DF_METHOD_PP,
				    nsite, site, sitehint, ndorder,
				    dorder, datahint, rout, rhint, cell_idx );
	  counter = 0;
	  for(int i = 0; i < templength; i++){
	    for(int j = 0; j < mlen; j++){
	      if(site[counter] <= mid) 
		jmask[i] += mask[j] * K0;
	      else{
		jmask[i] += mask[j] * rout[counter];}
	      counter++;
	    }
	  }
	}

	// Now evaluating the dot product
	for(int i = 0; i < siLength; i++){
		mlen = hb->mlength[si[i]];
		msupp = hb->msupport[si[i]];
		mask = hb->mask[si[i]];
		eval[i] = 0;
		for(int j = 0; j < mlen; j++){
			eval[i] += mask[j]*jmask[tempindex[msupp[j]]];
		}
	}
	

	}
}

PetscErrorCode KSPMonitorResidual(KSP ksp, PetscInt its, PetscReal rnorm, void *ctx)
{
	Vec x, b;
	Mat A;
	PetscErrorCode ierr;

	//printf("KSPMonitorResidual: KSPGetOperators\n");///
	ierr = KSPGetOperators( ksp, &A, 0, 0 ); CHKERRQ(ierr);
	ierr = MatGetVecs(A, &x, &b);
	ierr = KSPBuildSolution(ksp, PETSC_NULL, &x); CHKERRQ(ierr);
	ierr = MatMult(A,x,b);

	PetscReal tol, norm_res, norm_rhs;
	PetscScalar neg_one = -1.0;
	Vec residual;
	//printf("KSPMonitorResidual: VecCopy\n");///
	ierr = VecDuplicate(*(((Context*)ctx)->rhs),&residual); CHKERRQ(ierr);
	ierr = VecCopy(*(((Context*)ctx)->rhs),residual); CHKERRQ(ierr);
	ierr = VecNorm(residual,NORM_2, &norm_rhs); CHKERRQ(ierr);
	//printf("KSPMonitorResidual: rhs norm is %g\n",norm_rhs);///
	//printf("KSPMonitorResidual: VecAXPY\n");///
	ierr = VecAXPY(residual,neg_one,b); CHKERRQ(ierr);
	//printf("KSPMonitorResidual: VecNorm\n");///
	ierr = VecNorm(residual,NORM_2, &norm_res); CHKERRQ(ierr);
	tol = norm_res/norm_rhs;
	ierr = PetscPrintf(PETSC_COMM_SELF,"Residual at iteration %D is %G.\n",its,tol);
	((Context*)ctx)->actualResidual[its] = tol;

	return(0);
}
