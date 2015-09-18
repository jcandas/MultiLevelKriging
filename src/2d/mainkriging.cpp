#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include "FastInverse.h"
#include "ByteSwapping.h"
#include "IO.h"
#include "linearSolver.h"
//jj 
#include "matrixio.h"
#include "mpi.h"

#ifndef M_PI
#define M_PI   3.14159265358979323846
#endif

#include <time.h>
#ifdef _WIN32
#include <sys/types.h>
#include <sys/timeb.h>
#else
#include <sys/time.h>
#endif

/*****************************************************************************************************************/
/*                                                                                                               */
/*  Get the current time in seconds as a double value                                                            */
/*                                                                                                               */
/*****************************************************************************************************************/
double getTime()
{
#ifdef _WIN32
	time_t ltime;
	_timeb tstruct;
	time( &ltime );
	_ftime( &tstruct );
	return (double) (ltime + 1e-3*(tstruct.millitm));
#else
    struct timeval t;
    gettimeofday( &t, NULL );  
    return (double)(t.tv_sec + 1e-6*t.tv_usec);
#endif
      }

bool getNormalizedPositions( double* rbfCenterPositions, double* xk, double* yk, double* zk, int rbfWidth, int gridSize, int numCenters, int interpFuncExtent, int numFreq )
{
	if( !rbfCenterPositions || !xk || !yk || !zk || (rbfWidth<1) || (gridSize<1) || (numCenters<1) ) return false;
	
	if( rbfWidth > 2*gridSize ) rbfWidth = gridSize; // we dont care outside this interval
	
	// get the largest side
	double minX, maxX, minY, maxY, minZ, maxZ;
	minX = maxX = rbfCenterPositions[0*0+0];
	minY = maxY = rbfCenterPositions[0*0+1];
	minZ = maxZ = rbfCenterPositions[0*0+2];
	int i;
	for( i=1; i<numCenters; i++ )
    {
		if( minX > rbfCenterPositions[i*3+0] ) minX = rbfCenterPositions[i*3+0];
		if( maxX < rbfCenterPositions[i*3+0] ) maxX = rbfCenterPositions[i*3+0];
		
		if( minY > rbfCenterPositions[i*3+1] ) minY = rbfCenterPositions[i*3+1];
		if( maxY < rbfCenterPositions[i*3+1] ) maxY = rbfCenterPositions[i*3+1];
		
		if( minZ > rbfCenterPositions[i*3+2] ) minZ = rbfCenterPositions[i*3+2];
		if( maxZ < rbfCenterPositions[i*3+2] ) maxZ = rbfCenterPositions[i*3+2];
    }
        printf("The max of the molecule was [%lf %lf %lf]\n", maxZ, maxY, maxX);
        printf("The min of the molecule was [%lf %lf %lf]\n", minZ, minY, minX);
	printf("The size of the molecule was [%lf %lf %lf]\n", maxZ-minZ, maxY-minY, maxX-minX );
	double largestWidth = maxX - minX;
	if( maxY-minY > largestWidth ) largestWidth = maxY-minY;
	if( maxZ-minZ > largestWidth ) largestWidth = maxZ-minZ;
	
	double scale = (gridSize/((double)(gridSize+2*rbfWidth))) / (double)largestWidth;

	printf("gridSize = %d, rbfWidth=%d\n", gridSize, rbfWidth);
	
	printf("Largest side was %lf and scale was %lf\n", largestWidth, scale );
	
	{
		double xTrans = - minX - (maxX-minX)/2.0;
		double yTrans = - minY - (maxY-minY)/2.0;
		double zTrans = - minZ - (maxZ-minZ)/2.0;

		double scaleFactor = scale;
		
		for( i=0; i<numCenters; i++ )
		{
			xk[i] = rbfCenterPositions[i*3+0] + xTrans; // from -w/2 to +w/2
			xk[i] = xk[i] * scaleFactor; // within -0.5 to 0.5, including appropriate 0 padding
			
			yk[i] = rbfCenterPositions[i*3+1] + yTrans; // from -w/2 to +w/2
			yk[i] = yk[i] * scaleFactor; // within -0.5 to 0.5, including appropriate 0 padding
			
			zk[i] = rbfCenterPositions[i*3+2] + zTrans; // from -w/2 to +w/2
			zk[i] = zk[i] * scaleFactor; // within -0.5 to 0.5, including appropriate 0 padding
		}
	}
	
	printf("Translated by [%lf %lf %lf]\n", -minX - (maxX-minX)/2.0, - minY - (maxY-minY)/2.0, -minZ - (maxZ-minZ)/2.0);
	printf("Scaled by %lf\n", scale );
	
	{
		double scale = 1;
		for( i=0; i<numCenters; i++ )
		{
			xk[i] = xk[i] * scale; // within -0.5 to 0.5, including appropriate 0 padding
			
			yk[i] = yk[i] * scale; // within -0.5 to 0.5, including appropriate 0 padding
			
			zk[i] = zk[i] * scale; // within -0.5 to 0.5, including appropriate 0 padding
		}
	}
	return true;
}

/*****************************************************************************************************************/
/*                                                                                                               */
/*                                                                                                               */
/*      The inputs to the program are:                                                                           */
/*         alpha - the grid sampling factor. We create a grid of size alpha*numFreq                              */
/*         interpFuncExtent - the gaussian extent                                                                */
/*         numFreq - Number of frequencies to extract                                                            */
/*         gridSize - How large a grid to embed molecule in for full resolution search                           */
/*         inputFile1 - The input file of the molecule.                                                          */
/*         out.txt - output file name                                                                            */
/*         rbfFile.txt - rbf values in 1D                                                                        */
/*         Target.txt -  Target centers (kriging)                                                                */
/*                                                                                                               */         
/*****************************************************************************************************************/
int main( int argc, char* argv[] )
{

	double alpha;
	int interpFuncExtent;
	int numFreq;
	int gridSize;

	char inputFile[256];
	char inputTargetFile[256];
	char inputCoef[256];
	char outputFileName[256];
	char rbfFile[256];
	int van;
	int pvan;
	int minlevel;  
	double tol;
	// int max_iter = 60;
	int max_iter = 10000;
	int myid;
	double t0, t1;
	int d = 0; 
	int bitmap;

	MPI_Init(0,0);	
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
if(myid == 0){
	t0 = dsecnd();
	if( argc < 15 )
    {
		printf("\n\nUsage : %s alpha interpFuncExtent numFreq gridSize inputFile1(coords) inputFile2(coefs) \
                        outputFileName kifmmoptions (ftilde+1) (f+1) GMRES_tolerance [bitmap] Target \
                        distance-criterion-cte minlevel \n",argv[0]);
		printf("   bitmap 1 = use fast diagonal computation (diagonal preconditioner)\n");
		printf("   bitmap 2 = use PMult (SSOR preconditioner)\n");
		printf("   bitmap 4 = write diagonal blocks\n");
		printf("   bitmap 8 = read in the recomputed diagonal blocks\n");
		printf("   bitmap 16 = compute actual residual for each iteration\n");
		printf("   bitmap 32 = compute sparse Kw matrix and write it to file \n");
		return -1;
    }
	
	printf("Wall-clock time used: %f seconds \n",dsecnd() - t0);
	alpha = atof( argv[1] );
	interpFuncExtent = atoi( argv[2] );
	numFreq = atoi( argv[3] );
	gridSize = atoi( argv[4] );
	strcpy(inputFile, argv[5] );
	strcpy(inputCoef, argv[6]);
	strcpy( outputFileName, argv[7] );
	strcpy( rbfFile, argv[8] );
	van = atoi( argv[9] );
	pvan = atoi( argv[10] );
	sscanf(argv[11],"%lf",&tol);
	bitmap = 3;
	if(argc > 12) bitmap = atoi(argv[12]);
        strcpy(inputTargetFile, argv[13] );	

	// Distance criterion neighborhood length
	d = atoi( argv[14] );
        minlevel = atoi( argv[15] );

	printf("Wall-clock time used: %f seconds \n",dsecnd() - t0);

	// printf("Reading the rbf values\n");
	int rbfWidth = 0;
	
	printf("Wall-clock time used: %f seconds \n",dsecnd() - t0);
	fftw_complex* RBF3D = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * gridSize * gridSize * gridSize);
	printf("Wall-clock time (test) used: %f seconds \n",dsecnd() - t0);

	{//jj
		int rbfSize = gridSize ;
		readRawiv(rbfFile,RBF3D,rbfSize);
		rbfWidth = gridSize /2;
	}
				
	printf("Wall-clock time used: %f seconds \n",dsecnd() - t0);
	//printf("Reading the centers\n");
	fftw_complex* rbfCenterValues = 0;
	double* rhs = NULL;
	double* rbfCenterPositions = 0;
	double** coefs;
	int numCenters;  //numCenters is usually N+2
	float min[3];
	float max[3];
	int m,n;
	// read in the centers values and positions of source points
	{
		FILE* fp = fopen( inputFile, "r" );
		if( !fp ) 
		{
			printf("Could not open file %s\n", inputFile );
			fftw_free( RBF3D );
			return -1;
		}
		if( !fscanf( fp, "%d\n", &numCenters ) )
		{
			printf("Format error in file %s line 0\n", inputFile );
			fftw_free( RBF3D );
			return -1;
		}
		rbfCenterValues = new fftw_complex[numCenters];
		rhs = new double[numCenters];
		rbfCenterPositions = new double[numCenters*3];
		
		int i;
		{
			if( fscanf( fp, "%lf %lf %lf %lf %lf\n", &(rbfCenterValues[0][0]), &(rbfCenterValues[0][1]), 
                          &(rbfCenterPositions[0*3+0]), &(rbfCenterPositions[0*3+1]), &(rbfCenterPositions[0*3+2]) ) != 5 )
			{
				printf("Format error in file %s line %d\n", inputFile, 0+1 );
				fftw_free( RBF3D );
				delete [] rbfCenterValues;
				delete [] rbfCenterPositions;
				delete [] rhs;
				return -1;
			}
			min[0] = max[0] = rbfCenterPositions[0*3+0];
			min[1] = max[1] = rbfCenterPositions[0*3+1];
			min[2] = max[2] = rbfCenterPositions[0*3+2];
		}
		for( i=1; i<numCenters; i++ )
		{
			if( fscanf( fp, "%lf %lf %lf %lf %lf\n", &(rbfCenterValues[0][0]), &(rbfCenterValues[0][1]), 
                          &(rbfCenterPositions[i*3+0]), &(rbfCenterPositions[i*3+1]), &(rbfCenterPositions[i*3+2]) ) != 5 )
			{
				printf("Format error in file %s line %d\n", inputFile, i+1 );
				fftw_free( RBF3D );
				delete [] rbfCenterValues;
				delete [] rbfCenterPositions;
				delete [] rhs;
				return -1;
			}
			if( rbfCenterPositions[i*3+0] < min[0] ) min[0] = rbfCenterPositions[i*3+0];
			if( rbfCenterPositions[i*3+1] < min[1] ) min[1] = rbfCenterPositions[i*3+1];
			if( rbfCenterPositions[i*3+2] < min[2] ) min[2] = rbfCenterPositions[i*3+2];
			
			if( rbfCenterPositions[i*3+0] > max[0] ) max[0] = rbfCenterPositions[i*3+0];
			if( rbfCenterPositions[i*3+1] > max[1] ) max[1] = rbfCenterPositions[i*3+1];
			if( rbfCenterPositions[i*3+2] > max[2] ) max[2] = rbfCenterPositions[i*3+2];
			
		}
		printf("min , max : [%f %f %f], [%f %f %f]\n",min[0],min[1],min[2],max[0],max[1],max[2]);///
		fclose( fp );
		readMatrix(inputCoef,m,n,coefs);
		printf("Input Coefficients m = %d n = %d \n", m, n);
	}


	// Read position of Target Points
	double* rbfTargetCenterPositions = 0;
        int numTargetCenters;  //numCenters is usually N+2

	
	{
		FILE* fp = fopen( inputTargetFile, "r" );
		if( !fp ) 
		{
			printf("Could not open file %s\n", inputTargetFile );
			fftw_free( RBF3D );
			delete [] rbfCenterValues;
			delete [] rbfCenterPositions;
			delete [] rhs;
			return -1;
		}
		if( !fscanf( fp, "%d\n", &numTargetCenters ) )
		{
			printf("Format error in file %s line 0\n", inputTargetFile );
			fftw_free( RBF3D );
			delete [] rbfCenterValues;
			delete [] rbfCenterPositions;
			delete [] rhs;
			return -1;
		}
		rbfTargetCenterPositions = new double[numTargetCenters*3];
		
		{
			if( fscanf( fp, "%lf %lf %lf\n",
                          &(rbfTargetCenterPositions[0*3+0]), &(rbfTargetCenterPositions[0*3+1]), 
				    &(rbfTargetCenterPositions[0*3+2]) ) != 3 )
			{
				printf("Format error in file %s line %d\n", inputFile, 0+1 );
				fftw_free( RBF3D );
				
			delete [] rbfCenterValues;
			delete [] rbfCenterPositions;
			delete [] rhs;
			delete [] rbfTargetCenterPositions;
				return -1;
			}
		}
                int i;
		for( i=1; i<numTargetCenters; i++ )
		{
			if( fscanf( fp, "%lf %lf %lf\n", 
                          &(rbfTargetCenterPositions[i*3+0]), &(rbfTargetCenterPositions[i*3+1]), 
				    &(rbfTargetCenterPositions[i*3+2]) ) != 3 )
			{
				printf("Format error in file %s line %d\n", inputFile, i+1 );
				fftw_free( RBF3D );
			        delete [] rbfCenterValues;
				delete [] rbfCenterPositions;
				delete [] rhs;
				delete [] rbfTargetCenterPositions;
				return -1;
			}
		}

		fclose( fp );
        }


	printf("Number of Targets centers: %d \n", numTargetCenters);
       	
	printf("Finished inputing source and target files. Wall-clock time used: %f seconds \n",dsecnd() - t0);
       	
	// Change to multiple dimensions
	int nnei = van*(van+1)/2; //2D
	
	// initializing unknown x	
	double** x = new double*[m];
	for(int i = 0; i < m; i++ ){
		x[i] = new double[n+nnei];
		for(int j = 0; j < n+nnei; j++)
			x[i][j] = 0.0;
		}
	
	printf("Wall-clock time used: %f seconds \n",dsecnd() - t0);
	KSP solver;
	Mat mA, mD;
	PetscErrorCode ierr;
	PetscReal* rhistory;
	PetscInt rhistorySize;
	rhistory = new PetscReal[max_iter];
	printf("argc = %d, n = %d, numCenters = %d, gridSize = %d, tol = %g, max_iter = %d, van = %d, pvan = %d\n",
		argc,n,numCenters,gridSize,tol,max_iter,van,pvan);///
	printf("RBF3D[%d][%d][%d] = %g\n",0,0,0,RBF3D[0][0]);///
	printf("RBF3D[%d][%d][%d] = %g\n",gridSize-1,gridSize-1,gridSize-1,RBF3D[gridSize*gridSize*gridSize-1][0]);///
	printf("Wall-clock time used: %f seconds \n",dsecnd() - t0);


	// Load RHS
	for(int j = 0; j < n; j++){
	  rhs[j] = coefs[0][j];
	}

	ierr = gmres_init(argc, argv, n, numCenters, gridSize, alpha, d, interpFuncExtent, 
			  tol, max_iter, van, pvan, rbfCenterPositions,  RBF3D, solver, mA, mD, minlevel, rhs);
	t1 = dsecnd();
	printf("Initialization wall-clock time used: %f seconds \n",t1 - t0);

	if(!(bitmap & 32)){
    if(!(bitmap & 4) || (bitmap & 8)){ //If not just computing the diag blks
	for (int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			rhs[j] = coefs[i][j];
		}
		printf("Ready for GMRES outer iterations.\n");///
		ierr = KSPSetResidualHistory(solver, rhistory, max_iter,PETSC_TRUE); CHKERRQ(ierr);
		ierr = gmres(n,solver, mA, rhs, x[i], tol, max_iter,rbfCenterPositions,numCenters,
			     rbfTargetCenterPositions,numTargetCenters);
                CHKERRQ(ierr);
		printf("GMRES converged in %d iterations with residue %g.\n",max_iter,tol);
		printf("Wall-clock time used: %f seconds \n",dsecnd() - t0);
		
		void* ctx;
		ierr = MatShellGetContext(mA, &ctx); CHKERRQ(ierr);
		printf("GMRES residual history: \n");
		for(int j = 1; j < max_iter + 2; j++)
		{
			printf("Iter %d: %14g",j-1,rhistory[j]/rhistory[0]);
			if(((Context*)ctx)->printResidual)
				printf("      Actual: %14g\n", ((Context*)ctx)->actualResidual[j-1]);
			else
				printf("\n");		
		}
	}
	printf("Wall-clock time used: %f seconds \n",dsecnd() - t0);
	
	printf("Writing out solutions to %s\n",outputFileName);
	writeMatrix(outputFileName, m, n+nnei, x);
    }
	}
	ierr = gmres_finalize(solver, mA, mD);
	
	delete[] rhistory;
	for(int i = 0; i < m; i++)
		delete[] coefs[i];
	delete[] coefs;
	delete[] x;
	
	printf("Freeing up RBF3D..\n");	
	fftw_free( RBF3D );
	delete [] rbfCenterValues;
	delete [] rbfCenterPositions;
        delete [] rbfTargetCenterPositions; // Delete target points

	printf("Freeing up Everything else..\n");	
	t1 = dsecnd();
	printf("Total wall-clock time used: %f seconds \n",t1 - t0);
}
	printf("%d: waiting on the barrier..\n",myid);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	printf("%d: done and waiting..\n",myid);
	return 0;
}
