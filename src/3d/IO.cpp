#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include "fftw3.h"
#include "IO.h"

#include "ByteSwapping.h"
#include "matrixio.h"
#include "cs.h"
#include "linearSolver.h"

bool makeVariablesBigEndian( float* data, int size )
{
	if( !data || (size<1)) return false;
	
	if (data && isLittleEndian()) 
		swapByteOrder(data, size* size* size);
	return true;
}

bool createRawIVHeader(char *header, float* real_data, float span[3], float orig[3], int size, float min[3], float max[3] )
{
	if (!real_data)
		return false;
	
	unsigned int numverts, numcells;
	
	numverts = size* size* size;
	numcells = (size-1)*(size-1)*(size-1);
	
	// minX
	memcpy(header, &min[0], sizeof(float));
	header += 4;
	// minY
	memcpy(header, &min[1], sizeof(float));
	header += 4;
	// minZ
	memcpy(header, &min[2], sizeof(float));
	header += 4;
	
	// maxX
	memcpy(header, &max[0], sizeof(float));
	header += 4;
	// maxY
	memcpy(header, &max[1], sizeof(float));
	header += 4;
	// maxZ
	memcpy(header, &max[2], sizeof(float));
	header += 4;
	
	// numverts 
	memcpy(header, &numverts, sizeof(numverts));
	header += 4;
	
	// numcells 
	memcpy(header, &numcells, sizeof(numcells));
	header += 4;
	
	// dimX
	memcpy(header, &size, sizeof(unsigned int));
	header += 4;
	// dimY
	memcpy(header, &size, sizeof(unsigned int));
	header += 4;
	// dimZ
	memcpy(header, &size, sizeof(unsigned int));
	header += 4;
	
	// originX
	memcpy(header, &orig[0], sizeof(float));
	header += 4;
	// originY
	memcpy(header, &orig[1], sizeof(float));
	header += 4;
	// originZ
	memcpy(header, &orig[2], sizeof(float));
	header += 4;
	
	// spanX
	memcpy(header, &span[0], sizeof(float));
	header += 4;
	// spanY
	memcpy(header, &span[1], sizeof(float));
	header += 4;
	// spanZ 
	memcpy(header, &span[2], sizeof(float));
	
	header -= 64;
	if (isLittleEndian())
		swapByteOrder((unsigned int *)header, 17);
	
	return true;
}

bool print_im( fftw_complex* data, int size, char* fileName, float min[3], float max[3] )
{
	if( !data || (size<1) || !fileName ) return false;
	
	FILE *fp = fopen(fileName, "wb");
	if( !fp ) return false;
	
	char header[68];
	
	float* real_data = new float[size* size* size];
	
	int i, j, k;
	int norm = 1;

	for( i=0; i<size; i++)
    {
		for( j=0; j<size; j++ )
		{
			for( k=0; k<size; k++ )
			{
				real_data[i*size*size + j*size + k] = (float)(data[i*size*size + j*size + k][1])/norm;
			}
		}
    }
		
	{
		int i;
		double max_re = data[0][0]/norm;
		double min_re = data[0][0]/norm;
		double max_im = data[0][1]/norm;
		double min_im = data[0][1]/norm;
		
		for( i=0; i<size*size*size; i++ )
		{
			if( max_re  < data[i][0]/norm ) max_re = data[i][0]/norm;
			if( min_re > data[i][0]/norm ) min_re =data[i][0]/norm;
			
			if( max_im < data[i][1]/norm ) max_im = data[i][1]/norm;
			if( min_im > data[i][1]/norm ) min_im = data[i][1]/norm;
		}
		printf("In the end, the min max are = [%lf %lf]-[%lf %lf]\n", min_re, min_im, max_re, max_im );
	}
	
	// make the data big endian
	if( !makeVariablesBigEndian(real_data, size) )
    {
		printf("Could not make big endian\n");
		delete []real_data;
		fclose( fp );
		return false;
    }
	
	// create the header
	float span[3];
	span[0] = (max[0] - min[0]) / (float)size;
	span[1] = (max[1] - min[1]) / (float)size;
	span[2] = (max[2] - min[2]) / (float)size;
	
	printf("\tMin Max: [%f %f %f] - [%f %f %f]\n", min[0], min[1],min[2], max[0], max[1], max[2]);
	if( !createRawIVHeader(header, real_data, span, min, size, min, max ) )
    {
		printf("Could not create the header\n");
		delete []real_data;
		fclose(fp);
		return false;
    }
	
	// write the header
	fwrite(header, sizeof(header), 1, fp);
	
	
	// write the data
	fwrite(real_data, sizeof(float), size*size*size, fp);
	
	// clean up
	fclose(fp);
	delete []real_data;
	
	return true;
}

bool print_re( fftw_complex* data, int size, char* fileName, float min[3], float max[3] )
{
	if( !data || (size<1) || !fileName ) return false;
	
	FILE *fp = fopen(fileName, "wb");
	if( !fp ) return false;
	
	char header[68];
	
	float* real_data = new float[size* size* size];
	
	int i, j, k;
	int norm = 1.0;

	for( i=0; i<size; i++)
    {
		for( j=0; j<size; j++ )
		{
			for( k=0; k<size; k++ )
			{
				real_data[i*size*size + j*size + k] = (float)(data[i*size*size + j*size + k][0])/norm;
			}
		}
    }
	
	
	{
		int i;
		double max_re = data[0][0]/norm;
		double min_re = data[0][0]/norm;
		double max_im = data[0][1]/norm;
		double min_im = data[0][1]/norm;
		
		for( i=0; i<size*size*size; i++ )
		{
			if( max_re  < data[i][0]/norm ) max_re = data[i][0]/norm;
			if( min_re > data[i][0]/norm ) min_re =data[i][0]/norm;
			
			if( max_im < data[i][1]/norm ) max_im = data[i][1]/norm;
			if( min_im > data[i][1]/norm ) min_im = data[i][1]/norm;
		}
		printf("In the end, the min max are = [%lf %lf]-[%lf %lf]\n", min_re, min_im, max_re, max_im );
	}
	
	// make the data big endian
	if( !makeVariablesBigEndian(real_data, size) )
    {
		printf("Could not make big endian\n");
		delete []real_data;
		fclose( fp );
		return false;
    }
	
	// create the header
	float span[3];
	span[0] = (max[0] - min[0]) / (float)size;
	span[1] = (max[1] - min[1]) / (float)size;
	span[2] = (max[2] - min[2]) / (float)size;
	
	printf("\tMin Max: [%f %f %f] - [%f %f %f]\n", min[0], min[1],min[2], max[0], max[1], max[2]);
	if( !createRawIVHeader(header, real_data, span, min, size, min, max ) )
    {
		printf("Could not create the header\n");
		delete []real_data;
		fclose(fp);
		return false;
    }
	
	// write the header
	fwrite(header, sizeof(header), 1, fp);
	
	
	// write the data
	fwrite(real_data, sizeof(float), size*size*size, fp);
	
	// clean up
	fclose(fp);
	delete []real_data;
	
	return true;
}

bool print( fftw_complex* data, int size, char* fileName, float min[3], float max[3] )
{
	if( !data || (size<1) || !fileName ) return false;
	
	FILE *fp = fopen(fileName, "wb");
	if( !fp ) return false;
	
	char header[68];
	
	float* real_data = new float[size* size* size];
	
	int i, j, k;
	float norm = size*size*size;
	for( i=0; i<size; i++)
    {
		for( j=0; j<size; j++ )
		{
			for( k=0; k<size; k++ )
			{
				real_data[i*size*size + j*size + k] = (float)(data[i*size*size + j*size + k][0])/norm;
			}
		}
    }
		
	{
		int i;
		double max_re = data[0][0]/norm;
		double min_re = data[0][0]/norm;
		double max_im = data[0][1]/norm;
		double min_im = data[0][1]/norm;
		
		for( i=0; i<size*size*size; i++ )
		{
			if( max_re  < data[i][0]/norm ) max_re = data[i][0]/norm;
			if( min_re > data[i][0]/norm ) min_re =data[i][0]/norm;
			
			if( max_im < data[i][1]/norm ) max_im = data[i][1]/norm;
			if( min_im > data[i][1]/norm ) min_im = data[i][1]/norm;
		}
		printf("In the end, the min max are = [%lg %lg]-[%lg %lg]\n", min_re, min_im, max_re, max_im );
	}
	
	// make the data big endian
	if( !makeVariablesBigEndian(real_data, size) )
    {
		printf("Could not make big endian\n");
		delete []real_data;
		fclose( fp );
		return false;
    }
	
	// create the header
	float span[3];
	span[0] = (max[0] - min[0]) / (float)size;
	span[1] = (max[1] - min[1]) / (float)size;
	span[2] = (max[2] - min[2]) / (float)size;
	
	printf("\tMin Max: [%f %f %f] - [%f %f %f]\n", min[0], min[1],min[2], max[0], max[1], max[2]);
	if( !createRawIVHeader(header, real_data, span, min, size, min, max ) )
    {
		printf("Could not create the header\n");
		delete []real_data;
		fclose(fp);
		return false;
    }
	
	// write the header
	fwrite(header, sizeof(header), 1, fp);
	
	
	// write the data
	fwrite(real_data, sizeof(float), size*size*size, fp);
	
	// clean up
	fclose(fp);
	delete []real_data;
	
	return true;
}

// reads only the raw part of rawiv, assumes no header.
void readRawiv( const char* fileName, fftw_complex* data, int size )
{
	FILE *fp = fopen(fileName, "rb");
	if( !fp ) { printf("Could not open file\n"); exit(0); }
	
	float* real_data = new float[size* size* size];
	fread( real_data, sizeof(float), size*size*size, fp);
	makeVariablesBigEndian(real_data, size);
	
	{
		int i;
		for( i=0; i<size*size*size; i++ )
		{
			data[i][0] = real_data[i];
			data[i][1] = 0;
		}
	}
	
	delete []real_data;
	
	fclose( fp );
}

void read_diag_blk(int argc, char* argv[], Diag_Blk* db){
	int m,n,blk_n;
	double* blk;
	double* blk_lu;
	int* ipiv;
	cs* csblk;

	int van = atoi(argv[9]);
	int pvan = atoi(argv[10]);

	int nlevel = argc - 13;
	for(int ii = 0; ii < nlevel; ii++){

		if((ii < db->sparse_level && db->blk_size[ii] < 1e4) || 
		  (van > pvan && ii == nlevel - 1)){ // coarse-lifted level
			readVector(argv[13 + ii],n,blk);
			blk_n = sqrt(n);
			blk_lu = new double[blk_n * blk_n];
			ipiv = new int[blk_n];
			for(int ii = 0; ii < n; ii++)
				blk_lu[ii] = blk[ii];
			lu_init(blk_n,blk_n,blk_lu,ipiv);
			db->dense_blk[ii] = blk;
			db->dense_lu[ii] = blk_lu;
			db->ipiv[ii] = ipiv;
			db->blk_size[ii] = blk_n;
		}else{
			readcsMatrix(argv[13+ii], csblk);
			db->sparse_blk[ii - db->sparse_level] = cs_compress(csblk);
			db->blk_size[ii] = csblk->n;
			cs_spfree(csblk);
		}
	}
}

void write_diag_blk(int argc, char* argv[], HB* hb, Diag_Blk* db){
	int m,n,blk_n;
	double* blk;
	double* blk_lu;
	cs* csblk;

	int van = atoi(argv[9]);
	int pvan = atoi(argv[10]);

	int nlevel = hb->wavMaxLevel + 1;
	for(int li = 0; li < nlevel; li++)
		if((li < db->sparse_level && db->blk_size[li] < 1e4)||
		  (van > pvan && li == nlevel - 1)){ // coarse-lifted level
			blk_n = db->blk_size[li];
			blk = db->dense_blk[li];
			writeMatrix(argv[13 + li], 1, blk_n*blk_n, &blk);
		}else{
			csblk = db->sparse_blk[li - db->sparse_level];
			writecsMatrix(argv[13 + li], csblk);
		}
}
