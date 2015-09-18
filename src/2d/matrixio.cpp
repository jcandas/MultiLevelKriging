#include <iostream>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include "matrixio.h"
using namespace std;

void readTexMatrix(char *input,int m, int n, double **& data)
{
        ifstream in(input);
        data = new double*[m];
        for(int i = 0; i < m; i++)
                data[i] = new double[n];
        for(int i = 0; i < m; i++)
        for(int j = 0; j < n; j++)
                in >>data[i][j];
}


void readMatrix(const char* fname, int& m, int& n, double**& data){
	struct stat filestat;
	//int um,un;
	if (stat(fname, &filestat) < 0) {
		fprintf(stderr, "cannot find data file %s\n", fname);
		return;
		exit(-1);
	}
	int sz = filestat.st_size;
	//	cout << "opening file.." <<'\n';///
	FILE *fp = fopen(fname, "rb" );
	//	cout << "file opened" << '\n';///
	if (fp == NULL) {
		fprintf(stderr, "cannot open %s\n", fname);
		exit(-1);
	}
	fread(&m,sizeof(int),1,fp);
	fread(&n,sizeof(int),1,fp);
	data=new double*[m];
	for(int ii = 0; ii < m; ii++)
		data[ii]=new double[n];
	for(int jj = 0; jj < n; jj++)
		for(int ii = 0; ii < m; ii++)
			fread(&data[ii][jj], sizeof(double), 1, fp);
	fclose(fp);
}

void readVector(const char* fname, int& m, double*& data){
	double** pdata;
	int n;
	readMatrix(fname,n,m,pdata);
	data = pdata[0];
	delete[] pdata;
}
void writeMatrix(const char* fname, int m, int n, double** data){
	struct stat filestat;

	FILE *fp = fopen(fname, "wb");
	if (fp == NULL) {
		fprintf(stderr, "cannot open %s\n", fname);
		exit(-1);
	}
	fwrite(&m, sizeof(int), 1, fp);
	fwrite(&n, sizeof(int), 1, fp);
	for(int jj = 0; jj < n; jj++)
		for(int ii = 0; ii < m; ii++)
			fwrite(&data[ii][jj], sizeof(double), 1, fp);
	fclose(fp);
}

void writeIntMatrix(const char* fname, int m, int n, int** data){
	struct stat filestat;
	double t;

	FILE *fp = fopen(fname, "wb");
	if (fp == NULL) {
		fprintf(stderr, "cannot open %s\n", fname);
		exit(-1);
	}
	fwrite(&m, sizeof(int), 1, fp);
	fwrite(&n, sizeof(int), 1, fp);
	for(int jj = 0; jj < n; jj++)
		for(int ii = 0; ii < m; ii++){
			t = data[ii][jj];
			fwrite(&t, sizeof(double), 1, fp);
		}
	fclose(fp);
}

void writeTripletMatrix(const char* fname, int m, int n, int* length, int** supp,double** mask){
	int sparse_max = 145000000;
	if(m*n < sparse_max) sparse_max = m*n;
	double** sparse_mat = new double*[3];
	double* sparse_i = new double[sparse_max];
	double* sparse_j = new double[sparse_max];
	double* sparse_v = new double[sparse_max];
	sparse_mat[0] = sparse_i;
	sparse_mat[1] = sparse_j;
	sparse_mat[2] = sparse_v;
        int index = 0;
        for(int i = 0; i < m; i++){
                for(int j = 0; j < length[i]; j++){
                        sparse_i[index] = i+1;
                        sparse_j[index] = supp[i][j]+1;
                        if(mask) sparse_v[index] = mask[i][j];
			else sparse_v[index] = 1;
                        index++;
                        if(index >= sparse_max){
                                printf("Sparse_max too small!\n");
				exit(-1);
			}
                }
        }
	// printf("Writing out the sparse matrix (3 x %d) to %s\n",index,outputFileName);///
        writeMatrix(fname, 3, index, sparse_mat);
	delete[] sparse_i;
	delete[] sparse_j;
	delete[] sparse_v;
	delete[] sparse_mat;
}

void writeIntVector(const char* fname, int m, int* data){
	struct stat filestat;
	FILE *fp = fopen(fname, "wb");
	int n=1;
	double t;
	if (fp == NULL) {
		fprintf(stderr, "cannot open %s\n", fname);
		exit(-1);
	}
	fwrite(&m, sizeof(int), 1, fp);
	fwrite(&n, sizeof(int), 1, fp);
	for(int jj = 0; jj < n; jj++)
		for(int ii = 0; ii < m; ii++){
			t=data[ii];
			fwrite(&t, sizeof(double), 1, fp);
		}
	fclose(fp);
}


void genA(double** coord,const int N, double**& A){
	double r2;
	A = new double*[N];
	for(int ii = 0; ii < N; ii++)
		A[ii]=new double[N];
	for(int ii = 0; ii < N; ii++)
		A[ii][ii] = 0;
		
	for(int ii = 0; ii < N-1; ii++)
		for(int jj = ii+1; jj < N; jj++){
			r2 = dist2(coord[ii],coord[jj]);
			A[ii][jj] = r2*log(r2)/2;
			A[jj][ii] = A[ii][jj];
		}
}


double dist2(double* x,double* y){
	        return (x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2]);
}

double sub(double* x,double* y, int r){
  return (x[r]-y[r]);
}




void writecsMatrix(const char* fname, cs* A){
	FILE *fp = fopen(fname, "wb");
	int m = 3;
	int n;
	if (CS_CSC(A)) 
		n = A->p[A->n];
	if (CS_TRIPLET(A)) 
		n =A->nz;
	double t;
	if (fp == NULL) {
		fprintf(stderr, "cannot open %s\n", fname);
		exit(-1);
	}
	fwrite(&m, sizeof(int), 1, fp);
	fwrite(&n, sizeof(int), 1, fp);
	//for compressed cs matrix
	if(CS_CSC(A)){
		for(int jj = 0; jj < A->n; jj++)
			for(int ii = A->p[jj]; ii < A->p[jj+1]; ii++){
				t= A->i[ii]+1;
				fwrite(&t, sizeof(double), 1, fp);
				t= jj+1;
				fwrite(&t, sizeof(double), 1, fp);
				t= A->x[ii];
				fwrite(&t, sizeof(double), 1, fp);
			}
	}
	if(CS_TRIPLET(A)){
		for(int ii = 0; ii < A->nz; ii++){
			t= A->i[ii]+1;
			fwrite(&t, sizeof(double), 1, fp);
			t= A->p[ii]+1;
			fwrite(&t, sizeof(double), 1, fp);
			t= A->x[ii];
			fwrite(&t, sizeof(double), 1, fp);
		}
	}
	fclose(fp);
}

void readcsMatrix(const char* fname, cs* &A){
	printf("Entering readcsMatrix ...%s\n",fname);
	int m,n;
	int i,j;
	double x;

	double data[3];
	struct stat filestat;
	//int um,un;
	if (stat(fname, &filestat) < 0) {
		fprintf(stderr, "cannot find data file %s\n", fname);
		return;
		exit(-1);
	}
	int sz = filestat.st_size;
	//	cout << "opening file.." <<'\n';///
	FILE *fp = fopen(fname, "rb" );
	//	cout << "file opened" << '\n';///
	if (fp == NULL) {
		fprintf(stderr, "cannot open %s\n", fname);
		exit(-1);
	}
	fread(&m,sizeof(int),1,fp);
	fread(&n,sizeof(int),1,fp);
	A = cs_spalloc(0,0,1,1,1);
	for(int jj = 0; jj < n; jj++){
		fread(&data[0], sizeof(double), 1, fp);
		fread(&data[1], sizeof(double), 1, fp);
		fread(&data[2], sizeof(double), 1, fp);
		i = data[0] - 1;
		j = data[1] - 1;
		x = data[2];
		if (! cs_entry(A, i, j, x)) {
			cs_spfree(A);
			exit(-1); // out of memory
		}
	}
	fclose(fp);
	// cout << "file read" << '\n';///
	printf("readcsMatrix: Done.\n");
}
