#ifndef MATRIXIO_H
#define MATRIXIO_H
#include "cs.h"
void readTexMatrix(char *input,int m, int n, double **& data);
void readMatrix(const char* fname, int& m, int& n, double**& data);
void readVector(const char* fname, int& m, double*& data);
void readcsMatrix(const char* fname,cs* &A);
void writeMatrix(const char* fname, int m, int n, double** data);
void writeIntMatrix(const char* fname, int m, int n, int** data);
void writeTripletMatrix(const char* fname, int m, int n, int* length, int** supp,double** mask);
void writeIntVector(const char* fname, int m,  int* data);
void genA(double** coord,const int N, double**& A);
double dist2(double* x,double* y);
double sub(double* x,double* y, int r);
void writecsMatrix(const char* fname,cs* A);
#endif
