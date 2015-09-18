#include "hb.h"
#include "matrixio.h"
#include <stdio.h>

struct sortStr{
	CS_INT i;
	CS_ENTRY x;
};
int compare(const void* a, const void* b){
	return( *(int*)a - *(int*)b );
}
int sortcol_cmp(const void* a, const void* b){
	return( ((sortStr*)a)->i - ((sortStr*)b)->i );
}

HB::HB(int treesize, int n, int van){
	int MAX_LEVEL = 10;
	int MAX_LEVEL_SECOND = 25;
	int m = van*(van+1)*(van+2)/6;
	size = 0;
	N = n;

	blength = new int[treesize*2];
	support = new int*[treesize*2];
	coef = new double*[treesize*2];
	wavelet = new bool[treesize*2];
	cubenode = new int[treesize*2];
	treelink = new int*[treesize];
	for(int ii = 0; ii < treesize; ii++)
		treelink[ii] = new int[8*m]; //this can possibly be reduced
	linklength = new int[treesize];
	for(int ii = 0; ii < treesize; ii++)
		linklength[ii] = 0;
	moments = new double*[treesize*2];
	waveletSize = 0;
	waveletNodes = new int[n];
	waveletIndex = new int[treesize*2];
	hb2high = new int[treesize*2];
	wavlist = new int[n];
	wavlevels = new int*[n]; //using n for worst case
	for(int ii = 0; ii < n; ii++)
		wavlevels[ii] = new int[2];  //could allocate later 
	wavMaxLevel = 0;
	neighborHBLength = new int[treesize*2];
	neighborHBNodes = new int*[treesize*2];
	

	// Multi Level Neighborhoods

	MutliLevelneighborHBNodesLength = new int[MAX_LEVEL * treesize * 2]; 
	MutliLevelneighborHBNodes = new int*[MAX_LEVEL * treesize * 2];
	MutliLevelneighborHBNodesLengthTranspose = new int[MAX_LEVEL * MAX_LEVEL * treesize * 2]; 
	MutliLevelneighborHBNodesTranspose = new int*[MAX_LEVEL * MAX_LEVEL * treesize * 2];
	
	multilevelnodes = new int[MAX_LEVEL * treesize * 2];
	multilevelcount = 0;

	neighborHBCubesLength = new int[treesize*2];
	neighborHBCubes = new int*[treesize*2];

	neighborCubesLength = new int[treesize*2];
	neighborCubes = new int*[treesize*2];
	
	highTriple = new cs*[MAX_LEVEL_SECOND];


	high = new cs*[MAX_LEVEL_SECOND];
	for(int ii = 0; ii < MAX_LEVEL_SECOND; ii++){

		highTriple[ii] = NULL;

		high[ii] = NULL;
	}
	buffer = new int[n+10]; //10 additional units for debugging
	buffer[n] = 0;  /// track max i for availableHB 
	buffer[n+1] = 0;  ///track max j for availableHB

}

void HB::addHBNode(int cube,int* sup, double* coe, int blen,  bool wav, double* moment, int mlength)
	//Add a HB node to the list
{
	blength[size] = blen;
	support[size] = new int[blen];
	for(int ii = 0; ii < blen; ii++)
		support[size][ii] = sup[ii];
	coef[size] = new double[blen];
	for(int ii = 0; ii < blen; ii++)
		coef[size][ii] = coe[ii];
	wavelet[size] = wav;
	cubenode[size] = cube;
	treelink[cube][linklength[cube]] = size;
	linklength[cube]++;
	moments[size] = new double[mlength];
	for(int ii = 0; ii < mlength; ii++)
		moments[size][ii] = moment[ii];
	if(wav){
		waveletNodes[waveletSize] = size;
		waveletSize++;
	}
	size++;	
}

void HB::genMask(int waveletNode, double*& mask, int*& msupport, int& mlength){
//	printf("Entering genMask: %d\n",waveletNode);///
	if(waveletNode >= size) {
		fprintf(stderr,"Node not found in HB!\n");
		return;
	}
	int node = waveletNode;
	double*  mcoef = new double[size];
	for(int ii = 0; ii < size; ii ++) mcoef[ii] = 0;
	mlength = 0;
	for(int ii = 0; ii < blength[node]; ii++)
		if (support[node][ii] < 0) {
			mcoef[-support[node][ii]-1] = coef[node][ii];
		}else{
			double* subcoef;
			int * subsupp;
			int submlength = 0;
			genMask(support[node][ii],subcoef,subsupp,submlength);

			for(int jj = 0; jj < submlength; jj++){
				mcoef[subsupp[jj]] += coef[node][ii]*subcoef[jj];
			}
			delete[] subcoef;
			delete[] subsupp;
		}
	for(int ii = 0 ; ii < size; ii++) 
		if (mcoef[ii]!=0) mlength++;
	msupport = new int[mlength];
	mask = new double[mlength];
	int counter = 0;
	for(int ii = 0; ii < size; ii++){
		if (mcoef[ii]!=0){
			mask[counter] = mcoef[ii];
			msupport[counter] = ii;
			counter++;
		}
	}
	delete[] mcoef;
}

void HB::genAllMask(int waveletNode, double** mask, int** msupport, int* mlength){
//  Not a good method since it uses recursion.
//  Use FWTd3::wavMask or FWTd3::hbMask instead.
//  WaveletNode is the node index in HB
//	printf("Entering genAllMask: %d\n",waveletNode);///
	double* subcoef;
	int* subsupp;
	int submlength;
	int index;
	if(waveletNode >= size) {
		fprintf(stderr,"Node not found in HB!\n");
		return;
	}
	int node = waveletNode;
	double*  mcoef = new double[N];
	for(int ii = 0; ii < N; ii ++) mcoef[ii] = 0;
	for(int ii = 0; ii < blength[node]; ii++)
		if (support[node][ii] < 0) {
			mcoef[-support[node][ii]-1] = coef[node][ii];
		}else{
			genAllMask(support[node][ii],mask,msupport,mlength);
			index = waveletIndex[support[node][ii]];
			subcoef = mask[index];
			subsupp = msupport[index];
			submlength = mlength[index];
			for(int jj = 0; jj < submlength; jj++){
				mcoef[subsupp[jj]] += coef[node][ii]*subcoef[jj];
			}
		}
	index = waveletIndex[node];
	mlength[index] = 0;
	for(int ii = 0; ii < N; ii++)
                if (mcoef[ii]!=0) mlength[index]++;
	if(msupport[index]) delete[] msupport[index];
	if(mask[index]) delete[] mask[index];
	msupport[index] = new int[mlength[index]];
        mask[index] = new double[mlength[index]];
        int counter = 0;
        for(int ii = 0; ii < N; ii++){
                if (mcoef[ii]!=0){
                        mask[index][counter] = mcoef[ii];
                        msupport[index][counter] = ii;
                        counter++;
                }
        }
	delete[] mcoef;
}

//Merge supp2 into supp1, assuming that supp1 has memory available 
//for at least size N
void HB::mergeSupport(int* supp1, int& length1, int* supp2, int length2){
	int mlen = 0;

	//buffer pointers 
	int bi = 0; 
	int bj = 0; 

	//markers of supp1, supp2
	int i1 = 0; 
	int i2 = 0;
	
	// printf("mergeSupport: doing qsort... len1 = %d, len2 = %d\n",length1,length2);///
	// Not sorting the first list, assuming it is already sorted.

	qsort(supp2,length2,sizeof(int),compare);

	while(i2 < length2){
		if(bi == bj)// no buffered supp1
		{
			// insert supp2 and buffer supp1 if 

			while(i1 < length1 && supp1[i1] < supp2[i2]) i1++;
			if(i1 < length1 && supp1[i1] == supp2[i2]) 
				i2++;
			else{
				if(i1 < length1){
					buffer[bj] = supp1[i1];
					bj++;
				}
				supp1[i1] = supp2[i2];
				i1++;
				i2++;
			}
		}
		else //looking in the buffer first
		{
			while(bi < bj && buffer[bi] < supp2[i2]){
				if(i1 < length1){
					buffer[bj] = supp1[i1];
					bj++;
				}
				supp1[i1] = buffer[bi];
				bi++;
				i1++;
			}
			if(bi < bj){
			
				if (buffer[bi] == supp2[i2])
					i2++;
				else{
					if(i1 < length1){
						buffer[bj] = supp1[i1];
						bj++;
					}
					supp1[i1] = supp2[i2];
					i1++;
					i2++;
				}
			}	
		}
	}
	// buffer still non-empty
	if(bi < bj){
		if(i1 < length1){
			for(int i = i1; i < length1; i++){
				buffer[bj] = supp1[i];
				bj++;
			}
		}
		for(int i = bi; i < bj; i++){
			supp1[i1] = buffer[i];
			i1++;
		}
	}
	if(length1 < i1) length1 = i1;
	if(supp1[0] < 0) {
		printf("Negative support!");///
		exit(-1);
	}
}

void HB::trackAvailableHB(int si, int sj, int level){
	//storing the entries in the transposed triplet format
	// since the sortAvailableHB will do another transpose
//	if(buffer[N] < si) buffer[N] = si;///
//	if(buffer[N+1] < sj) buffer[N+1] = sj;///
//	cs_entry(availableHBTriple[level],si,sj,1); 
//	DavailableHB[si][sj] = 1; ///
}

// By sorting, we get efficient access to the computed HB entries
void HB::sortAvailableHB(int level){
	cs* T;
	cs* C;
	printf("Entering sortAvailableHB at level %d...\n",level);///
	// Sort availableHB first

    if(highTriple[level]){
	printf("sortAvailableHB: sort pointers... size = %d\n",highTriple[level]->nz);///

	// now sort high
	T = highTriple[level];
	if(high[level]) cs_spfree(high[level]);
	C = cs_compress(T);
	// sort high
	sortcols(C);
	
	// sort highTriple to have the same order as high
	decompress(C,highTriple[level]);
	high[level] = C;
    }else{
	highTriple[level] = cs_spalloc(N,N,N * N ,1,1);
	high[level] = cs_compress(highTriple[level]);
    }
    printf("sortAvailableHB: Done.\n");
}

//Always call sortAvailableHB first
int HB::testAvailableHB(int si, int sj,int level){
	si = hb2high[si];
	sj = hb2high[sj];
	cs* A = high[level];
	CS_INT ii = A->p[sj];
	CS_INT jj = A->p[sj+1];
	CS_INT* Ai = A->i;
	CS_INT kk = (ii + jj)/2;
	while(ii < kk){
		if (Ai[kk] < si) ii = kk;
		else if(Ai[kk] > si) jj = kk;
		else { 
			ii = kk; 
			break;
		}
		kk = (ii + jj)/2;
	}	
	if (Ai[ii] == si && A->p[sj] < A->p[sj+1]) return 1;
	else return 0;
}

double HB::highGet(int si, int sj, int level){
	cs* A = high[level];
	CS_INT ii = A->p[sj];
	CS_INT jj = A->p[sj+1];
	CS_INT* Ai = A->i;
	CS_INT kk = (ii + jj)/2;
	while(ii < kk){
		if (Ai[kk] < si) ii = kk;
		else if(Ai[kk] > si) jj = kk;
		else { 
			ii = kk; 
			break;
		}
		kk = (ii + jj)/2;
	}	
	if (Ai[ii] == si && A->p[sj] < A->p[sj+1]){
		return A->x[ii];
	}else {
		return 0;
	}
}

void HB::highSet(int si, int sj, int level, double val){

	cs* T = highTriple[level];
	cs* A = high[level];
	CS_INT ii = A->p[sj];
	CS_INT jj = A->p[sj+1];
	CS_INT* Ai = A->i;
	// use bisection to locate the element
	CS_INT kk = (ii + jj)/2;
	while(ii < kk){
		if (Ai[kk] < si) ii = kk;
		else if(Ai[kk] > si) jj = kk;
		else { 
			ii = kk; 
			break;
		}
		kk = (ii + jj)/2;
	}	
	// if element found, update value
	if (Ai[ii] == si && A->p[sj] < A->p[sj+1]){
		A->x[ii] = val; //may not be necessary
		T->x[ii] = val;
	}else{// insert a new element
		cs_entry(T,si,sj,val);
	}
}

void HB::sortcols(cs* A){
	sortStr* w = new sortStr[A->m];
	for(CS_INT j = 0; j < A->n; j++){
		CS_INT start = A->p[j];
		for(CS_INT ii = A->p[j]; ii < A->p[j+1]; ii++){
			w[ii - start].i = A->i[ii];
			w[ii - start].x = A->x[ii];
		}
	
		qsort(w, A->p[j+1] - A->p[j], sizeof(sortStr), sortcol_cmp);
	
		for(CS_INT ii = A->p[j]; ii < A->p[j+1]; ii++){
			A->i[ii] = w[ii - start].i;
			A->x[ii] = w[ii - start].x;
		}
	}
	delete[] w;
}

//Assuming that C and T has the same contents but in different order
void HB::decompress(cs* C, cs* T){
	for(CS_INT j = 0; j < C->n; j++)
		for(CS_INT ii = C->p[j]; ii < C->p[j+1]; ii++){
			T->i[ii] = C->i[ii];
			T->p[ii] = j;
			T->x[ii] = C->x[ii];
		}
}

double HB::cs_getEntry(cs* A,int si, int sj){
	CS_INT ii = A->p[sj];
	CS_INT jj = A->p[sj+1];
	CS_INT* Ai = A->i;
	CS_INT kk = (ii + jj)/2;
	while(ii < kk){
		if (Ai[kk] < si) ii = kk;
		else if(Ai[kk] > si) jj = kk;
		else { 
			ii = kk; 
			break;
		}
		kk = (ii + jj)/2;
	}	
	if (Ai[ii] == si && A->p[sj] < A->p[sj+1]){
		return A->x[ii];
	}else {
		return 0;
	}
}
