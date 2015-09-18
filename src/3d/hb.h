#ifndef HB_H
#define HB_H

#include "diag_blk.h"

#include "cs.h"

using namespace std;

class HB{
public:
	int size;
	//total number of nodes in HB
	
	int N;
	//total number of points refered to by coord	
	
	int arbolsize;
	//Treesize arbol is tree in spanish

	int * blength;
	//length of support (in terms of HB nodes) for each HB node
	
	int **support;
	//positive index for HB element; negative index for actual element
	//The negative value of the index is the negative of actual index to the point
	//refered to by coord -1
	
	double **coef; 
	//coefficients corresponding to support
	
	bool *wavelet; 
	//flag for scaling function or wavelet
	
	int *cubenode;
	//index to the cube in the tree (usually the cube on the highest level if linked by 
	//multiple cubes.

	int **treelink;
	//link from tree nodes to the scaling functions/wavelets in the box
	//It is possible for multiple tree nodes to link to the same scaling function, since
	//lifting only happens when the number of scaling function on one node exceeds # of moments

	int *linklength;
	//length of the treelink for each node
	
	double **moments;
	//moments, only for scaling functions

//	int *momementlink;
	//link from a scaling func HB node to it's moment
	
	int waveletSize;
	//size of the list for wavelets

	int primary_waveletSize;
	//size of the wavelets after 1st MultiLift

	int* waveletNodes;
	//A list for the wavelet nodes (scaling functions do not count)
	//except the top nnei number of scaling functions appended at the end
	

	int* waveletIndex;
	//Pointers from HB node to its position in the level-ordered list of wavelets
	
	int* hb2high;
	//Pointers from HB node to its natural position refered to by coord
	//this pointer value is unique within each level, but may overlap 
	//acrosslevels (Pointers for 2 HB nodes with seperate support is 
	//guaranteed to be different)

	int* wavlist;
	//Pointers from level-ordered wavelets to index in waveletNodes

	int** wavlevels;
	//pointers to the starting and ending index of each wavelet level in the level-ordered list

	//finest level for lifted wavelets
	//Could also be the coarse-lifted wavelets, if they exist.
	int wavMaxLevel;

	int** neighborHBNodes;
	int* neighborHBLength;

	//Multi Level Cubes
	int* multilevelnodes;
	int  multilevelcount;

	int** neighborHBCubes;
	int*  neighborHBCubesLength;

	int** neighborCubes;
	int*  neighborCubesLength;

	int** MutliLevelneighborHBNodes;
	int*  MutliLevelneighborHBNodesLength;

	int** MutliLevelneighborHBNodesTranspose;
	int*  MutliLevelneighborHBNodesLengthTranspose;
	
	//HB basis
	double** mask;
	int** msupport;
	int* mlength;

	Diag_Blk* db;

	cs** availableHBTriple;
	cs** availableHB;
	cs** high;//high in compressed format for all level
	cs** highTriple;//high in triplet format for all level
//	cs* highTripleIndex; //pointer from sparse high to sparse highTriple
	//the highTriple for coarser level are always copied from finer level

	///for debug only
//	int** DavailableHB;
//	double** Dhigh;
	///

	//buffer for sorting/merging functions (size N)
	int* buffer;

	char* options_file;
	int kernelType;
	double kernelParameter;

	HB(int treesize, int n, int van);
	
	void addHBNode(int cube,int* sup, double* coe, int blen,  bool wav, double* moment, int mlength);
	//Add a HB node to the list
	
	void genMask(int waveletNode, double*& mask, int*& msupport, int& mlength);
	void genAllMask(int waveletNode, double** mask, int** msupport, int* mlength);
	void mergeSupport(int* supp1, int& length1, int* supp2, int length2);
    	void trackAvailableHB(int si,int sj,int level);
    	void sortAvailableHB(int level);
    	int testAvailableHB(int si, int sj,int level);
	double highGet(int si, int sj,int level);
	void highSet(int si, int sj,int level, double val);
	double cs_getEntry(cs* A, int si, int sj);
private:
	void sortcols(cs* A);
	void decompress(cs* C, cs* T);
};
#endif //HB_H
