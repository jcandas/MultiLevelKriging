#ifndef DIAG_BLK_H
#define DIAG_BLK_H
#include "cs.h"
class Diag_Blk{
public:
	int max_level;
	int sparse_level;  //assume sparse rep for level >= sparse_level
	cs** sparse_blk;
	int* blk_size;
	double** dense_blk;
	double** dense_lu;
	int** ipiv;
};

#endif //DIAG_BLK_H
