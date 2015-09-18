#ifndef FWTD3_H
#define FWTD3_H
#include <vector>
#include "neighbor3.h"
#include "cs.h"
#include "hb.h"
#include "diag_blk.h"
#define FWT_FORWARD 1
#define FWT_BACKWARD -1
#define FWT_ALL_LEVELS -1


class FWTd3{
    double calcLiftingError(double** vande,vector<int>& phi,vector<double>& phi1,double* sii, int m, int k, int n);
    double** moments_y;
    double* mome;
public:
    Neighbor3* nei;
    double** high;
    double** s;
    int* slength;
    double liftingError;
    int vanishingMoment;
    int pvanishingMoment;
    int minKwlevel;
    HB* hb;
    FWTd3(double** coord, int N, int van,int pvan, int d, int minlevel);
    void trafo(double** K1); //Wavelet transform on matrix
    void trafo(double* d,int direction, int level); //Wavelet transform on vector
    void trafo_diag_blk(Diag_Blk* db);
    void multilevel_trafo_diag_blk(Diag_Blk* db);
    void precond(int minlevel, double**& p, double*& T, int& m);
    double wavdistLift3(int i, int j);
    double wavdist3(int i, int j);
    double wavdist3Linf(int i, int j);
    void multiLift(int cube,int van);
    double coarseLift(int cube,int van,int pvan);
    void wavMask(double**& w);
    void wavMask(double**& mask,int**& msupp,int*& mlength);
    void hbMask(double**& mask,int**& msupp,int*& mlength);
    void getsupport(int si, int*& tsupp, int& tlength, int* tindex);
    void findNeighborHBNodes(int d);
    void findMultiLevelNeighborHBNodes(int d);
    void findMultiLevelNeighborHBNodesT(int d);
    void DebuggingfindMultiLevel(int d, int N);
    void StoreLinkDataMultiLevel(int level, int multilevel, int i, int MAX_LEVEL, int counter, int *levelcubes, int *nodes);
    void multileveldistcritUpper(int level, int current, int *nodes, int *nodeslength);
    void Kwblockmatrixcolmninit(int lj, int li, int treemark, int nnei, double **temp, int *FMMkk, int *FMMsj, double *eval, int *FMMtemp, int *t, int *nodes, int *nodeslength, int wavMaxLevel);
    void Kwblockmatrixcolmn(int lj, int li, int treemark, int nnei, double **temp, int *FMMkk, int *FMMsj, double *eval, int *FMMtemp, int *t, int *nodes, int *nodeslength, int wavMaxLevel);
    void Kwblockmatrixrow(int li, int lj, int currenttreemark, int nnei, double **temp, int *FMMkk, int *FMMsj, double *eval, int *FMMtemp, int *t, int *nodes, int *nodeslength, int wavMaxLevel);
    void Kwblockmatrixrowinit(int li, int lj, int currenttreemark, int nnei, double **temp, int *FMMkk, int *FMMsj, double *eval, int *FMMtemp, int *t, int *nodes, int *nodeslength, int wavMaxLevel);
    void KwStore(cs* csKw, int li, int lj, int *nodes, int *nodeslength,int wavMaxLevel);
    void KwStoreInt(cs* csKw, int li, int lj, int *nodes, int *nodeslength,int wavMaxLevel);
    void MultilevelNeighbor(int lj, int si, int *nodes, int *nodeslength, int *levelcubes, int *subdividedlevelcubes);
};

#endif //FWTD3_H
