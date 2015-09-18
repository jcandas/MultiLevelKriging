#ifndef NEIGHBOR3_H
#define NEIGHBOR3_H
#include <vector>
using namespace std;
class Neighbor3{
    void addTreeNode(int firstNode, int lastNode,int parentNode, double left, double right, double down, double up, double back, double front);
    void addTreeNode(vector<int>& q, int parentNode, int quadrant,  double left, double right, double down, double up, double back, double front);
    void quick(vector<double>& dis,vector<int>& y,int start,int end);
    double dist(double x1,double y1,double z1,double x2,double y2,double z2);

public:
    int** tree;
    //Each tree node has: 0:firstNode, 1:lastNode, 2-9:children, 
    // 10:parent, 11:level
    int treesize;
    double** bound; //left right bottom top back front
    vector<int> level;
    vector<int> lset;
    vector<int> levelstep;  //Mark the start of each level in the tree
    int* treenode; //mark the cube where the node belongs in the octree
		   //no longer used with HB scheme
    int** neighbors;
    int N;
    double** coord; //N x 3, (X,Y,Z) coords

    Neighbor3(double** coord1,int n, int van);
    void getQuadrant3(double left, double right, double down, double up, double back, double front, int firstNode, int lastNode,vector<int>& y);
    void levelset3(int parent,vector<int>& y);
//quicksort sorts in ascending order for dis and record the order in terms of indices to the original position in y
    void quicksort(vector<double>& dis,int size,vector<int>& y);
    void quicksort(vector<int>& dis,int size,vector<int>& y);
    void cover3(int index, int clevel, vector<int>& y);
    int contained(double* bound1, double* bound2);
};
#endif //NEIGHBOR3_H
