//#include "fwtdata.h"
#include "neighbor3.h"
#include <math.h>
#include <iostream>
using namespace std;

Neighbor3::Neighbor3(double** coord1,int n, int van){
    cout <<"Entering Neighbor3" <<'\n';
    int ii,jj,kk;
    N = n;
    if(n == 0) return;
    cout <<"N = " <<N <<'\n';
    coord=coord1;
    double eps = 2.2204E-16;
    double left = coord[0][0];
    double right = coord[0][0];
    double down = coord[0][1];
    double up = coord[0][1];
    double back = coord[0][2];
    double front = coord[0][2];
    for(ii = 0; ii < n; ii++){
        if(left > coord[ii][0]) left = coord[ii][0];
        if(right < coord[ii][0]) right = coord[ii][0];
        if(down > coord[ii][1]) down = coord[ii][1];
        if(up < coord[ii][1]) up = coord[ii][1];
        if(back > coord[ii][2]) back = coord[ii][2];
        if(front < coord[ii][2]) front = coord[ii][2];
    }

    // Change to 2D, avoids starting with 0 width in the z direction.
    back = 0.1;
    front = 0.9;

    int logFactor = 2;
    tree = new int*[logFactor*n];
    bound = new double*[logFactor*n];
    treesize=0;

    cout <<"starting to build tree..[" <<logFactor*n <<"]\n";
    cout <<"  N = " <<n <<'\n';
    addTreeNode(0,n-1,0,left,right*(1+eps),down,up*(1+eps),back,front*(1+eps));
    for(ii = 0; ii < n; ii++)
	level.push_back(0) ;
    cout <<"level initialized" <<'\n';
    int i = 0;
    while(i < treesize){
        while(i < treesize && (tree[i][0] == tree[i][1])){
            i++;
	    //cout <<"   i = " <<i <<'\n';
        }
        if (i < treesize){
            double midx = (bound[i][0]+bound[i][1])/2;
            double midy = (bound[i][2]+bound[i][3])/2;
            double midz = (bound[i][4]+bound[i][5])/2;
            int midn = tree[i][0];
            while (coord[midn][0] < midx && midn<=tree[i][1]){
                midn++;
            }
            vector<int> q;
	    q.clear();
	    getQuadrant3(bound[i][0],midx,bound[i][2],midy,bound[i][4],midz,tree[i][0],midn-1,q);//ldb
            if (q.size()>0){
                addTreeNode(q,i,8,bound[i][0],midx,bound[i][2],midy,bound[i][4],midz);
            }
	    q.clear();
            getQuadrant3(bound[i][0],midx,midy,bound[i][3],bound[i][4],midz,tree[i][0],midn-1,q);//lub
            if (q.size()>0){
                addTreeNode(q,i,6,bound[i][0],midx,midy,bound[i][3],bound[i][4],midz);
            }
	    q.clear();
            getQuadrant3(bound[i][0],midx,midy,bound[i][3],midz,bound[i][5],tree[i][0],midn-1,q);//luf
            if (q.size()>0){
                addTreeNode(q,i,2,bound[i][0],midx,midy,bound[i][3],midz,bound[i][5]);
            }
	    q.clear();
            getQuadrant3(bound[i][0],midx,bound[i][2],midy,midz,bound[i][5],tree[i][0],midn-1,q);//ldf
            if (q.size()>0) {
                addTreeNode(q,i,4,bound[i][0],midx,bound[i][2],midy,midz,bound[i][5]);
            }
	    q.clear();
            getQuadrant3(midx,bound[i][1],bound[i][2],midy,bound[i][4],midz,midn,tree[i][1],q);//rdb
            if (q.size()>0){
                addTreeNode(q,i,9,midx,bound[i][1],bound[i][2],midy,bound[i][4],midz);
            }
	    q.clear();
            getQuadrant3(midx,bound[i][1],midy,bound[i][3],bound[i][4],midz,midn,tree[i][1],q);//rub
            if (q.size()>0){
                addTreeNode(q,i,7,midx,bound[i][1],midy,bound[i][3],bound[i][4],midz);
            }
	    q.clear();
            getQuadrant3(midx,bound[i][1],bound[i][2],midy,midz,bound[i][5],midn,tree[i][1],q);//rdf
            if (q.size()>0){
                addTreeNode(q,i,5,midx,bound[i][1],bound[i][2],midy,midz,bound[i][5]);
            }
	    q.clear();
            getQuadrant3(midx,bound[i][1],midy,bound[i][3],midz,bound[i][5],midn,tree[i][1],q);//ruf
            if(q.size()>0){
                addTreeNode(q,i,3,midx,bound[i][1],midy,bound[i][3],midz,bound[i][5]);
            }
            i++;
        }
    }
	cout <<"Final treesize: " <<treesize <<'(' <<tree[treesize-1][11] <<')' <<'\n';
	cout << "Done" <<'\n';   
	cout <<"Here 1" <<'\n';

    treenode = new int[n];
    for(ii = 0;  ii < n;  ii++){
        int node=0;
        while(tree[node][0] != ii) node++;
        level[ii] = tree[node][11];
        treenode[ii] = node;
    }
    



    // generate lset
    vector<int> lv = level;
    for (int ii = 0; ii < N; ii++){
        lv[ii] = lv[ii]*N + ii;
    }
    quicksort(lv,N,lset);
    lv.clear();
    
    // generate levelstep
    int minli=-1;
    for (ii = 0; ii < treesize; ii++)
        if (minli < tree[ii][11]){
            levelstep.push_back(ii);
            minli = tree[ii][11];
        }
    levelstep.push_back(treesize);

    cout <<"treenode generated" <<'\n';
}

void Neighbor3::addTreeNode(int firstNode, int lastNode,int parentNode, double left, double right, double down, double up, double back, double front){
    	tree[treesize] = new int[12];
        tree[treesize][0] = firstNode;
        tree[treesize][1] = lastNode;
        tree[treesize][2] = 0;
        tree[treesize][3] = 0;
        tree[treesize][4] = 0;
        tree[treesize][5] = 0;
        tree[treesize][6] = 0;
        tree[treesize][7] = 0;
        tree[treesize][8] = 0;
        tree[treesize][9] = 0;
        tree[treesize][10] = parentNode;
        if(treesize>0)
            tree[treesize][11] = tree[parentNode][11]+1;
        else
            tree[treesize][11] = 0;
		
        bound[treesize] = new double[6];
        bound[treesize][0]=left;
        bound[treesize][1]=right;
        bound[treesize][2]=down;
        bound[treesize][3]=up;
        bound[treesize][4]=back;
        bound[treesize][5]=front;
	
        treesize++;
}

void Neighbor3::addTreeNode(vector<int>& q,int parentNode, int quadrant, double left, double right, double down, double up, double back, double front){
	// cout <<"Tring to add node.." <<"\n  ";
	tree[treesize] = new int[12];
        tree[treesize][0] = q[0]; //first node
        tree[treesize][1] = q[q.size()-1]; //last node
        tree[treesize][2] = 0;
        tree[treesize][3] = 0;
        tree[treesize][4] = 0;
        tree[treesize][5] = 0;
        tree[treesize][6] = 0;
        tree[treesize][7] = 0;
        tree[treesize][8] = 0;
        tree[treesize][9] = 0;
        tree[treesize][10] = parentNode;
        if(treesize>0)
            tree[treesize][11] = tree[parentNode][11]+1;
        else
	  tree[treesize][11] = 0;
        tree[parentNode][quadrant] = treesize;

	
        bound[treesize] = new double[6];
        bound[treesize][0]=left;
        bound[treesize][1]=right;
        bound[treesize][2]=down;
        bound[treesize][3]=up;
        bound[treesize][4]=back;
        bound[treesize][5]=front;

        treesize++;
}

void Neighbor3::getQuadrant3(double left, double right, double down, double up, double back, double front, int firstNode, int lastNode, vector<int>& y){
        int node=firstNode;
	y.clear();
        while (node <= lastNode && coord[node][0] < left)
            node++;
        while (node <= lastNode && coord[node][0] < right){
            if (coord[node][1] >= down && coord[node][1] <up &&
            coord[node][2] >=back && coord[node][2]<front)
                y.push_back(node);
            node++;
        }
}

void Neighbor3::levelset3(int parent,vector<int>& y){
        int ii;
	y.clear();
        if(parent <= 0 || parent > treesize-1)
            return;
        int level = tree[parent][11];
        ii=1;
        while (ii < treesize && tree[ii][11] < level)
            ii++;
        int y1 = ii;
        while (ii < treesize && tree[ii][11] == level)
            ii++;
        int y2 = ii;
        for (ii = y1; ii < y2; ii++) y.push_back(ii);
}

void Neighbor3::quicksort(vector<double>& dis,int size,vector<int>& y){
        vector<double> dis1;
	int ii;
	y.clear();
        for (ii=0; ii < size; ii++){
            y.push_back(ii);
            dis1.push_back(dis[ii]);
        }
        quick(dis1,y,0,size-1);
	dis1.clear();
}

void Neighbor3::quicksort(vector<int>& dis,int size,vector<int>& y){
        vector<double> dis1;
	y.clear();
	int ii;
        for (ii=0; ii < size; ii++){
            y.push_back(ii);
            dis1.push_back(dis[ii]);
        }
        quick(dis1,y,0,size-1);
	dis1.clear();
}

void Neighbor3::quick(vector<double>& dis,vector<int>& y,int start,int end){
        int i = start;
        int j = end;
        double t = dis [j];
        int yt = y[j];
        while(i < j){
            while(i < j && dis[i] <= t) i++;
            if(i < j) {
                dis[j] = dis[i];
                y[j] = y[i];
                j--;
            }
            while(i < j && t <= dis[j]) j--;
            if (i < j){
                dis[i] = dis[j];
                y[i] = y[j];
                i++;
            }
        }
        dis[j] = t;
        y[j] = yt;
        if (start < i-1) quick(dis,y,start,i-1);
        if (i+1 < end) quick(dis,y,i+1,end);
    }

double Neighbor3::dist(double x1,double y1,double z1,double x2,double y2,double z2){
	return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
}

//test whether bound1 is contained in bound2
int Neighbor3::contained(double* bound1, double* bound2){
  if(bound1[0] >= bound2[0] && bound1[1] <= bound2[1] &&
	bound1[2] >= bound2[2] && bound1[3] <= bound2[3] &&
	bound1[4] >= bound2[4] && bound1[5] <= bound2[5])
	return 1;
  else
	return 0;
}

void Neighbor3::cover3(int index, int clevel, vector<int>& y){
	int node=0;
	y.clear();
	while (node < treesize && 
		(tree[node][0] != index || tree[node][11] != clevel))
		node++;
	if (node < treesize)
		getQuadrant3(bound[node][0],bound[node][1],bound[node][2],bound[node][3],bound[node][4],bound[node][5],tree[node][0],tree[node][1],y);
}
   
