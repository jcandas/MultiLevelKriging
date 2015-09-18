#include "dgesvd.h"
#include "fwtd3.h"
#include "neighbor3.h"
#include "vande3.h"
#include "matrixio.h"
#include "hb.h"
#include "linearSolver.h"
#include <math.h>
// extern "C" double dsecnd();

#define MAX_DIRECT_KERNEL 2*7200000


	vector<double> Vande3::cx, Vande3::cy, Vande3::cz;
	vector<double> Vande3::polyx, Vande3::polyy, Vande3::polyz;
	double Vande3::temp_a, Vande3::temp_b;

FWTd3::FWTd3(double** coord, int N, int van, int pvan, int d, int minlevel){
	cout << "Entering fwtd3 " << '\n';
        int ii,jj,kk;
	vanishingMoment = van;
	pvanishingMoment = pvan;

	printf("Here \n");
	minKwlevel = minlevel;

	printf("Min Level %d \n",minKwlevel);

	int nnei = van*(van+1)/2; //2D
	int pnnei = pvan*(pvan+1)/2; //2D


        vector<int> interp;
	cout << "FWTd3: nnei = " <<nnei << '\n';

	cout << "FWTd3: generating neighborhoods ..." <<'\n';
	
	moments_y = new double*[nnei];
	mome = new double[nnei];
	for(ii = 0; ii < nnei; ii++) moments_y[ii] = new double[1];

        nei = new Neighbor3(coord, N, nnei);
	cout << "FWTd3:  Done " << '\n';

        liftingError = 0;

	hb = new HB(nei->treesize,N,van);

	cout << "Start lifting.." <<N <<'\n';

	printf("nei->tree[%d][0] = %d\n",1431,nei->tree[1431][0]);///	
	printf("nei->tree[%d][1] = %d\n",1431,nei->tree[1431][1]);///	
	
	multiLift(0,van);
	
	hb->primary_waveletSize = hb->waveletSize;


	//append the unlifted nodes to the end of the waveletNodes
	for(ii = 0; ii < nnei && ii < hb->linklength[0]; ii++)
		hb->waveletNodes[hb->primary_waveletSize+ii] = hb->treelink[0][ii];
	
	cout << "FWTd3: lifting done (maxError = " <<liftingError 
		<<", HB size = " <<hb->size <<")" << '\n';

	vector<int> wavcube,wavlist;
	for(ii = 0; ii < hb->primary_waveletSize; ii ++) 
		wavcube.push_back(-hb->cubenode[hb->waveletNodes[ii]]);
	nei->quicksort(wavcube,hb->primary_waveletSize,wavlist);
	for(ii = 0; ii < nnei; ii++)
		wavlist.push_back(hb->primary_waveletSize+ii);
	printf("Wavelet Size: %d,%d,%d\n",hb->primary_waveletSize,
		hb->waveletSize,wavlist.size());

	for(ii = 0; ii < N; ii++)
		hb->wavlist[ii] = wavlist[ii];

	// initialize wavlevels
	int** wavlevels = hb->wavlevels;
	int minl = -1;
	int treemark;
	for(ii = hb->primary_waveletSize; ii > 0; ii--){
		treemark = hb->cubenode[hb->waveletNodes[wavlist[ii-1]]];
		// printf("(%d:%d) ",treemark,nei->tree[treemark][11]);///
		if(minl < nei->tree[treemark][11]){
			minl = nei->tree[treemark][11];
			printf("\nminl = %d; ii = %d\n",minl,ii);///
			wavlevels[minl][1] = ii;
			if(minl > 0)
				wavlevels[minl-1][0] = ii;
		}
	}
	for(ii = minl; ii < nei->levelstep.size(); ii++)
		wavlevels[ii][0] = 0; //The finest level starts first
	for(ii = minl+1; ii < nei->levelstep.size(); ii++)
		wavlevels[ii][1] = 0; //The finest level starts first, some unlifted levels are not visible after lifting
	hb->wavMaxLevel = minl;
	printf("hb->wavMaxLevel = %d\n",hb->wavMaxLevel);
		
		
	int* hb2high = hb->hb2high;
	int* waveletIndex = hb->waveletIndex;
	for(ii = 0; ii < hb->size; ii++) hb2high[ii] = -1;
	treemark = nei->treesize - 1;
	int maxLevel = nei->levelstep.size() - 1;
	minl = 1;
	for(int li = maxLevel-1; li >= minl-1; li--){
	  // cout <<"treemark at the end of level("<<li <<"): "<<treemark <<endl;
		while(treemark >=0 && nei->tree[treemark][11] == li){
			int si = hb->treelink[treemark][0];
			for(jj = 0; jj < hb->blength[si]; jj++)
				if(hb->support[si][jj] < 0)
					hb2high[hb->treelink[treemark][jj]] = -hb->support[si][jj] -1;
				else
					hb2high[hb->treelink[treemark][jj]] = hb2high[hb->support[si][jj]];
			treemark--;
		}
	}
	int* wavlist_1 = new int[N];
	for(ii = 0; ii < N; ii++)
		wavlist_1[hb->hb2high[hb->waveletNodes[wavlist[ii]]]] = ii;
	for(ii = 0; ii < hb->size; ii++)
		waveletIndex[ii] = wavlist_1[hb2high[ii]];//breakdown?
	delete[] wavlist_1;
	wavlist.clear();
	wavcube.clear();

	coarseLift(0,van,pvan);

	double** wm;
	int** ws;
	int* wl;

	hbMask(hb->mask,hb->msupport,hb->mlength);

	// Find Neighbor nodes. Set to d = 0 if only diagonal preconditioner used
	// Other set it to larger for sparse output of K_W matrix
	

	if (d == 0) findNeighborHBNodes(0);
	
	else{
	// Multi Level Sparse Representation distance criterion
	if (d > 0) {
	findMultiLevelNeighborHBNodes(d);
	findMultiLevelNeighborHBNodesT(d);
	}
	}
	
	if (d == -1){
	  findMultiLevelNeighborHBNodes(0);
	  findMultiLevelNeighborHBNodesT(0);
	}
}

void FWTd3::trafo(double** K1){
        int ii,jj,kk,ll;
	int N = nei->N;
	int van = vanishingMoment;
        // int nnei = van*(van+1)*(van+2)/6;
	int nnei = van*(van+1)/2; //2D
	
        vector<int> levelstep = nei->levelstep;
        cout << "levelstep initiated: \n";


	// determin the minimal lifted level
	int minl=1;
	// rearrange the un-lifted wavelets into one single scale level

	double** temp = new double*[nnei*8];
	for(ii = 0; ii < nnei*8; ii++)
		temp[ii] = new double[N];
	int* hb2high = hb->hb2high;

	cout << "Starting wavelet transform.." <<endl;	
        high = K1;
        int iter;
        for (iter = 0; iter < 2; iter++){
            int treemark = nei->treesize-1;
	    // cout <<"  treemark initialized : "<< treemark <<endl;
            int maxLevel = levelstep.size()-1;
            int li;
            for (li = maxLevel-1; li >= minl-1; li--){
	    	cout <<"treemark at the end of level("<<li <<"): "<<treemark <<endl;
                while (treemark >= 0 && nei->tree[treemark][11] == li){
			int si = hb->treelink[treemark][0];
			// size of support is the same as the total number of scaling func and wavelets
			// and the size of support is the same for all hb nodes in the cube
			for(jj = 0; jj < hb->blength[si]; jj++)

			if(hb->linklength[treemark]>nnei)
				for(kk = 0; kk < N; kk++)
					temp[jj][kk] = high[hb2high[hb->treelink[treemark][jj]]][kk];

			// compute scaling functions and wavelets
			// printf("treemark: %d(%d)\n",treemark,hb->linklength[treemark]);///
			if(hb->linklength[treemark]>nnei)
			for(jj = 0; jj < hb->linklength[treemark]; jj++){
				si = hb->treelink[treemark][jj];
				// printf("  Accessing hb->treelink[%d][%d] : %d -> %d\n",treemark,jj,si,hb2high[si]);
				// Assuming all subnodes have the exactly same set of support
				for(ll = 0; ll < N; ll++){
					high[hb2high[si]][ll] = 0;
					for(kk = 0; kk < hb->blength[si]; kk++)
						high[hb2high[si]][ll] += hb->coef[si][kk]*temp[kk][ll];
				}
			}
			treemark--;
	        }
            }
	    
            double t;
            for (ii = 0; ii < N-1; ii++)
                for (jj = ii+1; jj < N; jj++){
                    t = high[ii][jj];
                    high[ii][jj] = high[jj][ii];
                    high[jj][ii] = t;
                }
        }
	
	for(ii = 0; ii < nnei*8; ii++) delete[] temp[ii];
	delete[] temp;

	//? This step requres O(n^2) memory
	int* wavlist = hb->wavlist;
	high = new double*[N];
	double* t = new double[N];
	for(ii = 0; ii < N; ii++){
		for( jj = 0; jj < N; jj++)
			t[jj] = K1[hb2high[hb->waveletNodes[wavlist[ii]]]][hb2high[hb->waveletNodes[wavlist[jj]]]];
		high[ii] = t;
		t = K1[hb2high[hb->waveletNodes[wavlist[ii]]]];
	}
	delete[] t;
	delete[] K1;

	cout << "FWTd3: Wavelet transform done " << '\n';
}

void FWTd3::trafo_diag_blk(Diag_Blk* db){
        int ii,jj,kk,ll;
	int N = nei->N;
	int van = vanishingMoment;

	
	// Initialize vectors for Matern kernel summation
	int nsite = MAX_DIRECT_KERNEL;
	double* dist   = new double[nsite];
	double* site   = new double[nsite];                  // Array of interpolation sites
	double* rout   = new double[nsite];          // spline evaluation results
	
	int nnei = van*(van+1)/2; //2D

	int* t = new int[N]; // temp memory for HB support 
	int tlen;  // length of it
	double val;  // temporary memory for a double
	double t0, t1;	
	bool coarseLifted = false;
        printf("HERE vanishingMoment = %d, pvanishingMoment = %d\n",vanishingMoment, pvanishingMoment);///
	if(vanishingMoment > pvanishingMoment) coarseLifted = true;

        vector<int> levelstep = nei->levelstep;
        cout << "Trafo_diag_blk: levelstep initiated: \n";

	//determin the minimal lifted level
	int minl=1;
	
	// Make note here. 
	double** temp = new double*[nnei*8];
	for(ii = 0; ii < nnei*8; ii++)
		temp[ii] = new double[N];
	int* hb2high = hb->hb2high;

	cout << "Trafo_diag_blk: Starting wavelet transform.." <<endl;	
        int treemark = nei->treesize-1;
	
	// May have a problem if coarse-lifting exitsts
	int wavMaxLevel = hb->wavMaxLevel;

        printf("  HERE wavMaxLevel = %d\n",wavMaxLevel);///
	if(coarseLifted) wavMaxLevel--;
        printf("  HERE wavMaxLevel = %d\n",wavMaxLevel);///

	while(treemark >= 0 && nei->tree[treemark][11] > wavMaxLevel)
		treemark--;

	// cout <<"  treemark initialized : "<< treemark <<endl;
        int maxLevel = levelstep.size()-1;
        int li, li_mark;
	int FMMcount;
	int* FMMkk = new int[N]; // index in the neighborHBNodes list
	int* FMMsj = new int[N]; // index in the HB list
	double* eval = new double[N];  // evaluation of transformed entries
	int* FMMtemp = new int[N*2]; // temparory memory for FMMevalHB

	int* wavlist = hb->wavlist;
	double* blk;
	double* blk_lu;
	int* ipiv;
	cs* csblk;
	cs* cs_t;

	hb->sortAvailableHB(wavMaxLevel+1);

        printf("HERE wavMaxLevel = %d\n",wavMaxLevel);///
	for (li = wavMaxLevel; li >= minl-1; li--){
	    // record the position for transform in the other direction
		t0 = dsecnd();
		li_mark = treemark; 
		int iter;
	    // Need to apply the HB transform twice for 2 directions
	    	cout <<"trafo_diag_blk: treemark at the end of level("<<li <<"): "<<treemark <<endl;
		//sort only at the beginning of each level
		hb->sortAvailableHB(li);
		//loop through the cubes
                while (treemark >= 0 && nei->tree[treemark][11] == li){
		    int si = hb->treelink[treemark][0];
			// size of support is the same as the total number of scaling func and wavelets
			// and the size of support is the same for all hb nodes in the cube
			// only collecting neighboring finer level scaling functions
		        // printf("preparing temp...\n");///
			//Assuming that all HB nodes under the same cube has the same neighbors
		    if(hb->linklength[treemark]>nnei){
			// merge support for all neighbors of si
			tlen = 0;
			for(kk = 0; kk < hb->neighborHBLength[si]; kk++){
				int sk = hb->neighborHBNodes[si][kk];	      
			
				if (hb->support[sk][0] < 0)
					// sk is in the finest cube,
					// which should contain only 1 point
					hb->mergeSupport(t,tlen,&sk,1);
				else	
					hb->mergeSupport(t,tlen,hb->support[sk],hb->blength[sk]);
			}
			// loop through the support of si
			for(jj = 0; jj < hb->blength[si]; jj++){
				int sj = hb->support[si][jj];
				FMMcount = 0;
				// printf("Colmns: si %d hb->neighborHBLength[%d] = %d\n",si,sj,hb->neighborHBLength[sj]);///
				// kk should go through the support of all neighbors of sj
				for(kk = 0; kk < tlen; kk++){
					// use FMM for missing entries in the finer level
					int sk = t[kk];
					// printf("kk %d sk %d sj %d \n",kk,sk,sj); 
					if(hb->testAvailableHB(sj,sk,li+1))
					    temp[jj][kk] = hb->highGet(hb2high[sj],hb2high[sk],li+1);
					else{
					    //record the position to use FMMevalHB
					    FMMkk[FMMcount] = kk;
					    FMMsj[FMMcount] = sk;
					    FMMcount++;
					    // no highSet here, but will be set later after FMMevalHB
					    hb->trackAvailableHB(sj,sk,li);
					}
				}
				// Evaluate the interaction between sj and neighbors
				// printf("Ready for FMMevalHB\n");///
				if(FMMcount>0)
					// FMMevalHB should handle negative HB index
					// as pointer to coordinates
					// Here the indices are tranposed, but should still be correct assuming A is symmtric
				  FMMevalHB(hb,FMMsj,FMMcount,sj,eval,FMMtemp,nei,dist,rout,site);
				// printf("Returned from FMMevalHB\n");///
				for(kk = 0; kk < FMMcount; kk++){
					temp[jj][FMMkk[kk]] = eval[kk];
				}
			}
			// compute scaling functions and wavelets in the cube
			// printf("treemark: %d(%d)\n",treemark,hb->linklength[treemark]);///
			for(jj = 0; jj < hb->linklength[treemark]; jj++){
				int sj = hb->treelink[treemark][jj];
				// Assuming all subnodes have the exactly same set of support
				// Adding up for only the neighboring HB Nodes
				for(ll = 0; ll < tlen; ll++){
					val = 0;
					int sl = t[ll];
					for(kk = 0; kk < hb->blength[sj]; kk++)
						val += hb->coef[sj][kk]*temp[kk][ll];
					// printf("Cols: highset ...\n");///
					hb->highSet(hb2high[sj],hb2high[sl],li,val);					
					// printf("Cols: trackAvailableHB ...\n");///
				}
			}
		    }
		    treemark--;
		}
		treemark = li_mark;
		// now operating on rows
	    	cout <<"trafo_diag_blk: Rows: treemark at the end of level("<<li <<"): "<<treemark <<endl;
		hb->sortAvailableHB(li);

                while (treemark >= 0 && nei->tree[treemark][11] == li){
		    int si = hb->treelink[treemark][0];
			// size of support is the same as the total number of scaling func and wavelets
			// and the size of support is the same for all hb nodes in the cube
			// only collecting neighboring finer level scaling functions
		        // printf("Rows: preparing temp...\n");///
		    if(hb->linklength[treemark]>nnei){
			tlen = 0;
			for(kk = 0; kk < hb->neighborHBLength[si]; kk++){
				int sk = hb->neighborHBNodes[si][kk];
				if (hb->support[sk][0] < 0)
					// sk is in the finest cube,
					// which should contain only 1 point
					hb->mergeSupport(t,tlen,&sk,1);
				else	
					hb->mergeSupport(t,tlen,hb->support[sk],hb->blength[sk]);
				// Alternatively, we could get the support from neighboring cubes
			}
			if(hb->neighborHBLength[si] != tlen){
				exit(-1);
			}
			// loop through the support of si
			// collect entries for the HB transform
			for(jj = 0; jj < hb->blength[si]; jj++){
				// sj loops through the support of si
				int sj = hb->support[si][jj];
				// printf("Rows: si %d hb->neighborHBLength[%d] = %d\n",si,sj,hb->neighborHBLength[sj]);///
				// loop through the support of neighbors of si
				for(kk = 0; kk < tlen; kk++){
					// all entries should be available now
				        int sk = t[kk];
					//printf("kk %d sk %d sj %d (%d,%d) \n",kk,sk,sj,hb2high[sk],hb2high[sj]);
					if(!hb->testAvailableHB(sk,sj,li)){
					  printf("Something wrong with testAvailableHB: (%d, %d), Treemark %d kk %d tlen %d. \
                                                  Increase ptilde or decrease tau \n",
						  sk,sj,treemark, kk, tlen);
						exit(-1);
					}
					
					temp[jj][kk] = hb->highGet(hb2high[sk],hb2high[sj],li);
				}
			}
			// compute scaling functions and wavelets
			// printf("Rows: treemark: %d(%d)\n",treemark,hb->linklength[treemark]);///
			for(jj = 0; jj < hb->linklength[treemark]; jj++){
				// sj loops through the lifted HB nodes in the cube
				int sj = hb->treelink[treemark][jj];
				// Assuming all subnodes have the exactly same set of neighbors
				// Adding up for only the neighboring HB Nodes
				for(ll = 0; ll < tlen; ll++){
					val = 0;
					int sl = t[ll];
					for(kk = 0; kk < hb->blength[sj]; kk++)
						val += hb->coef[sj][kk]*temp[kk][ll];
					hb->highSet(hb2high[sl],hb2high[sj],li,val);
				}
				// The neighbors of sj is the same as those of si
				// The points in t and neighborHBNodes are for the same
				// locations in high but on differnt level and in different orders
				for(ll = 0; ll < hb->neighborHBLength[si]; ll++){
					int sl = hb->neighborHBNodes[si][ll];
					hb->trackAvailableHB(sl,sj,li);
				}
			}
		    }
		    treemark--;
		}
		t1 = dsecnd();
		printf("Pyramid scheme used %lf seconds for level %d.\n",t1 - t0,li);

		cs_spfree(hb->highTriple[li+1]);
		cs_spfree(hb->high[li+1]);

		hb->sortAvailableHB(li);

		db->blk_size[li] = hb->wavlevels[li][1] - hb->wavlevels[li][0]; 
		if(li < db->sparse_level && db->blk_size[li] < 1e4){
			printf("trafo_diag_blk: LU for dense blk... %d\n",db->blk_size[li]);///	
			int blk_n = db->blk_size[li];
			blk = new double[blk_n * blk_n];
			blk_lu = new double[blk_n * blk_n];
			ipiv = new int[blk_n];
			for(ii = 0 ; ii < blk_n; ii++){
				for(jj = 0; jj < blk_n; jj++)
					blk[ii + blk_n*jj] = hb->highGet(hb2high[hb->waveletNodes[hb->wavlist[ii+hb->wavlevels[li][0]]]],hb2high[hb->waveletNodes[hb->wavlist[jj+hb->wavlevels[li][0]]]],li);
			}
			for(ii = 0; ii < blk_n*blk_n; ii++)
				blk_lu[ii] = blk[ii];
			lu_init(blk_n,blk_n,blk_lu,ipiv);
			db->dense_blk[li] = blk;
			db->dense_lu[li] = blk_lu;
			db->ipiv[li] = ipiv;
		}else{
			printf("trafo_diag_blk: cs for sparse blk...\n");///	
			csblk = cs_spalloc(db->blk_size[li],db->blk_size[li],db->blk_size[li]*8,1,1);
			for(ii = 0; ii < db->blk_size[li]; ii++){
				int si = hb->waveletNodes[hb->wavlist[ii + hb->wavlevels[li][0]]];
				for(jj = 0; jj < hb->neighborHBLength[si]; jj++){
					int sj = hb->neighborHBNodes[si][jj];
					if(hb->wavelet[sj])
					  cs_entry(csblk, ii, hb->waveletIndex[sj] - hb->wavlevels[li][0], hb->highGet(hb2high[si],hb2high[sj],li)); 
					//(sorting not necessary) storing as transposed for the sorting
				}
			}
			db->sparse_blk[li - db->sparse_level] = cs_compress(csblk);
			cs_spfree(csblk);
		}
        }
	if(coarseLifted){
		li = wavMaxLevel + 1;
		db->blk_size[li] = hb->wavlevels[li][1] - hb->wavlevels[li][0]; 
		printf("trafo_diag_blk: LU for dense blk...%d\n",db->blk_size[li]);///	
		int blk_n = db->blk_size[li];
		blk = new double[blk_n * blk_n];
		blk_lu = new double[blk_n * blk_n];
		ipiv = new int[blk_n];
		for(ii = 0 ; ii < blk_n; ii++){
			for(jj = 0; jj < blk_n; jj++)
				blk[ii + blk_n*jj] = hb->highGet(hb2high[hb->waveletNodes[hb->wavlist[ii+hb->wavlevels[li][0]]]],hb2high[hb->waveletNodes[hb->wavlist[jj+hb->wavlevels[li][0]]]],minl-1);
		}
		for(ii = 0; ii < blk_n*blk_n; ii++)
			blk_lu[ii] = blk[ii];
		lu_init(blk_n,blk_n,blk_lu,ipiv);
		db->dense_blk[li] = blk;
		db->dense_lu[li] = blk_lu;
		db->ipiv[li] = ipiv;
	}
	    	
	for(ii = 0; ii < nnei*8; ii++) delete[] temp[ii];
	delete[] temp;
	delete[] t;

	// debug: check FMMevalHB
/*	double** blk1 = new double*[840];
	for(int i = 0; i < 840; i++){
		blk1[i] = new double[840];
		int si = hb->waveletNodes[hb->wavlist[i]];
		for(int j = 0; j < 840; j++){
			FMMcount = 840;
			FMMsj[j] = hb->waveletNodes[hb->wavlist[j]];
		}	
		FMMevalHB(hb,FMMsj,FMMcount,si,eval,FMMtemp,nei);
		for(int j = 0; j < 840; j++){
			blk1[i][j] = eval[j];
		}
	}
	writeMatrix("finerBlk.dat",840,840,blk1);///
	for(int i = 0; i < 840; i++)
		delete[] blk1[i];
	delete[] blk1; 
*/
	hb->sortAvailableHB(minl-1); ///

	// collect diag_blks from high for each level and rearrange
	printf("trafo_diag_blk: collecting and rearranging...\n");///	
	//? This step requres O(n^2) memory

	delete[] FMMkk; //index in the neighborHBNodes list
	delete[] FMMsj; //index in the HB list
	delete[] eval;  //evaluation of transformed entries
	delete[] FMMtemp; //temparory memory for FMMevalHB
	delete[] site;
	delete[] rout;
	delete[] dist;

	cout << "trafo_diag_blk: Done. " << '\n';
}


void FWTd3::multilevel_trafo_diag_blk(Diag_Blk* db){

        printf(" In Multi Level trafo diag blk *************************\n");
	int q,qm;
	int MAX_LEVEL = 20;
        int ii,jj,kk,ll;
	int N = nei->N;
	int van = vanishingMoment;
        // int nnei = van*(van+1)*(van+2)/6;
	int nnei = van*(van+1)/2; //2D
	int* t = new int[MAX_LEVEL * N]; //temp memory for HB support 
	int tlen;  //length of it
	double val;  //temporary memory for a double
	double t0, t1;	
	bool coarseLifted = false;
        printf("HERE vanishingMoment = %d, pvanishingMoment = %d\n",vanishingMoment, pvanishingMoment);
	
	if(vanishingMoment > pvanishingMoment) coarseLifted = true;

        vector<int> levelstep = nei->levelstep;
        cout << "Trafo_diag_blk: levelstep initiated: \n";

	// determin the minimal lifted level
	int minl = minKwlevel + 1;

	// Make note here. Might cause problems with big d.
	double** temp = new double*[nnei*8];
	for(ii = 0; ii < nnei*8; ii++)
		temp[ii] = new double[MAX_LEVEL * N];
	int* hb2high = hb->hb2high;

	cout << "Trafo_diag_blk: Starting wavelet transform.." <<endl;	
        int treemark = nei->treesize-1;
	
	// May have a problem if coarse-lifting exitsts
	int wavMaxLevel = hb->wavMaxLevel;

        printf("  HERE wavMaxLevel = %d\n",wavMaxLevel);///
	if(coarseLifted) wavMaxLevel--;
        printf("  HERE wavMaxLevel = %d\n",wavMaxLevel);///

	while(treemark >= 0 && nei->tree[treemark][11] > wavMaxLevel)
		treemark--;

        int maxLevel = levelstep.size()-1;
        int li, lj, li_mark, lj_mark;
	int FMMcount;
	int* nodeslength = new int;
	int* nodes = new int[2 * N]; //index in the neighborHBNodes list
	int* FMMkk = new int[MAX_LEVEL * N]; //index in the neighborHBNodes list
	int* FMMsj = new int[MAX_LEVEL * N]; //index in the HB list
	double* eval = new double[MAX_LEVEL * N];  //evaluation of transformed entries
	int* FMMtemp = new int[MAX_LEVEL * N*2]; //temparory memory for FMMevalHB
	int currenttreemark;

	int pause;
	int* wavlist = hb->wavlist;
	double* blk;
	double* blk_lu;
	int* ipiv;
	cs* csblk;
	cs* cs_t;
	cs* csKw;

	printf("Min Level %d \n",minKwlevel);
	int finestcube_mark = treemark;
        printf("  HERE wavMaxLevel = %d\n",wavMaxLevel);///

	// Intialize all levels
	for (li = ( (wavMaxLevel+1)*(wavMaxLevel+2)/2 - 1) + 1; li >= minl-1; li--)
	  hb->sortAvailableHB(li);

	// Reserve memory for sparse matrix
	csKw = cs_spalloc(N,N,N*100,1,1); // Sparse Matrix Kw to store all blocks

	// Initialize 
	// Sweep Rows from wavMaxLevel to 0 for the finest column level wavMaxLevel
	lj = wavMaxLevel;
	hb->sortAvailableHB(lj + 1);

	for (li = lj; li >= minl-1; li--){
	  // Operating on columns
	  cout <<"multi_trafo_diag_blk: Rows: treemark at the end of level("<<li <<"): "<<treemark <<endl;
	  printf("Treemark Level %d HB Size %d Tree Size % d \n",treemark,li,hb->size,nei->treesize);
	  hb->sortAvailableHB(li);
	  // Process Columns
	  t0 = dsecnd();
	  Kwblockmatrixcolmninit(lj, li, treemark, nnei, temp, FMMkk, FMMsj, eval, FMMtemp, t, nodes, nodeslength, wavMaxLevel);
	  t1 = dsecnd();
	  printf("Column Pyramid scheme used %lf seconds for Kw. Levels i %d Level j %d.\n",t1 - t0,li,lj);
	  
	  hb->sortAvailableHB(li);
	  printf("Leave sort \n");
	}

	for (li = lj; li >= minl-1; li--){
	  // operating on rows
	  cout <<"multi_trafo_diag_blk: Rows: treemark at the end of level("<<li <<"): "<<treemark <<endl;
	  printf("Treemark Level %d HB Size %d Tree Size % d \n",treemark,li,hb->size,nei->treesize);

	  // Process rows
	  t0 = dsecnd();
	  Kwblockmatrixrow(li, lj, treemark, nnei, temp, FMMkk, FMMsj, eval, FMMtemp, t, nodes, nodeslength, wavMaxLevel);
	  t1 = dsecnd();
	  printf("Rows Pyramid scheme used %lf seconds for Kw. Levels i %d Level j %d.\n",t1 - t0,li,lj);
	  
	  hb->sortAvailableHB(li);
	  printf("Leave sort \n");

	  // Store in Sparse Matrix
	  KwStore(csKw, li, lj, nodes, nodeslength, wavMaxLevel);
	}

	printf("Erase last block from column %d %d \n",li,lj);
	li = li + 1;
	q = (wavMaxLevel - lj);
	cs_spfree(hb->highTriple[li + q*(wavMaxLevel + 1) - (q - 1)*q/2]);	
	cs_spfree(hb->high[li + q*(wavMaxLevel + 1) - (q - 1)*q/2]);
		
	// Process the rest of the matrix Kw
	for (lj = wavMaxLevel - 1; lj >= minl - 1; lj--){

	  q  = (wavMaxLevel - lj); 
	  qm = (wavMaxLevel - (lj + 1) );
	  li = lj;
	  hb->sortAvailableHB(li + q*(wavMaxLevel + 1) - (q - 1)*q/2 );
	  t0 = dsecnd();
	  Kwblockmatrixcolmn(lj, li, treemark, nnei, temp, FMMkk, FMMsj, eval, FMMtemp, t, nodes, nodeslength, wavMaxLevel);
	  t1 = dsecnd();
	  printf("Column Pyramid scheme used %lf seconds for Kw. Levels i %d Level j %d.\n",t1 - t0,li,lj);
	  hb->sortAvailableHB(li + q*(wavMaxLevel + 1) - (q - 1)*q/2 );

	  cs_spfree(hb->highTriple[li + 1 + qm*(wavMaxLevel + 1) - (qm-1)*qm/2]);
	  cs_spfree(hb->high[li + 1 + qm*(wavMaxLevel + 1) - (qm-1)*qm/2]);


	  for (li = lj - 1; li >= minl-1; li--){
	    // operating on Columns
	    cout <<"multi_trafo_diag_blk: Rows: treemark at the end of level("<<li <<"): "<<treemark <<endl;
	    printf("Treemark Level %d HB Size %d Tree Size % d \n",treemark,li,hb->size,nei->treesize);
	    
	    q  = (wavMaxLevel - lj); 
	    hb->sortAvailableHB(li + q*(wavMaxLevel + 1) - (q - 1)*q/2 );

	    // Process colmns
	    t0 = dsecnd();
	    Kwblockmatrixcolmninit(lj, li, treemark, nnei, temp, FMMkk, FMMsj, eval, FMMtemp, t, nodes, nodeslength, wavMaxLevel);
	    t1 = dsecnd();
	    printf("Colns Pyramid scheme used %lf seconds for Kw. Levels i %d Level j %d.\n",t1 - t0,li,lj);
	    hb->sortAvailableHB(li + q*(wavMaxLevel + 1) - (q - 1)*q/2 );
	    printf("Leave sort \n");

	    cs_spfree(hb->highTriple[li + 1 + qm*(wavMaxLevel + 1) - (qm-1)*qm/2]);
	    cs_spfree(hb->high[li + 1 + qm*(wavMaxLevel + 1) - (qm-1)*qm/2]);

	  }

	  
	  for (li = lj; li >= minl-1; li--){
	    // operating on Rows
	    cout <<"multi_trafo_diag_blk: Rows: treemark at the end of level("<<li <<"): "<<treemark <<endl;
	    printf("Treemark Level %d HB Size %d Tree Size % d \n",treemark,li,hb->size,nei->treesize);
	    
	    q  = (wavMaxLevel - lj);
	    qm = (wavMaxLevel - (lj + 1) );
	    hb->sortAvailableHB(li + q*(wavMaxLevel + 1) - (q - 1)*q/2 );

	    // Process rows
	    t0 = dsecnd();
	    Kwblockmatrixrow(li, lj, treemark, nnei, temp, FMMkk, FMMsj, eval, FMMtemp, t, nodes, nodeslength, wavMaxLevel);
	    t1 = dsecnd();
	    printf("Row Pyramid scheme used %lf seconds for Kw. Levels i %d Level j %d.\n",t1 - t0,li,lj);
	    hb->sortAvailableHB(li + q*(wavMaxLevel + 1) - (q - 1)*q/2 );
	    printf("Leave sort \n"); 

	    // Store in Sparse Matrix
	    KwStoreInt(csKw, li, lj, nodes, nodeslength, wavMaxLevel);
	  }

	  // Delete Old Unsed Matrix
	  printf("Erase last block %d %d \n",li,lj);
	  li  = 0;

	  q = (wavMaxLevel - lj);
	  cs_spfree(hb->highTriple[li + q*(wavMaxLevel + 1) - (q - 1)*q/2]);
	  cs_spfree(hb->high[li + q*(wavMaxLevel + 1) - (q - 1)*q/2]);

	}

       // Save matrix to a file
       char fname[200];
       sprintf(fname,"/home/julio/Projects/Data/Kriging/Kw.dat");
       writecsMatrix(fname,csKw);
       cs_spfree(csKw);

	for(ii = 0; ii < nnei*8; ii++) delete[] temp[ii];
	delete[] temp;
	delete[] t;
	delete nodeslength;

	delete[] FMMkk; //index in the neighborHBNodes list
	delete[] FMMsj; //index in the HB list
	delete[] eval;  //evaluation of transformed entries
	delete[] FMMtemp; //temparory memory for FMMevalHB
	cout << "trafo_diag_blk: Done. " << '\n';
}

// Use level = -1 for full transform on all levels.
// If level >=0, only the wavelet part at the corresponding level will be computed
//	and scaling function information will be ignored.
void FWTd3::trafo(double* d,int direction, int level){
	int ii,jj,kk,ll;
	int N = nei->N;
	int van = vanishingMoment;
	int pvan = pvanishingMoment;
	// int nnei = van*(van+1)*(van+2)/6;
	int nnei = van*(van+1)/2; //2D

	// If pvan < van, 
	// the finest level is hb->wavMaxlevel-1
	// hb->wavMaxLevel is a shaddow level that coexists with the coarsest level
	if (pvan < van && level == hb->wavMaxLevel)
		level = -1;  // necessary since wavelets on the shaddow level are 
			     // indeed lifted scaling functions

	vector<int> levelstep = nei->levelstep;

	int minl = 1;
	double* temp = new double[nnei*8];
	int* hb2high = hb->hb2high;
	int* wavlist = hb->wavlist;
	int* waveletIndex = hb->waveletIndex;
	int maxLevel = levelstep.size() - 1;
	
	double* t = new double[N]; //temperary storage for the transformed vector
	for(ii = 0; ii < N; ii++) t[ii] = 0;
	if(direction == FWT_BACKWARD){
	  // printf( " Starting inverse wavelet transform for level %d..\n",level);
	   int treemark;
	   if(level>=0) 
		treemark = levelstep[level];
	   else
		treemark = 0;
	   for(int li = level; li < maxLevel; li++){
	     // cout <<"treemark at the beginning of level("<<li <<"): "<<treemark <<endl;
  		while(treemark < levelstep[maxLevel] && nei->tree[treemark][11] == li){
			if(hb->linklength[treemark] > nnei){
			  int sj,sk;
			  for(jj = 0; jj < hb->linklength[treemark]; jj++){
				sj = hb->treelink[treemark][jj];
		      		temp[jj] = d[waveletIndex[sj]];
				d[waveletIndex[sj]] = 0;
			  }
			  // clear out the scaling function coefs on specified level
			  if(level>=0 && li==level)
				for(jj = 0; jj < nnei; jj++)
					temp[jj] = 0;
			  if(level < 0 || li == level || hb->linklength[treemark] < nnei)
			    for(jj = 0; jj < hb->linklength[treemark]; jj++){
				sj = hb->treelink[treemark][jj];
				for(kk = 0; kk < hb->blength[sj]; kk++){
					sk = hb->treelink[treemark][kk];
					d[waveletIndex[sk]] += hb->coef[sj][kk]*temp[jj];
				}
				for(kk = 0; kk < hb->blength[sj]; kk++){
					sk = hb->treelink[treemark][kk];
					t[waveletIndex[sk]] = d[waveletIndex[sk]];
				}
			    }
			  else
			    // Only necessary to resolve scaling functions to the finest scale if level>=0
			    for(jj = 0; jj < nnei; jj++){
				sj = hb->treelink[treemark][jj];
				for(kk = 0; kk < hb->blength[sj]; kk++){
					sk = hb->treelink[treemark][kk];
					d[waveletIndex[sk]] += hb->coef[sj][kk]*temp[jj];
				}
				for(kk = 0; kk < hb->blength[sj]; kk++){
					sk = hb->treelink[treemark][kk];
					t[waveletIndex[sk]] = d[waveletIndex[sk]];
				}
			    }
			}
			treemark++;
	        }
	   }
	   for(ii = 0; ii < N; ii++) d[hb->hb2high[hb->waveletNodes[hb->wavlist[ii]]]] = t[ii];
	}
	else if(direction == FWT_FORWARD){
	  //	   printf(" Starting forward wavelet transform for level %d..\n",level);
	   int treemark = nei->treesize - 1;
	   for(int li = maxLevel - 1; li >= level && li >= 0; li--){
	     //		cout <<"treemark at the end of level(" <<li <<"): "<<treemark <<endl;
		while(treemark >=0 && nei->tree[treemark][11] == li){
			if(hb->linklength[treemark] > nnei){
			  // si is the index to the first HB nodes linked to the treenode
			  int si = hb->treelink[treemark][0];
			  // Assuming all HB nodes linked to a treenode have the same support
 			   for(jj = 0; jj < hb->blength[si]; jj++)
			     temp[jj] = d[hb2high[hb->treelink[treemark][jj]]];
			   if(level < 0 || li == level || hb->linklength[treemark] < nnei){
			     for(jj = 0; jj < hb->linklength[treemark] && jj < nnei; jj++){
				si = hb->treelink[treemark][jj];
				d[hb2high[si]] = 0;
				for(kk = 0; kk < hb->blength[si]; kk++)
					d[hb2high[si]] += hb->coef[si][kk]*temp[kk];
				// not necessarry to record scaling functions unless level<0
				if(level < 0) 
					t[hb2high[si]] = d[hb2high[si]];
			     }
			     for(jj = nnei; jj < hb->linklength[treemark]; jj++){
				si = hb->treelink[treemark][jj];
				d[hb2high[si]] = 0;
				for(kk = 0; kk < hb->blength[si]; kk++)
					d[hb2high[si]] += hb->coef[si][kk]*temp[kk];
				// always record the wavelets that are lifted
				t[hb2high[si]] = d[hb2high[si]];
			     }
			   }else
			     // only necessary to evaluate the scaling function on finer levels for level>=0
			     for(jj = 0; jj < nnei; jj++){
				si = hb->treelink[treemark][jj];
				d[hb2high[si]] = 0;
				for(kk = 0; kk < hb->blength[si]; kk++)
					d[hb2high[si]] += hb->coef[si][kk]*temp[kk];
			     }	
			}
			treemark--;
		}
	   }

	   for(ii = 0; ii < N; ii++)
		d[ii] = t[hb2high[hb->waveletNodes[wavlist[ii]]]];
	}
	delete[] t;
}

double FWTd3::calcLiftingError(double** vande,vector<int>& phi,vector<double>& phi1,double* sii, int m, int k, int n){
	//  vande: m x k
	//  phi: k x n ( (n+1) x 1 actually )
	//  phi1: k x 1
	//  sii: n x 1
	
        // cout <<"FWTd3::calcLiftingError: Entering ..\n";
        // cout <<"m = " <<m <<"\nn = " <<n <<"\nk = "<<k <<'\n';
	double y = 0;
	double y1 = 0;
	double* phi_sii = new double[k];
	double* vande_phi = new double[m];
	for (int ii = 0; ii < k; ii++){
		phi_sii[ii] = 0;
	}
	//	cout <<"1\n";
	for (int ii = 0; ii < n; ii++){
		for (int jj = phi[ii]; jj < phi[ii+1]; jj++)
			phi_sii[jj] += sii[ii];
	}
	//	cout <<"2\n";
	for (int ii = 0; ii < k; ii++){
		phi_sii[ii] -= phi1[ii];
	}
	//	cout <<"3\n";
	for (int ii = 0; ii < m; ii++){
		vande_phi[ii] = 0;
		for (int jj = 0; jj < k; jj++)
			vande_phi[ii] += vande[ii][jj]*phi_sii[jj];
	}
	//	cout <<"4\n";
	for (int ii = 0; ii < m; ii++)
		y += vande_phi[ii]*vande_phi[ii];
	//	cout <<"5\n";
	int jj = 0;
	while (phi1[jj] == 0 ) jj++;
	//	cout <<"6\n";
	//	cout <<" jj: "<<jj <<endl;
	for (int ii = 0; ii < m; ii++)
		y1 += vande[ii][jj]*vande[ii][jj];
	y = sqrt(y/y1);
	//	cout <<"Ready to return :" <<y <<endl;
	delete[] phi_sii;
	delete[] vande_phi;
	//	cout <<" junk deleted" <<endl;
	return y;	
}

void FWTd3::precond(int oi, double**& p, double*& T,int&m){
	cout <<"Entering precond" <<endl;///
	
	int n = nei->N;
	m = n-oi;
	T = new double[m];
	for(int ii = 0; ii < m; ii++){
		if(high[ii][ii]>0)
			T[ii]=1/sqrt(high[ii][ii]);
		else T[ii] = 1;
	}
	//	cout <<"Preconditioner ready.." <<endl;///
	p = new double*[m];
	for (int ii = 0; ii < m; ii++)
		p[ii] = new double[m];
	for (int ii = 0; ii < m; ii++)
		for (int jj = 0; jj < m; jj++)
			p[ii][jj] = high[ii][jj]*T[ii];
	for (int ii = 0; ii < m; ii++)
		for (int jj = 0; jj < m; jj++)
			p[jj][ii] = p[jj][ii]*T[ii];
	//	cout <<"  T: " <<T <<endl;///
	//	cout <<"  p: " <<p <<endl;///
	//	writeMatrix("p.dat",m,m,p);///
}

// L_2 distance between cubes
double FWTd3::wavdist3(int i, int j){
	double r = 0;
	double bi[6],bj[6];
	double t;
	double ri, rj;
	double d;
	for(int ii = 0; ii < 6; ii++) bi[ii] = nei->bound[hb->cubenode[hb->waveletNodes[i]]][ii];
	for(int jj = 0; jj < 6; jj++) bj[jj] = nei->bound[hb->cubenode[hb->waveletNodes[j]]][jj];

	ri = 0;
	rj = 0;
	for (int ii = 0; ii < 3; ii++){
		t = bi[ii*2+1] - bi[ii*2];
		ri += t * t;  //diameter of B_a
		t = bj[ii*2+1] - bj[ii*2];
		rj += t * t;  //diameter of B_b
		t = (bi[ii*2+1] + bi[ii*2]);
		t -= (bj[ii*2+1] + bj[ii*2]);
		r += t * t;  //twice the distance between centers squared
	}
	d = sqrt(r) - sqrt(ri) - sqrt(rj);
	if (d < 0) d = 0;
	return d * 0.5;
}

// L_inf distance between cubes
double FWTd3::wavdist3Linf(int i, int j){
	double r[3];
	double bi[6],bj[6];
	double t[4];
	double d;
	for(int ii = 0; ii < 6; ii++) bi[ii] = nei->bound[hb->cubenode[hb->waveletNodes[i]]][ii];
	for(int jj = 0; jj < 6; jj++) bj[jj] = nei->bound[hb->cubenode[hb->waveletNodes[j]]][jj];

	for (int ii = 0; ii < 3; ii++){
		//min_i - min_j
		t[0] = bi[ii*2] - bj[ii*2];
		//min_i - max_j
		t[1] = bi[ii*2] - bj[ii*2+1];
		//max_i - min_j	
		t[2] = bi[ii*2+1] - bj[ii*2];
		//max_i - max_j
		t[3] = bi[ii*2+1] - bj[ii*2+1];
		if (t[0]*t[1] <= 0 || t[2]*t[3] <= 0
				|| t[0]*t[2] <= 0 || t[1]*t[3] <=0 )
			r[ii] = 0;
		else{
			r[ii] = t[1]>0 ? t[1]:-t[2];
			if (r[ii]<0) printf("Error in wavdist3Linf!\n");
		}
	}
	d = r[0];
	for(int ii = 1; ii < 3; ii++)
		if (d < r[ii]) d = r[ii];
	return d;
}

//L_2 distance between expanded neighbors
double FWTd3::wavdistLift3(int i, int j){
	double r[3];
	double bi[6],bj[6];
	double t[4];
	for(int ii = 0; ii < 6; ii++) bi[ii] = nei->bound[nei->treenode[i]][ii];
	for(int jj = 0; jj < 6; jj++) bj[jj] = nei->bound[nei->treenode[j]][jj];

	for(int ii = 0; ii < slength[ii]; ii++){
		for (int jj = 0; jj < 6; jj+=2)
			if(bi[jj]>nei->bound[nei->treenode[nei->neighbors[i][ii]]][jj])
				bi[jj] = nei->bound[nei->treenode[nei->neighbors[i][ii]]][jj];
		for (int jj = 1; jj < 6; jj+=2)
			if(bi[jj]>nei->bound[nei->treenode[nei->neighbors[i][ii]]][jj])
				bi[jj] = nei->bound[nei->treenode[nei->neighbors[i][ii]]][jj];
		for (int jj = 0; jj < 6; jj+=2)
			if(bj[jj]>nei->bound[nei->treenode[nei->neighbors[j][ii]]][jj])
				bj[jj] = nei->bound[nei->treenode[nei->neighbors[j][ii]]][jj];
		for (int jj = 1; jj < 6; jj+=2)
			if(bj[jj]>nei->bound[nei->treenode[nei->neighbors[j][ii]]][jj])
				bj[jj] = nei->bound[nei->treenode[nei->neighbors[j][ii]]][jj];

	}
	for (int ii = 0; ii < 3; ii++){
		t[0] = bi[ii*2] - bj[ii*2];
		t[1] = bi[ii*2] - bj[ii*2+1];
		t[2] = bi[ii*2+1] - bj[ii*2];
		t[3] = bi[ii*2+1] - bj[ii*2+1];
		if (t[0]*t[1] <= 0 || t[2]*t[3] <= 0
				|| t[0]*t[2] <= 0 || t[1]*t[3] <=0 )
			r[ii] = 0;
		else{
			r[ii] = fabs(t[0]);
			for (int jj = 1; jj < 4; jj++)
				if(r[ii] > fabs(t[jj])) r[ii] = fabs(t[jj]);
		}
	}
	return sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
}

void FWTd3::multiLift(int cube,int van){
	// printf("Entering MultiLift (cube = %d, van=%d)\n",cube,van);///
	// int m = van*(van+1)*(van+2)/6;
	int m = van*(van+1)/2; //2D
	
	// printf("   Accessing nei->tree[%d]\n",cube);///
	if(nei->tree[cube][0] == nei->tree[cube][1]){
		int pt = -nei->tree[cube][0]-1;
		double coef = 1.0;
		// printf("Claiming Vande3...\n");///
		//-pt-1 is the index to the actual center
		Vande3::get_vandemonde3(&(nei->coord[-pt-1]),van,1,moments_y);
		for(int ii = 0; ii < m; ii++) mome[ii] = moments_y[ii][0];
		hb->addHBNode(cube, &pt, &coef, 1, false,mome,m);
		return ;	
	}
	// printf("MultiLift: A\n");///
	int scalingFuncSize = 0;
	for(int ii = 2; ii < 10; ii++){
		if (nei->tree[cube][ii] > 0){
			multiLift(nei->tree[cube][ii],van);
			if(hb->linklength[nei->tree[cube][ii]] < m)
				scalingFuncSize += hb->linklength[nei->tree[cube][ii]];
			else scalingFuncSize += m;
		}
	}
	// printf("MultiLift: B\n");///
	if(m < scalingFuncSize){
		double** vpt = new double*[scalingFuncSize];
		int* support = new int[scalingFuncSize];
		int counter = 0;
		for(int ii = 2; ii < 10; ii++)
			if(nei->tree[cube][ii]>0)
			for(int jj = 0; jj < hb->linklength[nei->tree[cube][ii]] && jj < m; jj++){
				support[counter] = hb->treelink[nei->tree[cube][ii]][jj]; 
				vpt[counter] = hb->moments[support[counter]];
				counter++;
			}
		if(counter!=scalingFuncSize) printf("  Lifting Error: scalingFuncSize Mismatch!!\n");///
		double** vp = new double*[m];
		for(int ii = 0; ii < m; ii++){
			vp[ii] = new double[scalingFuncSize];
			for(int jj = 0; jj < scalingFuncSize; jj++)
				vp[ii][jj] = vpt[jj][ii];
		}
		double* y = new double[m];
		double** vt = new double*[scalingFuncSize];
		for(int ii = 0; ii < scalingFuncSize; ii++)
			vt[ii] = new double[scalingFuncSize];
		double** u = NULL;
		interfacedgesvd(vp,m,scalingFuncSize,y,u,vt);
		double* moment = new double[m];
		double err;
		for(int ii = 0; ii < m; ii++){
			for(int kk = 0; kk < m; kk++)
				moment[kk] = 0;
			// compute moments for the scaling function
			for(int jj = 0; jj < scalingFuncSize; jj++)
				for(int kk = 0; kk < m; kk++)
					moment[kk] += vt[ii][jj]*vp[kk][jj];
			hb->addHBNode(cube, support, vt[ii], scalingFuncSize, false,
				       moment,m);	
		}
		for(int ii = m; ii < scalingFuncSize; ii++){
			for(int kk = 0; kk < m; kk++)
				moment[kk] = 0;
			err = 0;
			// compute lifting error for the wavelets
			for(int jj = 0; jj < scalingFuncSize; jj++)
				for(int kk = 0; kk < m; kk++)
					moment[kk] += vt[ii][jj]*vp[kk][jj];
			for (int kk = 0; kk < m; kk++)
				err += moment[kk]*moment[kk];
			err = sqrt(err);
			if(liftingError<err) liftingError = err;
			hb->addHBNode(cube, support, vt[ii], scalingFuncSize, true,
					NULL,0);
		}
		delete[] y;
		for(int ii = 0; ii < scalingFuncSize; ii++) delete[] vt[ii];
		delete[] vt;
		delete[] moment;
		delete[] vpt;
		for(int ii = 0; ii < m; ii++) delete[] vp[ii];
		delete[] vp;
		delete[] support;
	}
	else{
		hb->linklength[cube] = scalingFuncSize;
		hb->treelink[cube] = new int[scalingFuncSize];
		int counter = 0;
		for(int ii = 2; ii < 10; ii++)
			if(nei->tree[cube][ii]>0)
			for(int jj = 0; jj < hb->linklength[nei->tree[cube][ii]]; jj++){
				hb->treelink[cube][counter] = hb->treelink[nei->tree[cube][ii]][jj];
				counter++;
			}
		return ;
	}
	// printf("MultiLift: Done\n",cube,van);///
	return ;
}

double FWTd3::coarseLift(int cube, int van, int pvan){
// van is the wavelet vanishing moment,
// pvan is the polynomial degree for the coarse lifting (pvan < van)
	printf("Entering coarseLift (cube = %d, van = %d, pvan = %d)\n",cube,van,pvan);///
	hb->primary_waveletSize = hb->waveletSize;
	if (pvan >= van) {
		printf("No coarse lifting to do. Done.\n");
		return 0;
	}

	int m = van*(van+1)/2; //2D
	int pm = pvan*(pvan+1)/2; //2D

	Vande3* vande3 = new Vande3;
	double** vv = vande3->vandemonde3(nei->coord,pvan,nei->N);
	int n = nei->N;
	int wn = hb->primary_waveletSize;
	int scalingFuncSize =  n - wn;
	int baseSize = hb->linklength[cube]; //Size of the bases for the root level
	printf("scalingFuncSize = %d, pm = %d, linklength = %d\n",
		scalingFuncSize,pm,hb->linklength[cube]);///
	printf("waveletSize at the beginning of coarselift is %d\n",
		hb->waveletSize);///

	if(pm < scalingFuncSize){
        	for(int i = 0; i < pm; i++){
                	trafo(vv[i],FWT_FORWARD,-1);  // transform to level 0 is optimal
							// -1 for now
        	}
		double** vp = new double*[pm];
		for(int ii = 0; ii < pm; ii++){
			vp[ii] = new double[scalingFuncSize];
			for(int jj = 0; jj < scalingFuncSize; jj++)
				vp[ii][jj] = vv[ii][wn+jj];
		}
		double* y = new double[pm];
		double** vt = new double*[scalingFuncSize];
		for(int ii = 0; ii < scalingFuncSize; ii++)
			vt[ii] = new double[scalingFuncSize];
		double** u = NULL;
		printf("Starting coarse lift...\n");
		interfacedgesvd(vp,pm,scalingFuncSize,y,u,vt);
		printf("SVD done.\n");
		double** moment = new double*[scalingFuncSize];
		double* mome = new double[pm];
		for(int ii = 0; ii < scalingFuncSize; ii ++)
			moment[ii] = new double[baseSize];
		// Assuming that the first scalingFuncSize HB nodes are old scaling functions
		double err;
		// printf("Computing moments ...\n");///
		int shift = 0;
		for(int ii = 0; ii < scalingFuncSize; ii++){
			// recompute coefs for the scaling function
			for(int jj = 0; jj < baseSize; jj++){
				moment[ii][jj] = 0;
				for(int kk = 0; kk < scalingFuncSize; kk++)
					// a shift is not necessary since we are only
					// lifting the scaling function part and it is 
					// at the beginning of the treelink[cube]
					moment[ii][jj] += vt[ii][kk] * hb->coef[hb->treelink[cube][shift+kk]][jj];
			}
		}
		// printf("New moments computed.\n");///
		shift = scalingFuncSize - pm;
		// Shift scaling functions to the end
		for(int ii = 0; ii < pm; ii++)
			for(int jj = 0; jj < baseSize; jj++)
				hb->coef[hb->treelink[cube][ii+shift]][jj] = moment[ii][jj];
		for(int ii = pm; ii < scalingFuncSize; ii++){
			for(int jj = 0; jj < baseSize; jj++)
				hb->coef[hb->treelink[cube][ii-pm]][jj] = moment[ii][jj];
			hb->wavelet[hb->treelink[cube][ii-pm]] = true;
		}
		hb->waveletSize += shift;
		// printf("Wavelet nodes adjusted.\n");///
		// Check lifting error
		for(int ii = pm; ii < scalingFuncSize; ii++){
			for(int kk = 0; kk < pm; kk++)
				mome[kk] = 0;
			err = 0;
			// compute lifting error for the wavelets
			for(int jj = 0; jj < scalingFuncSize; jj++)
				for(int kk = 0; kk < pm; kk++)
					mome[kk] += vt[ii][jj]*vp[kk][jj];
			for (int kk = 0; kk < pm; kk++)
				err += mome[kk]*mome[kk];
			err = sqrt(err);
			if(liftingError<err) liftingError = err;
		}
		printf("Coarse lift done.\n");///
		delete[] y;
		for(int ii = 0; ii < scalingFuncSize; ii++){
			delete[] vt[ii];
			delete[] moment[ii];
		}
		delete[] vt;
		delete[] moment;
		delete[] mome;
		for(int ii = 0; ii < pm; ii++) delete[] vp[ii];
		delete[] vp;
		// printf("Leaving coarse lift.\n");///
	}
	else{
		printf("Not enough centers for the specified polynomial degree.\n");///
	}

	// Assuming that the scaling functions are sorted to the very end
	int minl = hb->wavMaxLevel;
	hb->wavlevels[minl+1][0] = hb->primary_waveletSize;
	hb->wavlevels[minl+1][1] = hb->waveletSize;
	minl++;
	hb->wavMaxLevel = minl;
	for(int ii = 0; ii <= minl; ii++)
		printf("level %d: %d - %d\n",ii,hb->wavlevels[ii][0],hb->wavlevels[ii][1]);///

	for(int i = 0; i < pm; i++) delete[] vv[i];
	delete[] vv;
	delete vande3; //consider making it static ?
	return liftingError;
}

void FWTd3::wavMask(double**& w){
        int mlen;
        int* supp;
        double* mask;
        int N = nei->N;
        w = new double*[N];
        for(int ii = 0; ii < N; ii++){
                w[ii] = new double[N];
                for(int jj = 0; jj < N; jj ++)
                        w[ii][jj] = 0;
        }
        for(int ii = 0; ii < N; ii++){
	  // printf("Trying to genMask for HB node %d(%d)\n",hb->waveletNodes[ii],hb->blength[hb->waveletNodes[ii]]);///
                hb->genMask(hb->waveletNodes[hb->wavlist[ii]],mask,supp,mlen);
		// printf(".. Done. mlen=%d\n",mlen);
                for(int jj = 0; jj < mlen; jj ++)
                        w[ii][supp[jj]] = mask[jj];
                delete[] mask;
                delete[] supp;
        }
}

void FWTd3::wavMask(double**& mask_out, int**& msupport_out, int*& mlength_out){

// Adding this as a way to generate the sparse wavelet coef matrix
// Idea: apply forward wavelet transform to a sparse identity matrix

	printf("Entering wavMask...\n");

	int ii,jj,kk,ll;
	int N = nei->N;
	int van = vanishingMoment;
	int pvan = pvanishingMoment;
	int nnei = van*(van+1)/2; //2D


	vector<int> levelstep = nei->levelstep;
	int minl = 1;
	double* temp = new double[nnei*8];
	int* hb2high = hb->hb2high;
	int* wavlist = hb->wavlist;
	int* waveletIndex = hb->waveletIndex;
	int maxLevel = levelstep.size() - 1;

	// pointers from actual support to the indexed(compressed) support
	int* tindex = new int[N];
	double* tmask;
	int* tsupp;
	int tlength;

	// allocate space for the output variables in sorted order
	mask_out = new double*[N];
	msupport_out = new int*[N];
	mlength_out = new int[N];

	// allocate space for the output variabls in unsorted order
	double** mask = new double*[N];
	int** msupport = new int*[N];
	int* mlength = new int[N];
	for(int ii = 0; ii < N; ii++){
		mask[ii] = NULL;
		msupport[ii] = NULL;
		mlength[ii] = 0;
	}

	//pointers to scaling functions before lifting
	int** supp0 = new int*[nnei*8];
	double** mask0 = new double*[nnei*8];
	int* length0 = new int[nnei*8];

	// Need to initialize all the mask, msupport, mlength to appropriate values
	// so that they can be interpreted correctly for finest hb nodes.
	

	int treemark = nei->treesize - 1;
	// computing the wavelet coefficient all the way to 
	// the coarsest level
	for(int li = maxLevel - 1; li >= 0; li--){
		cout <<"treemark at the end of level(" <<li <<"): "<<treemark <<endl;
		while(treemark >=0 && nei->tree[treemark][11] == li){
			if(hb->linklength[treemark] > nnei){
			// si is the index to the first HB nodes linked to the treenode
			     int si = hb->treelink[treemark][0];
			     // initialize the temporary holders for the finer
			     // level HB nodes
			     for(jj = 0; jj < hb->linklength[treemark];jj++){
				int sj = hb->treelink[treemark][jj];
				supp0[jj] = msupport[hb2high[sj]];
				mask0[jj] = mask[hb2high[sj]];
				length0[jj] = mlength[hb2high[sj]];
			     }
			     // Assuming all HB nodes linked to a treenode have the same support
			     for(jj = 0; jj < hb->linklength[treemark]; jj++){
				int sj = hb->treelink[treemark][jj];
				// Just collecting the coeffcients for node sj
				getsupport(sj,tsupp,tlength,tindex);
				// it maybe more efficient to copy support from si
				tmask = new double[tlength];
				for(kk = 0; kk < tlength; kk++)
					tmask[kk] = 0;
				for(kk = 0; kk < hb->blength[sj]; kk++){
					int sk = hb->treelink[treemark][kk];
					if(mask0[kk])
					// when mask has already been generated for the finer level node
					  for(ll = 0; ll < length0[kk]; ll++)
					// Adding the contribution from finer level hb nodes
						tmask[tindex[supp0[kk][ll]]] += hb->coef[sj][kk]*mask0[kk][ll];
					else
					//when the finest level nodes are reached
						tmask[tindex[hb2high[sk]]] += hb->coef[sj][kk];
				}
				// Replace records for the hb nodes mapped by hb2high.
 			        mlength[hb2high[sj]] = tlength;
				mask[hb2high[sj]] = tmask;
				msupport[hb2high[sj]] = tsupp;
			     }
			     for(jj = 0; jj < hb->linklength[treemark]; jj++)
				if(mask0[jj]){
					delete[] mask0[jj];
					delete[] supp0[jj];
				}
			}
			treemark--;
		}
	}
	for(ii = 0; ii < N; ii++){
		mlength_out[ii] = mlength[hb2high[hb->waveletNodes[wavlist[ii]]]];
		mask_out[ii] = mask[hb2high[hb->waveletNodes[wavlist[ii]]]];
		msupport_out[ii] = msupport[hb2high[hb->waveletNodes[wavlist[ii]]]];
	}
	printf("wavMask generated. Done.\n");///
}

void FWTd3::hbMask(double**& mask_out, int**& msupport_out, int*& mlength_out){

// Adding this as a way to generate the sparse wavelet coef matrix
// Idea: apply forward wavelet transform to a sparse identity matrix

	printf("Entering hbMask...\n");

	int ii,jj,kk,ll;
	int N = nei->N;
	int van = vanishingMoment;
	int pvan = pvanishingMoment;
	int nnei = van*(van+1)/2; //2D


	vector<int> levelstep = nei->levelstep;
	int minl = 1;
	double* temp = new double[nnei*8];
	int* hb2high = hb->hb2high;
	int* wavlist = hb->wavlist;
	int* waveletIndex = hb->waveletIndex;
	int maxLevel = levelstep.size() - 1;

	// pointers from actual support to the indexed(compressed) support
	int* tindex = new int[N];
	double* tmask;
	int* tsupp;
	int tlength;

	// allocate space for the output variables in sorted order
	mask_out = new double*[hb->size];
	msupport_out = new int*[hb->size];
	mlength_out = new int[hb->size];

	// allocate space for the output variabls in unsorted order
	double** mask = new double*[N];
	int** msupport = new int*[N];
	int* mlength = new int[N];
	for(int ii = 0; ii < N; ii++){
		mask[ii] = NULL;
		msupport[ii] = NULL;
		mlength[ii] = 0;
	}

	// pointers to scaling functions before lifting
	int** supp0 = new int*[nnei*8];
	double** mask0 = new double*[nnei*8];
	int* length0 = new int[nnei*8];

	//Need to initialize all the mask, msupport, mlength to appropriate values
	// so that they can be interpreted correctly for finest hb nodes.

	int treemark = nei->treesize - 1;
	// computing the wavelet coefficient all the way to 
	// the coarsest level
	for(int li = maxLevel - 1; li >= 0; li--){
		cout <<"treemark at the end of level(" <<li <<"): "<<treemark <<endl;
		while(treemark >=0 && nei->tree[treemark][11] == li){
			if(hb->linklength[treemark] > nnei){ //If lifted,
			// si is the index to the first HB nodes linked to the treenode
			     int si = hb->treelink[treemark][0];
			     // initialize the temporary holders for the finer
			     // level HB nodes
			     for(jj = 0; jj < hb->linklength[treemark];jj++){
				int sj = hb->treelink[treemark][jj];
				if(sj < 0) {
					printf("Negative index!\n");
					exit(-1);
				}
				supp0[jj] = msupport[hb2high[sj]];
				mask0[jj] = mask[hb2high[sj]];
				length0[jj] = mlength[hb2high[sj]];
			     }
			// Assuming all HB nodes linked to a treenode have the same support
			     for(jj = 0; jj < hb->linklength[treemark]; jj++){
				int sj = hb->treelink[treemark][jj];
				// Just collecting the coeffcients for node sj
				getsupport(sj,tsupp,tlength,tindex);
				// it maybe more efficient to copy support from si
				tmask = new double[tlength];
				// Initialize the mask to 0
				for(kk = 0; kk < tlength; kk++)
					tmask[kk] = 0;
				// loop through HB support of sj
				for(kk = 0; kk < hb->blength[sj]; kk++){
					int sk = hb->treelink[treemark][kk];
					if(mask0[kk])
					// when mask has already been generated for the finer level node
					  for(ll = 0; ll < length0[kk]; ll++)
					// Adding the contribution from finer level hb nodes
						tmask[tindex[supp0[kk][ll]]] += hb->coef[sj][kk]*mask0[kk][ll];
					else
					// when the finest level nodes are reached
						tmask[tindex[hb2high[sk]]] += hb->coef[sj][kk];
				}
			// Replace records for the hb nodes mapped by hb2high.
 			        mlength[hb2high[sj]] = tlength;
				mask[hb2high[sj]] = tmask;
				msupport[hb2high[sj]] = tsupp;
			// Save the mask for hb nodes
				mlength_out[sj] = tlength;
				mask_out[sj] = tmask;
				msupport_out[sj] = tsupp;
			     }
			}else{ // If unlifted,
			       // Loop through the HB nodes
			     for(jj = 0; jj < hb->linklength[treemark];jj++){
				int sj = hb->treelink[treemark][jj];
				if(mask[hb2high[sj]]){
					// Here we may have multiple HB nodes 
					// linked to the same mask/support
					mlength_out[sj] = mlength[hb2high[sj]];
					mask_out[sj] = mask[hb2high[sj]];
					msupport_out[sj] = msupport[hb2high[sj]];
				}else{
					// We should have unlifted scaling functions here
					getsupport(sj,tsupp,tlength,tindex);
					if(tlength!=1){
						printf("Sth. wrong here. (tlength=%d)\n",tlength);
						exit(-1);
					}
					mlength_out[sj] = tlength;
					mask_out[sj] = new double[1];
					mask_out[sj][0] = 1;
					msupport_out[sj] = tsupp;
				}
			    }
			}
			treemark--;
		}
	}
	printf("hbMask generated. Done.\n");///
}

void FWTd3::getsupport(int si, int*& tsupp, int& tlength, int* tindex){
//si: index to HB node
//tsupp: returned support
//tlength: size of tsupp
//tindex: pointers from node index to tsupp index
	vector<int> y;
	int node = hb->cubenode[si];
	double** bound = nei->bound;
	int** tree = nei->tree;
	nei->getQuadrant3(bound[node][0],bound[node][1],bound[node][2],bound[node][3],bound[node][4],
			  bound[node][5],tree[node][0],tree[node][1],y);
	tlength = y.size();
	tsupp = new int[tlength];
	for(int i = 0; i < tlength; i++){
		tsupp[i] = y[i];
		tindex[y[i]] = i;
	}
	y.clear();
}

void FWTd3::findNeighborHBNodes(int d){
	int cube;
	double* bound;
	double* bound1 = new double[6];
	double* bound2;
	double stepX,stepY,stepZ;
	int treemark;
	int level;
	int van = vanishingMoment;
	int nnei = van*(van+1)/2; //2D

	int maxLength = (2*d + 1);
	maxLength = maxLength*maxLength*maxLength * hb->size; //maximum # of neighboring HB nodes
	int* nodes = new int[maxLength];
	int length;
	double eps = 2.2204E-16;
	printf("Entering findNeighborHBNodes...\n");///
	printf("tau = %d \n", d);///
	for(int i = 0; i < hb->size; i++){
		cube = hb->cubenode[i];
		level = nei->tree[cube][11];
		bound = nei->bound[cube];
		stepX = bound[1] - bound[0];
		stepY = bound[3] - bound[2];
		stepZ = bound[5] - bound[4];
		// printf("stepX %f stepY %f stepZ %f \n", stepX, stepY, stepZ);
		length = 0;
	     if(level>0)
		for(int ii = -d; ii < d+1; ii++)
		  for(int jj = -d; jj < d+1; jj++)
		    for(int kk = -d; kk < d+1; kk++)
		      // continue;
		    {
		      // printf("findNeighborHBNodes: testing HB node #%d (%d,%d,%d) cube # %d at level %d\n",
		      //		i,ii,jj,kk,cube,level);///
			bound1[0] = bound[0] + stepX*ii + eps;
			bound1[1] = bound[1] + stepX*ii - eps;
			bound1[2] = bound[2] + stepY*jj + eps;
			bound1[3] = bound[3] + stepY*jj - eps;
			bound1[4] = bound[4] + stepZ*kk + eps;
			bound1[5] = bound[5] + stepZ*kk - eps;
			// printf("bounds: [%f, %f, %f] - [%f, %f, %f]\n",bound1[0], bound1[2], 
			// bound1[4],bound1[1],bound1[3],bound1[5]);///
			treemark = 0;
			// efficiently locate the neighboring node
			// printf("nei->tree[treemark][11] = %d\n",nei->tree[treemark][11]);///
			int found = 0;
			while(nei->tree[treemark][11] < level){
				found = 0;
				for(int j = 2; j < 10; j++){
					if(nei->tree[treemark][j] > 0){
						bound2 = nei->bound[nei->tree[treemark][j]];
						if (nei->contained(bound1,bound2)){
							treemark = nei->tree[treemark][j];
							found = 1;
							break;
						}
					}
				}
				if (found == 0) break;
			}

			if(found && nei->tree[treemark][11] == level)
			  // printf("treemark: %d at level %d\n",treemark,nei->tree[treemark][11]);///
				for(int j = 0; j < hb->linklength[treemark]; j++){
				    // add the neighbor only if it is on the same level
				    // i.e. the highest level of its linked cube matches "level"
				    if(nei->tree[hb->cubenode[hb->treelink[treemark][j]]][11] == level){
					nodes[length] = hb->treelink[treemark][j];
					length++;
				    }
				}
		    }
	      else//for level 0	
		for(int j = 0; j < hb->linklength[0]; j++){
			nodes[length] = hb->treelink[0][j];
			length++;
		}
		if (length==0){
			printf("neighborhood failed!\n");
			printf("i = %d, level = %d\n",i,level);///
			exit(-1);
		}
		hb->neighborHBLength[i] = length;
		hb->neighborHBNodes[i] = new int[length];
		for(int j = 0; j < length; j++) {
			hb->neighborHBNodes[i][j] = nodes[j];
		}

	}
	printf("findNeighborHBNodes: Done.\n");///
	delete[] bound1;
	delete[] nodes;
}

void FWTd3::DebuggingfindMultiLevel(int d, int N){
// Test neighborhood

  int cube, level, node, length, out;
  int* nodeslength = new int;
  int* nodes = new int[2 * N]; //index in the neighborHBNodes list
  int treemark;

  printf("Tests ----------------------------------------------------------\n");

	for(int si = 0; si < hb->size; si++){
		cube = hb->cubenode[si];
		level = nei->tree[cube][11];
		length = hb->neighborHBLength[si];

		if (!(length==*nodeslength)) printf("First Test: Test broke Node %d Level %d Cube %d Original \
                                             Length %d Second Length %d \n",si,level,cube,length,*nodeslength);

		// Second Test
		int flag = 0;
		for(int j = 0; j < length; j++) {
		  node = hb->neighborHBNodes[si][j];
		  if (!(node == nodes[j])) printf("Second Test: Test broke j %d Cube %d Current Node %d Level %d \
                                                   Original Length %d Second Length %d Original Node %d Comparison Node %d  \n", 
						   j, cube, si, level, length, *nodeslength, node, nodes[j]);
		  flag = 1;
		}
	}
}




// Find Multi Level Neighborhood Nodes for the HB
void FWTd3::findMultiLevelNeighborHBNodes(int d){
        int MAX_LEVEL = 20;
	int cube;
	double* bound;
	double* bound1 = new double[6];
	double* bound2;
	double stepX,stepY,stepZ;
	int treemark;
	int level;
	int van = vanishingMoment;
        // int nnei = van*(van+1)*(van+2)/6;
	int nnei = van*(van+1)/2; //2D

	int maxLength = (2*d + 1);
	maxLength = maxLength*maxLength*maxLength * hb->size; //maximum # of neighboring HB nodes
	int* nodes = new int[maxLength];
	int length;
	double eps = 2.2204E-16;
	int wavMaxLevel = hb->wavMaxLevel;
	printf("Entering Multi Level findNeighborHBNodes...\n");///
	for(int i = 0; i < hb->size; i++){
		cube = hb->cubenode[i];
		level = nei->tree[cube][11];
		// printf("Cube %d Level %d \n", cube, level);
		bound = nei->bound[cube];
		stepX = bound[1] - bound[0];
		stepY = bound[3] - bound[2];
		stepZ = bound[5] - bound[4];
		length = 0;
	     if(level>0)
		for(int ii = -d; ii < d+1; ii++)
		  for(int jj = -d; jj < d+1; jj++)
		    for(int kk = -d; kk < d+1; kk++)
		    {
			bound1[0] = bound[0] + stepX*ii + eps;
			bound1[1] = bound[1] + stepX*ii - eps;
			bound1[2] = bound[2] + stepY*jj + eps;
			bound1[3] = bound[3] + stepY*jj - eps;
			bound1[4] = bound[4] + stepZ*kk + eps;
			bound1[5] = bound[5] + stepZ*kk - eps;
			// printf("bounds: [%f, %f, %f] - [%f, %f, %f]\n",bound1[0], bound1[2], \
			// bound[4],bound1[1],bound1[3],bound1[5]);///
			treemark = 0;
			// efficiently locate the neighboring node
			int found = 0;
			while(nei->tree[treemark][11] < level){
				found = 0;
				for(int j = 2; j < 10; j++){
					if(nei->tree[treemark][j] > 0){
						bound2 = nei->bound[nei->tree[treemark][j]];
						if (nei->contained(bound1,bound2)){
							treemark = nei->tree[treemark][j];
							found = 1;
							break;
						}
					}
				}
				if (found == 0) break;
				// printf("treemark: %d at level %d\n",treemark,nei->tree[treemark][11]);///
			}
			if(found && nei->tree[treemark][11] == level) {
				    // add the neighbor cube only if it is on the same level
				    // i.e. the highest level of its linked cube matches "level"
					nodes[length] = treemark;
					length++;
			}				
		    }
	      else//for level 0	
		{
			nodes[length] = 0;
			length++;
		}
		if (length==0){
			printf("neighborhood failed!\n");
			printf("i = %d, level = %d\n",i,level);///
			exit(-1);
		}
		hb->neighborHBCubesLength[i] = length;
		hb->neighborHBCubes[i] = new int[length];
		for(int j = 0; j < length; j++)
			hb->neighborHBCubes[i][j] = nodes[j];
	}


	delete[] bound1;
	delete[] nodes;
}



// Find Multi Level Neighborhood Nodes for the HB
// Compute transpose of multilevel distance criterion for i > j
void FWTd3::findMultiLevelNeighborHBNodesT(int d){

        // Store in hb->neighborHBCubesLength[cube] and hb->neighborHBCubes[cube] all neighborhood cubes linked functions
        // Note that we are redefining hb->neighborHBCubesLength[cube] and hb->neighborHBCubes[cube]
        int MAX_LEVEL = 20;
	// int cube;
	double* bound;
	double* bound1 = new double[6];
	double* bound2;
	double stepX,stepY,stepZ;
	int treemark;
	int level;
	int van = vanishingMoment;
	int nnei = van*(van+1)/2; //2D

	int maxLength = (2*d + 1);
	maxLength = maxLength*maxLength*maxLength * hb->size; //maximum # of neighboring HB nodes
	int* nodes = new int[maxLength];
	int length;
	double eps = 2.2204E-16;
	printf("Entering Multi Level findNeighborHBNodes...\n");///
	for(int cube = 0; cube < nei->treesize; cube++){
	        // cube = hb->cubenode[i];
		level = nei->tree[cube][11];
		// printf("Cube %d Level %d \n", cube, level);
		bound = nei->bound[cube];
		stepX = bound[1] - bound[0];
		stepY = bound[3] - bound[2];
		stepZ = bound[5] - bound[4];
		length = 0;
	     if(level>0)
		for(int ii = -d; ii < d+1; ii++)
		  for(int jj = -d; jj < d+1; jj++)
		    for(int kk = -d; kk < d+1; kk++)
		    {
			bound1[0] = bound[0] + stepX*ii + eps;
			bound1[1] = bound[1] + stepX*ii - eps;
			bound1[2] = bound[2] + stepY*jj + eps;
			bound1[3] = bound[3] + stepY*jj - eps;
			bound1[4] = bound[4] + stepZ*kk + eps;
			bound1[5] = bound[5] + stepZ*kk - eps;
			// printf("bounds: [%f, %f, %f] - [%f, %f, %f]\n",bound1[0], bound1[2], bound[4], \
			// bound1[1],bound1[3],bound1[5]);///
			treemark = 0;
			//efficiently locate the neighboring node
			int found = 0;
			while(nei->tree[treemark][11] < level){
				found = 0;
				for(int j = 2; j < 10; j++){
					if(nei->tree[treemark][j] > 0){
						bound2 = nei->bound[nei->tree[treemark][j]];
						if (nei->contained(bound1,bound2)){
							treemark = nei->tree[treemark][j];
							found = 1;
							break;
						}
					}
				}
				if (found == 0) break;
//				printf("treemark: %d at level %d\n",treemark,nei->tree[treemark][11]);///
			}
			if(found && nei->tree[treemark][11] == level)
				for(int j = 0; j < hb->linklength[treemark]; j++){
				    //add the neighbor only if it is on the same level
				    // i.e. the highest level of its linked cube matches "level"
				    if(nei->tree[hb->cubenode[hb->treelink[treemark][j]]][11] == level){
					nodes[length] = hb->treelink[treemark][j];
					length++;
				    }
				}
		    }
	      else//for level 0	
		for(int j = 0; j < hb->linklength[0]; j++){
			nodes[length] = hb->treelink[0][j];
			length++;
		}
		if (length==0){
		  // printf("No neighborhood \n");
		  // printf("cube = %d, level = %d\n",cube,level);///
		  hb->neighborCubesLength[cube] = length;
		  hb->neighborCubes[cube] = new int[1];
		  hb->neighborCubes[cube][0] = -1;
		  // exit(-1);
		}
		else{
		hb->neighborCubesLength[cube] = length;
		hb->neighborCubes[cube] = new int[length];
		// printf("Cube %d Length %d 2nd Length %d Level %d \n",cube,length,hb->neighborCubesLength[cube], level);
		for(int j = 0; j < length; j++){
			hb->neighborCubes[cube][j] = nodes[j];
			// printf("%d ",	hb->neighborCubes[cube][j]);
		}
		}
		// printf("\n");
	}
	
        printf("Compute i > j distance criterion neighborhoods \n");
	
	int multilevel = 0;
	int multilevelcounter = 0;
	int index = 0;

	for(int i = 0; i < hb->size; i++){
	    treemark = hb->cubenode[i];
	    level = nei->tree[treemark][11];
	    multilevel = level;
	    multilevelcounter = level;
	    hb->multilevelcount = 0;
	    
	    while(multilevelcounter > -1) {
	      length = 0;

	      // Store data in temporary multilevel array level
	      
	      hb->multilevelnodes[hb->multilevelcount] = multilevel; // 1st entry: level
	      hb->multilevelcount++;
	      hb->multilevelnodes[hb->multilevelcount] = treemark;   // 2nd entry: Treemark
	      hb->multilevelcount++;
	    
	      // get Parent and update multilevel
	      treemark = nei->tree[treemark][10];
	      // printf("Parent %d \n", treemark);
	      
	      multilevel = nei->tree[treemark][11];
	      multilevelcounter--;
	    }
	    hb->MutliLevelneighborHBNodesTranspose[i] = new int[hb->multilevelcount];
	    for(int f = 0; f < hb->multilevelcount; f++){
	      hb->MutliLevelneighborHBNodesTranspose[i][f] = hb->multilevelnodes[f];
	    }
	    hb->MutliLevelneighborHBNodesLengthTranspose[i] = hb->multilevelcount;
	}


	    printf("Multi Level distance criterion findNeighborHBNodes: Done.\n");///
	    delete[] bound1;
	    delete[] nodes;	    
}
void FWTd3::StoreLinkDataMultiLevel(int level, int multilevel, int i, int MAX_LEVEL, int counter, int *levelcubes, int *nodes){
  // Store scaling function neighborhs

  // Debug
  int debugq = 0;
  int valueq = 0;

            // Find scaling function neighborhs (at multilevel) for cube at level. 
            int totalcounter = 0;    
	    // printf("Nodes \n");
	    for(int q = 0; q < counter; q++){

	      for(int f = 0; f < hb->linklength[levelcubes[q]]; f++){
		// Only store neighbors functions at the multilevel level.
		if (nei->tree[hb->cubenode[hb->treelink[levelcubes[q]][f]]][11] == multilevel) {
		  nodes[totalcounter] = hb->treelink[levelcubes[q]][f];
		  totalcounter++;
		}
	      }
	    }
	    //printf("\n");
	    
	    // Store neighbors 
	    // Next 2 lines is to store the level and number of functions at that level
	    // This is used to find the neighbor functions
	    // Commented out put all the nodes together
	    hb->multilevelnodes[hb->multilevelcount] = multilevel; // 1st entry: level
	    hb->multilevelcount++;
	    hb->multilevelnodes[hb->multilevelcount] = totalcounter; // 2nd entry: length of nodes
	    hb->multilevelcount++;
	    
	    int shift = hb->multilevelcount;
	    // printf("Shift %d Nodes \n", hb->multilevelcount);
	    for(int f = 0; f < totalcounter; f++){
	      hb->multilevelnodes[f + shift] = nodes[f];
	      // printf("%d ", nodes[f]);
	      hb->multilevelcount++;
	    }
}

void FWTd3::multileveldistcritUpper(int level, int currentnode, int *nodes, int *nodeslength){

            // Define neighborhood between 
            // li = level of node si 
            // level >= li

            for(int ii = 0; ii < hb->MutliLevelneighborHBNodesLength[currentnode]; ii++){
	      // printf("%d ", hb->MutliLevelneighborHBNodes[currentnode][ii]);
	    }
	    // printf("\n");

            int locallength;
            int seek = 0;
	    int seeklevel = hb->MutliLevelneighborHBNodes[currentnode][seek];
	    seek++;
	    int seekcounter = hb->MutliLevelneighborHBNodes[currentnode][seek];
	    seek++;
	    while(seek > -1){
	      if (seeklevel == level) {
		locallength = 0;
		for(int ii = seek; ii < seek + seekcounter; ii++){
		  nodes[locallength] = hb->MutliLevelneighborHBNodes[currentnode][ii];
		  // printf("%d %d ", locallength, nodes[locallength]);
		  locallength++;
		}
		break;
	      }
	      seek = seek + seekcounter;
	      seeklevel = hb->MutliLevelneighborHBNodes[currentnode][seek];
	      seek++;
	      seekcounter = hb->MutliLevelneighborHBNodes[currentnode][seek];
	      seek++;
	    }
	    
	    *nodeslength = locallength;
}

void FWTd3::Kwblockmatrixcolmninit(int lj, int li, int treemark, int nnei, double **temp, int *FMMkk, int *FMMsj, double *eval, int *FMMtemp, int *t, int *nodes, int *nodeslength, int wavMaxLevel){

         int tlen;
	 // int currenttreemark;
	 int FMMcount;
	 int jj;
	 int ll;
	 double val;
	 int kk;
	 int q = (wavMaxLevel - lj );

	 // Declare cube placeholders to pass to distance criterion. Makes it faster.
	 int maxLengthcubes = 2 * nei->treesize;
	 int* levelcubes = new int[maxLengthcubes];
	 int* subdividedlevelcubes = new int[maxLengthcubes];	
	 int* hb2high = hb->hb2high;
	 
	 // Initialize vectors for Matern kernel summation
	 int nsite = MAX_DIRECT_KERNEL;
	 double* dist   = new double[nsite];
	 double* site   = new double[nsite];                  // Array of interpolation sites
	 double* rout   = new double[nsite];          // spline evaluation results

	 // loop through the cubes
	 printf("Columns Treemark %d lj %d\n",treemark, lj);

	 while (treemark > -1){
	   if (nei->tree[treemark][11] == li) {
	     int si = hb->treelink[treemark][0];
	     // size of support is the same as the total number of scaling func and wavelets
	     // and the size of support is the same for all hb nodes in the cube
	     // only collecting neighboring finer level scaling functions
	     // printf("preparing temp...\n");///
	     // Assuming that all HB nodes under the same cube has the same neighbors
	     if(hb->linklength[treemark]>nnei){
	       // merge support for all neighbors of si
	       tlen = 0;	  
	    
	       MultilevelNeighbor(lj, si, nodes, nodeslength, levelcubes, subdividedlevelcubes);
	    
	       for(kk = 0; kk < *nodeslength; kk++){
		   int sk = nodes[kk];
	       
		   if (sk == -1 ) {
		     printf("Error ********************************************* \n",sk);
		     exit(-1);}
		 if (hb->support[sk][0] < 0)
		   // sk is in the finest cube,
		   // which should contain only 1 point
		   hb->mergeSupport(t,tlen,&sk,1);
		 else	
		   hb->mergeSupport(t,tlen,hb->support[sk],hb->blength[sk]);
		 // Alternatively, we could get the support from neighboring cubes
	       }
	       // loop through the support of si
	       for(jj = 0; jj < hb->blength[si]; jj++){
		 int sj = hb->support[si][jj];
		 FMMcount = 0;
		 // printf("hb->neighborHBLength[%d] = %d\n",sj,hb->neighborHBLength[sj]);///
		 // kk should go through the support of all neighbors of sj
		 for(kk = 0; kk < tlen; kk++){
		   // use FMM for missing entries in the finer level
		   // printf("kk %d \n", kk); 
		   int sk = t[kk];
		   // None is available
		   
		   if(hb->testAvailableHB(sj,sk, li + 1 + q*(wavMaxLevel + 1) - (q-1)*q/2)){
		     temp[jj][kk] = hb->highGet(hb2high[sj],hb2high[sk], li + 1 + q*(wavMaxLevel + 1) - (q-1)*q/2);
		     }
		   else{
		     FMMkk[FMMcount] = kk;
		     FMMsj[FMMcount] = sk;
		     FMMcount++;
		     // no highSet here, but will be set later after FMMevalHB
		   }
		 }
		 // Evaluate the interaction between sj and 
		 // printf("Ready for FMMevalHB\n");///
		 if(FMMcount>0)
		   // FMMevalHB should handle negative HB index
		   // as pointer to coordinates
		   // Here the indices are tranposed, but should still be correct assuming A is symmtric
		   // printf("Before FMM sj %d FMMcount %d \n", sj, FMMcount);
		   FMMevalHB(hb,FMMsj,FMMcount,sj,eval,FMMtemp,nei,dist,rout,site);
	         // printf("After FMM \n");
	         // printf("Returned from FMMevalHB\n");///
	         for(kk = 0; kk < FMMcount; kk++){				  
		   temp[jj][FMMkk[kk]] = eval[kk];
		 }
	       }
	       // compute scaling functions and wavelets in the cube
	       // printf("treemark: %d(%d)\n",treemark,hb->linklength[treemark]);///
	       for(jj = 0; jj < hb->linklength[treemark]; jj++){
		 int sj = hb->treelink[treemark][jj];
		 // printf("  Accessing hb->treelink[%d][%d] : %d -> %d\n",treemark,jj,sj,hb2high[sj]);
		 // Assuming all subnodes have the exactly same set of support
		 // Adding up for only the neighboring HB Nodes
		 for(ll = 0; ll < tlen; ll++){
		   val = 0;
		   int sl = t[ll];
		   for(kk = 0; kk < hb->blength[sj]; kk++){
		     val += hb->coef[sj][kk]*temp[kk][ll];
		     // printf("Something wrong %d %d %f %f %f \n",sj,kk,temp[kk][ll],hb->coef[sj][kk],val);
		   }
		   // printf("Cols: highset ...\n");
		   hb->highSet(hb2high[sj],hb2high[sl], li + q*(wavMaxLevel + 1) - (q-1)*q/2, val);
		   // printf("Cols: trackAvailableHB ...\n");///
		 }
	       }
	     }
	   }
	   treemark--;
	 }

	 delete[] levelcubes;
	 delete[] subdividedlevelcubes;
	 delete[] site;
	 delete[] rout;
	 delete[] dist;
}


void FWTd3::Kwblockmatrixcolmn(int lj, int li, int treemark, int nnei, double **temp, int *FMMkk, int *FMMsj, double *eval, int *FMMtemp, int *t, int *nodes, int *nodeslength, int wavMaxLevel){

         int tlen;
	 // int currenttreemark;
	 int FMMcount;
	 int jj;
	 int ll;
	 double val;
	 int kk;
	 int q = (wavMaxLevel - lj );
	 int qm = (wavMaxLevel - (lj + 1) );

	 // Initialize vectors for Matern kernel summation
	 int nsite = MAX_DIRECT_KERNEL;
	 double* dist   = new double[nsite];
	 double* site   = new double[nsite];                  // Array of interpolation sites
	 double* rout   = new double[nsite];          // spline evaluation results

	 //Declare cube placeholders to pass to distance criterion. Makes it faster.
	 int maxLengthcubes = 2 * nei->treesize;
	 int* levelcubes = new int[maxLengthcubes];
	 int* subdividedlevelcubes = new int[maxLengthcubes];	

	 int* hb2high = hb->hb2high;
	 
	 // loop through the cubes
	 printf("Columns Treemark %d lj %d\n",treemark, lj);

	 while (treemark > -1){
	   if (nei->tree[treemark][11] == li) {
	     int si = hb->treelink[treemark][0];
	     // size of support is the same as the total number of scaling func and wavelets
	     // and the size of support is the same for all hb nodes in the cube
	     // only collecting neighboring finer level scaling functions
	     // printf("preparing temp...\n");///
	     // Assuming that all HB nodes under the same cube has the same neighbors
	     if(hb->linklength[treemark]>nnei){
	       // merge support for all neighbors of si
	       tlen = 0;	       

	       MultilevelNeighbor(lj, si, nodes, nodeslength, levelcubes, subdividedlevelcubes);

	       for(kk = 0; kk < *nodeslength; kk++){
		   int sk = nodes[kk];
	       
		   if (sk == -1 ) {
		     printf("Error ********************************************* \n",sk);
		     exit(-1);}
		 if (hb->support[sk][0] < 0)
		   // sk is in the finest cube,
		   // which should contain only 1 point
		   hb->mergeSupport(t,tlen,&sk,1);
		 else	
		   hb->mergeSupport(t,tlen,hb->support[sk],hb->blength[sk]);
		 // Alternatively, we could get the support from neighboring cubes
	       }
	       // loop through the support of si
	       for(jj = 0; jj < hb->blength[si]; jj++){
		 int sj = hb->support[si][jj];
		 FMMcount = 0;
		 // printf("hb->neighborHBLength[%d] = %d\n",sj,hb->neighborHBLength[sj]);///
		 // kk should go through the support of all neighbors of sj
		 for(kk = 0; kk < tlen; kk++){
		   //use FMM for missing entries in the finer level
		   // printf("kk %d \n", kk); 
		   int sk = t[kk];
		   // None is available
		   if(hb->testAvailableHB(sj,sk, li + 1 + qm*(wavMaxLevel + 1) - (qm-1)*qm/2)){
		       temp[jj][kk] = hb->highGet(hb2high[sj],hb2high[sk], li + 1 + qm*(wavMaxLevel + 1) - (qm-1)*qm/2);
		   }
		   else{
		     //record the position to use FMMevalHB
		     FMMkk[FMMcount] = kk;
		     FMMsj[FMMcount] = sk;
		     FMMcount++;
		     //no highSet here, but will be set later after FMMevalHB
		   }
		 }
		 // Evaluate the interaction between sj and 
		 // printf("Ready for FMMevalHB\n");///
		 if(FMMcount>0)
		   // FMMevalHB should handle negative HB index
		   // as pointer to coordinates
		   // Here the indices are tranposed, but should still be correct assuming A is symmtric
		   // printf("Before FMM sj %d FMMcount %d \n", sj, FMMcount);
		   FMMevalHB(hb,FMMsj,FMMcount,sj,eval,FMMtemp,nei,dist,rout,site);
	         // printf("After FMM \n");
	         // printf("Returned from FMMevalHB\n");///
	         for(kk = 0; kk < FMMcount; kk++){				  
		   temp[jj][FMMkk[kk]] = eval[kk];
		 }
	       }
	       // compute scaling functions and wavelets in the cube
	       // printf("treemark: %d(%d)\n",treemark,hb->linklength[treemark]);///
	       for(jj = 0; jj < hb->linklength[treemark]; jj++){
		 int sj = hb->treelink[treemark][jj];
		 // printf("  Accessing hb->treelink[%d][%d] : %d -> %d\n",treemark,jj,sj,hb2high[sj]);
		 // Assuming all subnodes have the exactly same set of support
		 // Adding up for only the neighboring HB Nodes
		 for(ll = 0; ll < tlen; ll++){
		   val = 0;
		   int sl = t[ll];
		   for(kk = 0; kk < hb->blength[sj]; kk++){
		     val += hb->coef[sj][kk]*temp[kk][ll];
		   }
		   // printf("Cols: highset ...\n");
		   hb->highSet(hb2high[sj],hb2high[sl], li + q*(wavMaxLevel + 1) - (q-1)*q/2, val);
		   // printf("Cols: trackAvailableHB ...\n");///
		 }
	       }
	     }
	   }
	   treemark--;
	 }
	delete[] levelcubes;
	delete[] subdividedlevelcubes;
	delete[] site;
	delete[] rout;
	delete[] dist;
}

void FWTd3::Kwblockmatrixrow(int li, int lj, int currenttremark, int nnei, double **temp, int *FMMkk, int *FMMsj, double *eval, int *FMMtemp, int *t, int *nodes, int *nodeslength, int wavMaxLevel){
            
            // Process block Kw matrix row       
            int tlen;
	    int FMMcount;
	    int jj;
	    int ll;
	    double val;
	    int kk;
	    int treemark;
	    int qm = (wavMaxLevel - (lj + 1) );
	    int q  = (wavMaxLevel - lj);
	    int* hb2high = hb->hb2high;
	    // printf("Currenttremark %d li %d \n",currenttremark,li);
	    // printf("Currenttremark %d level %d \n",1,nei->tree[1][11]);

	    // Initialize vectors for Matern kernel summation
	    int nsite = MAX_DIRECT_KERNEL;
	    double* dist   = new double[nsite];
	    double* site   = new double[nsite];                  // Array of interpolation sites
	    double* rout   = new double[nsite];          // spline evaluation results

	    while (currenttremark > -1){

	      if (nei->tree[currenttremark][11] == lj) {

	          int si = hb->treelink[currenttremark][0];
			// size of support is the same as the total number of scaling func and wavelets
			// and the size of support is the same for all hb nodes in the cube
			//only collecting neighboring finer level scaling functions
		  if(hb->linklength[currenttremark]>nnei){
		    tlen = 0;
		    
		    int out = 0;
		    for(int i = 0; i < (hb->MutliLevelneighborHBNodesLengthTranspose[si]/2); i++){
		      if (hb->MutliLevelneighborHBNodesTranspose[si][2*i] == li) {
			treemark = hb->MutliLevelneighborHBNodesTranspose[si][2*i + 1];
			out = 1;
			break;
		      }
		    }
		    
		    if (out == 0) {
		      printf("Something went wrong \n");
		      exit(-1);
		    }

		    // printf("Column Current Treemark %d Level %d \n",treemark, nei->tree[treemark][11]);
		    // for(kk = 0; kk < *nodeslength; kk++){
		    //  int sk = nodes[kk];

		     for(kk = 0; kk < hb->neighborCubesLength[treemark]    ; kk++){
		      int sk = hb->neighborCubes[treemark][kk];
			// printf("currenttremark %d si %d Length %d kk % d sk %d \n",currenttremark, si, *nodeslength, kk, sk);
			if (hb->support[sk][0] < 0){
			// sk is in the finest cube,
			// which should contain only 1 point
			  hb->mergeSupport(t,tlen,&sk,1); }
			else
			    hb->mergeSupport(t,tlen,hb->support[sk],hb->blength[sk]);
			    // Alternatively, we could get the support from neighboring cubes
		     }	    
		    // loop through the support of si
		    // collect entries for the HB transform
		    for(jj = 0; jj < hb->blength[si]; jj++){
		       // sj loops through the support of si
		       int sj = hb->support[si][jj];
		       FMMcount = 0;
		       // printf("Rows: hb->neighborHBLength[%d] = %d\n",sj,hb->neighborHBLength[sj]);///
		       // loop through the support of neighbors of si
		       for(kk = 0; kk < tlen; kk++){
			  // all entries should be available now
			  int sk = t[kk];
			  // Check if pre computed
			 
			  // hb->sortAvailableHB(li + q*wavMaxLevel - (q-1)*q/2 );
			  if(hb->testAvailableHB(sk,sj, ( (li) + q*(wavMaxLevel + 1) - (q-1)*q/2 )  ))   { 
			    // exit(-1);
			      temp[jj][kk] = hb->highGet(hb2high[sk],hb2high[sj],( (li) + q*(wavMaxLevel + 1) - (q-1)*q/2) );
			  }
			  else {
			    printf("Something wrong testAvailableHB: (%d, %d). Increase ptilde or decrease tau. \n",sk,sj);
			    exit(-1);
			    temp[jj][kk] = 0;
			    
			  }
		       }
		       if(FMMcount>0)
			 FMMevalHB(hb,FMMsj,FMMcount,sj,eval,FMMtemp,nei,dist,rout,site);
		       // printf("Returned from FMMevalHB\n");///
		       for(kk = 0; kk < FMMcount; kk++){				  
			 temp[jj][FMMkk[kk]] = eval[kk];
		       }
		    }
		    //compute scaling functions and wavelets
		    for(jj = 0; jj < hb->linklength[currenttremark]; jj++){
		       // sj loops through the lifted HB nodes in the cube
		       int sj = hb->treelink[currenttremark][jj];
		       // Assuming all subnodes have the exactly same set of neighbors
		       // Adding up for only the neighboring HB Nodes
		       for(ll = 0; ll < tlen; ll++){
			  val = 0;
			  int sl = t[ll];
			  for(kk = 0; kk < hb->blength[sj]; kk++)
			    val += hb->coef[sj][kk]*temp[kk][ll];
			  
			  hb->highSet(hb2high[sl],hb2high[sj], (li + q*(wavMaxLevel + 1) - (q-1)*q/2 ),val);

			     if (li == lj){
			    }

		       }
		       // The neighbors of sj is the same as those of si
		       // The points in t and neighborHBNodes are for the same
		       // locations in high but on differnt level and in different orders
		    }
		  }		  
	      }
	      currenttremark--;
	    }

	delete[] site;
	delete[] rout;
	delete[] dist;
}


void FWTd3::Kwblockmatrixrowinit(int li, int lj, int currenttremark, int nnei, double **temp, int *FMMkk, int *FMMsj, double *eval, int *FMMtemp, int *t, int *nodes, int *nodeslength, int wavMaxLevel){
            
            // Process block Kw matrix row       
            int tlen;
	    int FMMcount;
	    int jj;
	    int ll;
	    double val;
	    int kk;
	    int treemark;

	    int* hb2high = hb->hb2high;

	    // Initialize vectors for Matern kernel summation
	    int nsite = MAX_DIRECT_KERNEL;
	    double* dist   = new double[nsite];
	    double* site   = new double[nsite];                  // Array of interpolation sites
	    double* rout   = new double[nsite];          // spline evaluation results

	    while (currenttremark > -1){

	      if (nei->tree[currenttremark][11] == lj) {

	          int si = hb->treelink[currenttremark][0];
			// size of support is the same as the total number of scaling func and wavelets
			// and the size of support is the same for all hb nodes in the cube
			// only collecting neighboring finer level scaling functions
		  if(hb->linklength[currenttremark]>nnei){
		    printf("Rows Currenttremark %d si %d li %d lj %d\n",currenttremark, si, li, lj);
		    tlen = 0;
		    
		    // Multilevel Distance criterion Upper levels
		    printf("Colmn si = %d \n", si);

		    for(int i = 0; i < hb->MutliLevelneighborHBNodesLengthTranspose[si]; i++){
		    }
		    
		    int out = 0;
		    for(int i = 0; i < (hb->MutliLevelneighborHBNodesLengthTranspose[si]/2); i++){
		      if (hb->MutliLevelneighborHBNodesTranspose[si][2*i] == li) {
			out = 1;
			treemark = hb->MutliLevelneighborHBNodesTranspose[si][2*i + 1];
			break;
		      }
		    }

		    if (out == 0) printf("Something wrong \n");

		       for(kk = 0; kk < hb->neighborCubesLength[treemark]    ; kk++){
			 int sk = hb->neighborCubes[treemark][kk];
			if (hb->support[sk][0] < 0){
			// sk is in the finest cube,
			// which should contain only 1 point
			  hb->mergeSupport(t,tlen,&sk,1); }
			else
			    hb->mergeSupport(t,tlen,hb->support[sk],hb->blength[sk]);
			    //Alternatively, we could get the support from neighboring cubes
		    }	    
		    //loop through the support of si
		    //collect entries for the HB transform
		    for(jj = 0; jj < hb->blength[si]; jj++){
		       //sj loops through the support of si
		       int sj = hb->support[si][jj];
		       FMMcount = 0;
		       // printf("Rows: hb->neighborHBLength[%d] = %d\n",sj,hb->neighborHBLength[sj]);///
		       // loop through the support of neighbors of si
		       for(kk = 0; kk < tlen; kk++){
			  //all entries should be available now
			  int sk = t[kk];
			  // Check if pre computed
			  if(hb->testAvailableHB(sk,sj,li+1) && li < -1){
			    temp[jj][kk] = hb->highGet(hb2high[sk],hb2high[sj],li+1);
			  }
			  else {
			    FMMkk[FMMcount] = kk;
			    FMMsj[FMMcount] = sk;
			    FMMcount++;
			  }
		       }
		       if(FMMcount>0)
			 FMMevalHB(hb,FMMsj,FMMcount,sj,eval,FMMtemp,nei,dist,rout,site);
		       // printf("Returned from FMMevalHB\n");///
		       for(kk = 0; kk < FMMcount; kk++){				  
			 temp[jj][FMMkk[kk]] = eval[kk];
		       }
		    }
		    // printf("FMM count %d \n",FMMcount);
		    //compute scaling functions and wavelets
		    // printf("Rows: currenttremark: %d(%d)\n",currenttremark,hb->linklength[currenttremark]);///
		    for(jj = 0; jj < hb->linklength[currenttremark]; jj++){
		       // sj loops through the lifted HB nodes in the cube
		       int sj = hb->treelink[currenttremark][jj];
		       // printf("  Accessing hb->treelink[%d][%d] : %d -> %d\n",currenttremark,jj,sj,hb2high[sj]);///
		       // Assuming all subnodes have the exactly same set of neighbors
		       // Adding up for only the neighboring HB Nodes
		       for(ll = 0; ll < tlen; ll++){
			  val = 0;
			  int sl = t[ll];
			  for(kk = 0; kk < hb->blength[sj]; kk++)
			    val += hb->coef[sj][kk]*temp[kk][ll];
			  
			  hb->highSet(hb2high[sl],hb2high[sj],li,val);
		       }
		       //The neighbors of sj is the same as those of si
		       //The points in t and neighborHBNodes are for the same
		       // locations in high but on differnt level and in different orders
		    }
		  }		  
		  // printf("Rows again: currenttremark: %d(%d)\n",currenttremark,hb->linklength[currenttremark]);///
	      }
	      currenttremark--;
	    }
	    delete[] site;
	delete[] rout;
	delete[] dist;
	// printf("Leaving Currenttremark %d \n", currenttremark);
}


void FWTd3::KwStore(cs* csKw, int li, int lj, int *nodes, int *nodeslength,int wavMaxLevel){

       // Store sparse matrix from block li and lj when li = lj
       // Only lower triangular entries are saved since the matrix is symmetrical
  
       //Declare cube placeholders to pass to distance criterion. Makes it faster.
       int maxLengthcubes = 2 * nei->treesize;
       int* levelcubes = new int[maxLengthcubes];
       int* subdividedlevelcubes = new int[maxLengthcubes];	

       int si,ii,jj;
       int q = (wavMaxLevel - lj);
       int* hb2high = hb->hb2high;
       printf("Store sparse matrix Initialize \n");
    


       int Kwrow_n = hb->wavlevels[li][1] - hb->wavlevels[li][0];
       printf("Kwrow_n %d Level i %d Level j %d \n", Kwrow_n, li, lj);

       for(ii = 0; ii < Kwrow_n; ii++){
	 si = hb->waveletNodes[hb->wavlist[ii + hb->wavlevels[li][0]]];
	 MultilevelNeighbor(lj, si, nodes, nodeslength, levelcubes, subdividedlevelcubes);
	 for(jj = 0; jj < *nodeslength; jj++){
	   int sj = nodes[jj];	   
	   if(hb->wavelet[sj]){
	     // Save lower triangular only
	     if(hb->waveletIndex[sj] <=  ii + hb->wavlevels[li][0]){
	     cs_entry(csKw, ii + hb->wavlevels[li][0], hb->waveletIndex[sj], hb->highGet(hb2high[si],
			   hb2high[sj], li + q*(wavMaxLevel + 1) - (q - 1)*q/2) ); 
	     }
	   }
	 }
       }

       delete[] levelcubes;
       delete[] subdividedlevelcubes;
}

void FWTd3::KwStoreInt(cs* csKw, int li, int lj, int *nodes, int *nodeslength,int wavMaxLevel){

       // Store sparse matrix from block li and lj
  
       // Declare cube placeholders to pass to distance criterion. Makes it faster.
       int maxLengthcubes = 2 * nei->treesize;
       int* levelcubes = new int[maxLengthcubes];
       int* subdividedlevelcubes = new int[maxLengthcubes];	

       int si,sj,ii,jj;
       int qm = (wavMaxLevel - (lj + 1) );
       int q = (wavMaxLevel - lj);
       int* hb2high = hb->hb2high;
       printf("Store sparse matrix Interior\n");
    
       int Kwrow_n = hb->wavlevels[li][1] - hb->wavlevels[li][0];
       printf("Kwrow_n %d Level i %d Level j %d \n", Kwrow_n, li, lj);


       for(ii = 0; ii < Kwrow_n; ii++){
	 si = hb->waveletNodes[hb->wavlist[ii + hb->wavlevels[li][0]]];

	 MultilevelNeighbor(lj, si, nodes, nodeslength, levelcubes, subdividedlevelcubes);

	 for(jj = 0; jj < *nodeslength; jj++){
	   sj = nodes[jj];

	   if(hb->wavelet[sj]){
	     // Save lower triangular only
	     if(hb->waveletIndex[sj] <=  ii + hb->wavlevels[li][0]){
	     cs_entry(csKw, ii + hb->wavlevels[li][0], hb->waveletIndex[sj], hb->highGet(hb2high[si],
		      hb2high[sj], li + q*(wavMaxLevel + 1) - (q - 1)*q/2) );
	     }
	   }
	 }
       }
		 
       delete[] levelcubes;
       delete[] subdividedlevelcubes;
}

	       
void FWTd3::MultilevelNeighbor(int lj, int si, int *nodes, int *nodeslength, int *levelcubes, int *subdividedlevelcubes){

	// Building Distance criterion for j > i levels
	// ------------------------------------------
	int  counter    = 0;
	int  precounter = 0;
        int totalcounter = 0;
        int multilevel;
	int cube;
	int level;
	int treemark;
	int maxLength = 2 * nei->treesize;

	// Leave here for debugging purposes. Then move to be declared out and pass the pointers
        // printf("Compute j > i distance criterion neighborhoods \n");       	
        cube = hb->cubenode[si];
	int silevel = nei->tree[cube][11];
        precounter = 0;
	   	 
	// Collect starting level neighborhoods cubes
	// printf("Current level cubes \n");
	for(int j = 0; j < hb->neighborHBCubesLength[si]; j++){
	    treemark = hb->neighborHBCubes[si][j];
	    level = nei->tree[treemark][11];
	    levelcubes[j] = treemark;
	    // printf("%d \n", levelcubes[j]);
	    precounter++;
	}
	// printf("\n");

        multilevel = level;
        while (lj > multilevel){ 	      
	// Subdivide current level and collect children cubes
	//  printf("Subdivided cubes precounter %d \n", precounter);
	  counter = 0;
	  for(int q = 0; q < precounter; q++){
	    for(int p = 2; p < 10; p++){
	      if(nei->tree[levelcubes[q]][p] > 0){
		subdividedlevelcubes[counter] = nei->tree[levelcubes[q]][p];
		counter++;
	      }
	    }
	  }
	  // printf("\n");
	  precounter = counter;	  
	  for(int q = 0; q < precounter; q++)
	    {
	      levelcubes[q]=subdividedlevelcubes[q];	      
	    }
	  multilevel++;
	}
	// Store nodes
	    for(int q = 0; q < precounter; q++){
	      for(int f = 0; f < hb->linklength[levelcubes[q]]; f++){
		// Only store neighbors functions at the multilevel level.
		if (nei->tree[hb->cubenode[hb->treelink[levelcubes[q]][f]]][11] == multilevel) {
		  nodes[totalcounter] = hb->treelink[levelcubes[q]][f];
		  totalcounter++;
		}
	      }
	    }
	    *nodeslength = totalcounter;
}
