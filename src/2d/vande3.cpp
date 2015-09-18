#include <math.h>
#include <stdio.h>
#include "vande3.h"
using namespace std;
//The input coordinates are expected to be within [-10,10]
double** Vande3::vandemonde3(double** coord, int van, int n){

		int ii,jj,kk,ll;
//		printf("Entering vande3 ... (van=%d,n=%d)\n",van,n);///
		double m1 = 100;
		double m2 = 100;
		double m3 = 100;

		cx.clear();
		cy.clear();
		cz.clear();

		for (ii = 0; ii < n; ii++){
			cx.push_back(coord[ii][0]/m1);
			cy.push_back(coord[ii][1]/m2);
			cz.push_back(coord[ii][2]/m3);
		}

		int m = van*(van+1)/2; //2D

		double** y = new double*[m];

		for (ii = 0; ii < m; ii++) {
			y[ii] = new double[n];
		}
				
		int cc=0;
		polyx.clear();
		polyy.clear();
		polyz.clear();
		for (kk = 0; kk < van; kk++){
		        chebpoly(kk,cx,polyx);
			for (jj = 0; jj < van-kk; jj++){
				chebpoly(jj,cy,polyy);
				// for (ii = 0; ii < van-kk-jj; ii++){
				//	chebpoly(ii,cz,polyz);
					for (ll = 0; ll < n; ll++)
						y[cc][ll] = polyx[ll]*polyy[ll];
					cc++;
			        //	polyz.clear();
				// }
				polyy.clear();
			}
			polyx.clear();
		}
		//assert(cc==m);
		cx.clear();
		cy.clear();
		cz.clear();
		// printf("cc = %d \n",cc);///
		// printf("vande3: Done\n");///
		return y;
		
	}
	
void Vande3::get_vandemonde3(double** coord, int van, int n,double** y){
		int ii,jj,kk,ll;
		// printf("Entering vande3 ... (van=%d,n=%d)\n",van,n);///
		double m1 = 100;
		double m2 = 100;
		double m3 = 100;

		cx.clear();
		cy.clear();
		cz.clear();
		for (ii = 0; ii < n; ii++){
			cx.push_back(coord[ii][0]/m1);
			cy.push_back(coord[ii][1]/m2);
			cz.push_back(coord[ii][2]/m3);
		}

		int m = van*(van+1)/2; //2D
		
		int cc=0;
		polyx.clear();
		polyy.clear();
		polyz.clear();
		for (kk = 0; kk < van; kk++){
		  //			printf("vande3 : loop kk = %d\n",kk);///
		        chebpoly_1(kk,cx,polyx);
			for (jj = 0; jj < van-kk; jj++){
			  //				printf("vande3 : loop jj = %d\n",jj);///
				chebpoly_1(jj,cy,polyy);

					for (ll = 0; ll < n; ll++)
						y[cc][ll] = polyx[ll]*polyy[ll];
					cc++;
					// polyz.clear();
					//}
				polyy.clear();
			}
			polyx.clear();
		}
		// assert(cc==m);
		cx.clear();
		cy.clear();
		cz.clear();
	}

void Vande3::chebpoly(int n, vector<double>& x, vector<double>& c){
		int ii,jj,kk;
		c.clear();
		if (n == 0)
			for (ii = 0; ii < x.size(); ii++)
				c.push_back(1);
		else if (n == 1)
			for (ii = 0; ii < x.size(); ii++)
				c.push_back(x[ii]);
		else{
		  // printf("chebypoly: a = new double[%d]\n",x.size());///
			double* a = new double[x.size()];
			// printf("chebypoly: b = new double[%d]\n",x.size());///
			double* b = new double[x.size()];
			for (ii = 0; ii < x.size(); ii++)
				b[ii] = x[ii];
			for (ii = 0; ii < x.size(); ii++)
				c.push_back(0);
			for (kk = 1; kk < n; kk++)
				for (jj = 0; jj < x.size(); jj++){
					c[jj] = 2 * x[jj] * b[jj] - a[jj];
					a[jj] = b[jj];
					b[jj] = c[jj];
				}
			delete[] a;
			delete[] b;
		}	
	}
	
void Vande3::chebpoly_1(int n, vector<double>& x, vector<double>& c){
		int ii,jj,kk;
		c.clear();
		if (n == 0)
			for (ii = 0; ii < x.size(); ii++)
				c.push_back(1);
		else if (n == 1)
			for (ii = 0; ii < x.size(); ii++)
				c.push_back(x[ii]);
		else{
		  // printf("chebypoly: a = new double[%d]\n",x.size());///
			temp_a = 1;
			// printf("chebypoly: b = new double[%d]\n",x.size());///
			temp_b = x[0];
			c.push_back(0);
			for (kk = 1; kk < n; kk++){
				c[0] = 2 * x[0] * temp_b - temp_a;
				temp_a = temp_b;
				temp_b = c[0];
				}
		}	
	}
	

