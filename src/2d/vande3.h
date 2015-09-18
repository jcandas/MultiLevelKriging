#ifndef VANDE3_H
#define VANDE3_H
#include <vector>
using namespace std;
class Vande3{
	static vector<double> cx,cy,cz;
	static vector<double> polyx,polyy,polyz;
	static double temp_a,temp_b;
	public:
	static double** vandemonde3(double** coord, int van, int n);
	static void get_vandemonde3(double** coord, int van, int n, double** y);
	static void chebpoly(int n, vector<double>& x, vector<double>& c);
	static void chebpoly_1(int n, vector<double>& x, vector<double>& c);
};
	      
#endif
