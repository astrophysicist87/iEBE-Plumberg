#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include <vector>
#include <stdlib.h>
#include <cstdarg>

#include<gsl/gsl_sf_bessel.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#include "Arsenal.h"

using namespace std;

#define TINY 1.0e-25

unsigned long int random_seed()
{
  unsigned long int seed = 0;
  unsigned long int *seed_ptr = &seed;

  ifstream dev_urandom ("/dev/urandom", ios::in | ios::binary);
  
  dev_urandom.read((char *) seed_ptr, sizeof (long int));

  dev_urandom.close();
  return(*seed_ptr);
}

int sgn(double val)
{
    return (0 < val) - (val < 0);
}

void logspace(double * x, double a, double b, int n)
{
//returns vector x of n logarithmically spaced points, with a and b as endpoints
	double b_by_a = b/a;
	double exponent = 1./(double)(n-1);
	double gamma = pow(b_by_a, exponent);

	//assume x has length n already
	x[0] = a;
	for (int i = 1; i < n; i++)
		x[i] = gamma*x[i-1];

	return;
}

void linspace(double * x, double a, double b, int n)
{
//returns vector x of n linearly spaced points, with a and b as endpoints
	double Del_x = (b-a)/(double)(n-1);

	//assume x has length n already
	for (int i = 0; i < n; i++)
		x[i] = a + Del_x * (double)i;

	return;
}

void linspace(vector<double> & x, double a, double b)
{
//returns vector x of n linearly spaced points, with a and b as endpoints
	int n = x.size();
	double Del_x = (b-a)/(double)(n-1);

	//assume x has length n already
	for (int i = 0; i < n; i++)
		x[i] = a + Del_x * (double)i;

	return;
}


void stratify_npts(double a, double b, int n1, int npts, double * x)
{
	// take a and b to be positive
	// region 1 is "inner" stratum, region 2 is "outer" stratum
	// a - scale of inner stratum
	// b - scale of outer stratum
	if (n1%2==0) exit(1);
	if (a < 0.0) a *= -1.;
	if (b < 0.0) b *= -1.;
	int m = npts - n1;
	
	double Delta_x_1 = 2.*a / double(n1-1);
	double Delta_x_2 = 2.*(b-a) / (double)m;
	
	for (int i = 0; i < (m/2); i++)
	{
		x[i] = -b + (double)i * Delta_x_2;
		x[npts - i - 1] = b - (double)i * Delta_x_2;
	}
	for (int i = 0; i < n1; i++)
		x[i+(m/2)] = -a + (double)i * Delta_x_1;

	return;
}

void scalepoints(double * x, double a, double b, double scale, int n)
{
// n is length of x
// returns x with x[0] == a, x[n-1] == b
// for even n: 50% of points (linearly spaced) above scale, 50% (linearly spaced) below
// for odd n: center point located at scale, 50% of rest (linearly spaced) above and 50% (linearly spaced) below

	double * dummy;
	if (n%2 == 0)	// if n is even
	{
		dummy = new double [n/2];
		linspace(dummy, a, scale, n/2);
		for (int i=0; i<(n/2); i++)
			x[i] = dummy[i];
		linspace(dummy, scale, b, n/2);
		for (int i=(n/2); i<n; i++)
			x[i] = dummy[i-(n/2)];
	}
	else	// if n is odd
	{
		dummy = new double [(n+1)/2];
		linspace(dummy, a, scale, (n+1)/2);
		for (int i=0; i<((n+1)/2); i++)
			x[i] = dummy[i];
		linspace(dummy, scale, b, (n+1)/2);
		if (abs(dummy[0] - x[(n-1)/2]) >= 1.e-10)
		{
			cerr << "Scalepoints: mismatch of indices" << endl
				<< "dummy[0] = " << dummy[0] << ", x[(n+1)/2] = " << x[(n+1)/2] << endl;
			exit(1);
		}
		for (int i=((n+1)/2); i<n; i++)
			x[i] = dummy[i - ((n-1)/2)];
	}

	return;
}

static int wt[16][16] = {
{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
{-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0},
{2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0},
{0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
{0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1},
{0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1},
{-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0},
{9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2},
{-6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2},
{2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0},
{-6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1},
{4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1}
};


//**********************************************************************
long binarySearch(double * A, int length, double value, bool skip_out_of_range /*== true*/, bool verbose /*== false*/)
// Return the index of the largest number less than value in the list A
// using binary search. Index starts with 0.
// If skip_out_of_range is set to true, then it will return -1 for those
// samples that are out of the table range (default is true).
{
   //int length = A->size();
   int idx_i, idx_f, idx;
   idx_i = 0;
   idx_f = length-1;

   if(value > A[idx_f])
   {
      if (verbose) cerr << "binarySearch: desired value is too large, exceeding the end of the table: value = " << value << " and A[idx_f] = " << A[idx_f] << endl;
      if (skip_out_of_range) return -1;
      exit(1);
   }
   if(value < A[idx_i])
   {
      if (verbose) cerr << "binarySearch: desired value is too small, exceeding the beginning of table: value = " << value << " and A[idx_i] = " << A[idx_i] << endl;
      if (skip_out_of_range) return -1;
      exit(1);
   }
   idx = (int) floor((idx_f+idx_i)/2.);
   //if (verbose) cerr << "Start: idx = " << idx << endl;
   while((idx_f-idx_i) > 1)
   {
     if(A[idx] < value)
     {
        idx_i = idx;
		//if (verbose) cerr << "idx_i = " << idx_i << endl;
     }
     else
     {
        idx_f = idx;
		//if (verbose) cerr << "idx_f = " << idx_f << endl;
     }
     idx = (int) floor((idx_f+idx_i)/2.);
	//if (verbose) cerr << "End: idx = " << idx << endl;
   }
   return(idx_i);
}


void get_1D_derivatives(double * x, double * f, double * derivatives, int length, double default_edge_fill = 0.0)
{
	for (int i = 1; i <= length-2; i++)
		derivatives[i] = (f[i+1] - f[i-1]) / (x[i+1] - x[i-1]);	//use centered differencing
	derivatives[0] = default_edge_fill;
	derivatives[length-1] = default_edge_fill;

	return;
}

void get_2D_derivatives(double * x, double * y, double ** f, double ** f1, double ** f2, double ** f12, int x_length, int y_length, double default_edge_fill = 0.0)
{
	//for the points not on the edge, get derivatives from surrounding points
	for (int ix = 1; ix <= x_length-2; ix++)
	for (int iy = 1; iy <= y_length-2; iy++)
	{//use centered differencing
		f1[ix][iy] = (f[ix+1][iy] - f[ix-1][iy]) / (x[ix+1] - x[ix-1]);
		f2[ix][iy] = (f[ix][iy+1] - f[ix][iy-1]) / (y[iy+1] - y[iy-1]);
		f12[ix][iy] = (f[ix+1][iy+1] - f[ix+1][iy-1] - f[ix-1][iy+1] + f[ix-1][iy-1]) / ((x[ix+1] - x[ix-1]) * (y[iy+1] - y[iy-1]));
	}
	for (int ix = 0; ix <= x_length-1; ix++)
	{
		f1[ix][0] = default_edge_fill;
		f2[ix][0] = default_edge_fill;
		f12[ix][0] = default_edge_fill;
		f1[ix][y_length-1] = default_edge_fill;
		f2[ix][y_length-1] = default_edge_fill;
		f12[ix][y_length-1] = default_edge_fill;
	}
	for (int iy = 0; iy <= y_length-1; iy++)
	{
		f1[0][iy] = default_edge_fill;
		f2[0][iy] = default_edge_fill;
		f12[0][iy] = default_edge_fill;
		f1[x_length-1][iy] = default_edge_fill;
		f2[x_length-1][iy] = default_edge_fill;
		f12[x_length-1][iy] = default_edge_fill;
	}

	return;
}

void bcucof(double * y, double * y1, double * y2, double * y12, double d1, double d2, double ** c)
{
//Given arrays y[0..3], y1[0..3], y2[0..3], and y12[0..3], containing the function, gradients,
//and cross-derivative at the four grid points of a rectangular grid cell (numbered counterclockwise
//from the lower left), and given d1 and d2, the length of the grid cell in the 1 and 2 directions, this
//routine returns the table c[0..3][0..3] that is used by routine bcuint for bicubic interpolation.
	int l,k,j,i;
	double xx, d1d2=d1*d2;
	double * cl = new double [16];
	double * x = new double [16];
	for (i=0;i<4;i++)
	{
		x[i]=y[i];
		x[i+4]=y1[i]*d1;
		x[i+8]=y2[i]*d2;
		x[i+12]=y12[i]*d1d2;
	}
	for (i=0;i<16;i++)
	{
		xx=0.0;
		for (k=0;k<16;k++) xx += wt[i][k]*x[k];
		cl[i]=xx;
	}
	l=0;
	for (i=0;i<4;i++)
	for (j=0;j<4;j++)
		c[i][j]=cl[l++];

	return;
}

/////////////////////////////////////////////////////////////////////////////////////////
void bcuint(double * y, double * y1, double * y2, double * y12,
double x1l, double x1u, double x2l, double x2u,
double x1, double x2, double &ansy, double &ansy1, double &ansy2)
{
//Bicubic interpolation within a grid square. Input quantities are y,y1,y2,y12 (as described in
//bcucof); x1l and x1u, the lower and upper coordinates of the grid square in the 1 direction;
//x2l and x2u likewise for the 2 direction; and x1,x2, the coordinates of the desired point for
//the interpolation. The interpolated function value is returned as ansy, and the interpolated
//gradient values as ansy1 and ansy2. This routine calls bcucof.
	int i;
	double t,u, d1=x1u-x1l, d2=x2u-x2l;
	//MatDoub c(4,4);
	double ** c = new double * [4];
	for (int ic = 0; ic<4; ic++)
		c[ic] = new double [4];
	bcucof(y,y1,y2,y12,d1,d2,c);
	//Get the c’s.
	if (x1u == x1l || x2u == x2l)
		throw("Bad input in routine bcuint");
	t=(x1-x1l)/d1;
	//Equation (3.6.4).
	u=(x2-x2l)/d2;
	ansy=0.0;
	ansy2=0.0;
	ansy1=0.0;
	for (i=3;i>=0;i--)
	{
		//Equation (3.6.6).
		ansy=t*ansy+((c[i][3]*u+c[i][2])*u+c[i][1])*u+c[i][0];
		ansy2=t*ansy2+(3.0*c[i][3]*u+2.0*c[i][2])*u+c[i][1];
		ansy1=u*ansy1+(3.0*c[3][i]*t+2.0*c[2][i])*t+c[1][i];
	}
	//for (int ic = 0; ic<4; ic++)
	//for (int icp = 0; icp<4; icp++)
	//	cout << "DEBUG: " << ic << "   " << icp << "   " << c[ic][icp] << endl;
	ansy1 /= d1;
	ansy2 /= d2;

	return;
}

/////////////////////////////////////////////////////////////////////////////////////////

void ratint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	int m,i,ns=1;
	double w,t,hh,h,dd,*c,*d;
	c = new double [n];
	d = new double [n];
	hh=fabs(x-xa[0]);
	for (i=0;i<n;i++)
	{
		h=fabs(x-xa[i]);
		if (h == 0.0)
		{
			*y=ya[i];
			*dy=0.0;
			return;
		}
		else if (h < hh)
		{
			ns=i+1;
			hh=h;
		}
		c[i]=ya[i];
		d[i]=ya[i]+TINY;
	}
	*y=ya[(ns--) - 1];
	for (m=1;m<n;m++)
	{
		for (i=1;i<=n-m;i++)
		{
			w=c[i]-d[i-1];
			h=xa[i+m-1]-x;
			t=(xa[i-1]-x)*d[i-1]/h;
			dd=t-c[i];
			if (dd == 0.0) cerr << "Error in routine ratint" << endl;
			dd=w/dd;
			d[i-1]=c[i]*dd;
			c[i-1]=t*dd;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns] : d[(ns--) - 1]));
	}
	
	delete [] c;
	delete [] d;

	return;
}


/////////////////////////////////////////////////////////////////////////////////////////

double ratint(vector<double> & xa, vector<double> & ya, double x, double *dy)
{
	double y = 0.0;
	int m,i,ns=1;
	const int n = xa.size();
	double w,t,hh,h,dd;
	//double *c,*d;
	//c = new double [n];
	//d = new double [n];
	double c[n];
	double d[n];
	hh=fabs(x-xa[0]);
	for (i=0;i<n;i++)
	{
		h=fabs(x-xa[i]);
		if (h == 0.0)
		{
//			*y=ya[i];
			y=ya[i];
			*dy=0.0;
			return(y);
		}
		else if (h < hh)
		{
			ns=i+1;
			hh=h;
		}
		c[i]=ya[i];
		d[i]=ya[i]+TINY;
	}
//	*y=ya[(ns--) - 1];
	y=ya[(ns--) - 1];
	for (m=1;m<n;m++)
	{
		for (i=1;i<=n-m;i++)
		{
			w=c[i]-d[i-1];
			h=xa[i+m-1]-x;
			t=(xa[i-1]-x)*d[i-1]/h;
			dd=t-c[i];
			if (dd == 0.0) cerr << "Error in routine ratint" << endl;
			dd=w/dd;
			d[i-1]=c[i]*dd;
			c[i-1]=t*dd;
		}
//		*y += (*dy=(2*ns < (n-m) ? c[ns] : d[(ns--) - 1]));
		y += (*dy=(2*ns < (n-m) ? c[ns] : d[(ns--) - 1]));
	}
	
	//delete [] c;
	//delete [] d;

	return(y);
}


/////////////////////////////////////////////////////////////////////////////////////////


void polint(double xa[], double ya[], long n, double x, double *y, double *dy)
//Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns a value y, and
//an error estimate dy. If P(x) is the polynomial of degree n − 1 such that P(xai) = yai
//, i = 1,...,n, then the returned value y = P(x).
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double * c = new double [n];
	double * d = new double [n];
	dif=fabs(x-xa[0]);
	//c=vector(1,n);
	//d=vector(1,n);
	for (i=1;i<=n;i++)
	{ //Here we ﬁnd the index ns of the closest table entry,
		if ( (dift=fabs(x-xa[i-1])) < dif)
		{
			ns=i;
			dif=dift;
		}
		c[i-1]=ya[i-1]; //and initialize the tableau of c’s and d’s.
		d[i-1]=ya[i-1];
	}
	*y=ya[(ns--)-1]; //This is the initial approximation to y.
	for (m=1;m<n;m++)
	{ //For each column of the tableau,
		for (i=1;i<=n-m;i++)
		{ //we loop over the current c’s and d’s and update them.
			ho=xa[i-1]-x;
			hp=xa[i+m-1]-x;
			w=c[i]-d[i-1];
			if ( (den=ho-hp) == 0.0) cerr << "Error in routine polint" << endl;
			//This error can occur only if two input xa’s are (to within roundoff) identical.
			den=w/den;
			d[i-1]=hp*den; //Here the c’s and d’s are updated.
			c[i-1]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1-1] : d[(ns--)-1]));
		//After each column in the tableau is completed, we decide which correction, c or d,
		//we want to add to our accumulating value of y, i.e., which path to take through the
		//tableau—forking up or down. We do this in such a way as to take the most “straight
		//line” route through the tableau to its apex, updating ns accordingly to keep track of
		//where we are. This route keeps the partial approximations centered (insofar as possible)
		//on the target x. The last dy added is thus the error indication.
	}
}


/////////////////////////////////////////////////////////////////////////////////////////

void polin2(double * x1a, double * x2a, double ** ya, long m, long n, double x1, double x2, double *y, double *dy)
//Given arrays x1a[1..m] and x2a[1..n] of independent variables, and a submatrix of function
//values ya[1..m][1..n], tabulated at the grid points deﬁned by x1a and x2a; and given values
//x1 and x2 of the independent variables; this routine returns an interpolated function value y,
//and an accuracy indication dy (based only on the interpolation in the x1 direction, however).
{
bool do_x2_first = true;
	if (do_x2_first)
	{
		int j;
		double * ymtmp = new double [m];
		//ymtmp=vector(1,m);
		for (j=0;j<m;j++)	//Loop over rows.
		{
			//cout << "Made it to loop #" << j << endl;
			polint(x2a, ya[j], n, x2, &ymtmp[j], dy); //Interpolate answer into temporary storage.
		}
		//cout << "Made it through!" << endl;
		polint(x1a, ymtmp, m, x1, y, dy); //Do the ﬁnal interpolation.
		//cout << "Finished everything!" << endl;
	}
	else
	{
		//define and set transpose of data values
		double ** yaT = new double * [n];
		for (int ix2 = 0; ix2 < n; ix2++)
		{
			yaT[ix2] = new double [m];
			for (int ix1 = 0; ix1 < m; ix1++)
				yaT[ix2][ix1] = ya[ix1][ix2];
		}

		int j;
		double * ymtmp = new double [n];
		for (j=0;j<n;j++)
			polint(x1a, yaT[j], m, x1, &ymtmp[j], dy);
		polint(x2a, ymtmp, n, x2, y, dy);
	}
}


/////////////////////////////////////////////////////////////////////////////////////////


//**********************************************************************
double interpLinearDirect(double * x, double * y, double x0, long size, bool returnflag /*= false*/, double default_return_value /* = 0*/)
// Returns the interpreted value of y=y(x) at x=x0 using linear interpolation method.
// -- x,y: the independent and dependent tables; x is assumed to be equal spaced and increasing
// -- x0: where the interpolation should be performed
{
	//long size = x->size();
	if (size==1) {cout<<"interpLinearDirect warning: table size = 1"<<endl; return y[0];}
	double dx = x[1]-x[0]; // increment in x

	// if close to left end:
	if (abs(x0-x[0])<dx*1e-30) return y[0];

	// find x's integer index
	long idx = floor((x0-x[0])/dx);

	if (idx<0 || idx>=size-1)
	{
		if (!returnflag)	//i.e., if returnflag is false, exit
		{
			cout    << "interpLinearDirect: x0 out of bounds." << endl
				<< "x ranges from " << x[0] << " to " << x[size-1] << ", "
				<< "x0=" << x0 << ", " << "dx=" << dx << ", " << "idx=" << idx << endl;
			exit(1);
		}
		else
		{
			idx = (idx<0) ? 0 : size-2;	//uses extrapolation
		}
		//else return default_return_value;
	}

  return y[idx] + (y[idx+1]-y[idx])/dx*(x0-x[idx]);
}

//**********************************************************************
double interpLinearNondirect(double * x, double * y, double x0, long size, bool returnflag /*= false*/, double default_return_value /* = 0*/)
// Returns the interpreted value of y=y(x) at x=x0 using linear interpolation method.
// -- x,y: the independent and dependent tables; x is assumed to be increasing but not equal spaced
// -- x0: where the interpolation should be performed
{
	//long size = x->size();
	if (size==1) {cout<<"interpLinearNondirect warning: table size = 1"<<endl; return y[0];}
	double dx = x[1]-x[0]; // increment in x

	// if close to left end:
	if (abs(x0-x[0])<dx*1e-30) return y[0];

	// find x's integer index
	//long idx = floor((x0-x[0])/dx);
	long idx = binarySearch(x, size, x0, true);
	if (idx<0 || idx>=size-1)
	{
cout << "returnflag = " << returnflag << endl;
//debugger(__LINE__, __FILE__);
		if (!returnflag)	//i.e., if returnflag is false, exit
		{
			cout    << "interpLinearNondirect: x0 out of bounds." << endl
				<< "x ranges from " << x[0] << " to " << x[size-1] << ", "
				<< "x0=" << x0 << ", " << "dx=" << dx << ", " << "idx=" << idx << endl;
			exit(1);
		}
		else
		{
			idx = (idx<0) ? 0 : size-2;	//uses extrapolation
		}
		//else return default_return_value;
	}

	return y[idx] + (y[idx+1]-y[idx])/(x[idx+1]-x[idx])*(x0-x[idx]);
}


//**********************************************************************
double interpBiLinearDirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size, bool returnflag /*= false*/, double default_return_value /* = 0*/)
{
	//long size = x->size();
	if (x_size==1 && y_size==1) {cout<<"interpBiLinearDirect warning: table size = 1"<<endl; return z[0][0];}
	double dx = x[1]-x[0]; // increment in x
	//double dy = y[1]-y[0]; // increment in y

	// assume not close to edges for now...
	// find x's integer index
	long xidx = floor((x0-x[0])/dx);
	if (xidx < 0 || xidx >= x_size-1)
	{
		if (!returnflag)	//i.e., if returnflag is false, exit
		{
			cerr << "interpBiLinearDirect(): index out of range!  Aborting!" << endl
				<< "interpBiLinearDirect(): x_size = " << x_size << ", x0 = " << x0 << ", xidx = " << xidx << endl;
			exit(1);
		}
		else return (default_return_value);
	}

	double xidxINT = interpLinearDirect(y, z[xidx], y0, y_size, returnflag);
	double xidxp1INT = interpLinearDirect(y, z[xidx+1], y0, y_size, returnflag);

	return xidxINT + (xidxp1INT-xidxINT)/dx*(x0-x[xidx]);
}


//**********************************************************************
double interpBiLinearNondirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size, bool returnflag /*= false*/, double default_return_value /* = 0*/)
{
	//long size = x->size();
	if (x_size==1 && y_size==1) {cout<<"interpBiLinearNondirect warning: table size = 1"<<endl; return z[0][0];}
	double dx = x[1]-x[0]; // increment in x
	//double dy = y[1]-y[0]; // increment in y

	// assume not close to edges for now...
	// find x's integer index
	//long xidx = floor((x0-x[0])/dx);
	long xidx = binarySearch(x, x_size, x0, true);
	if (xidx < 0 || xidx >= x_size-1)
	{
		if (!returnflag)	//i.e., if returnflag is false, exit
		{
			cerr << "interpBiLinearNondirect(): index out of range!  Aborting!" << endl
				<< "interpBiLinearNonDirect(): x_size = " << x_size << ", x0 = " << x0 << ", xidx = " << xidx << endl;
			exit(1);
		}
		else return (default_return_value);
	}

	double xidxINT = interpLinearNondirect(y, z[xidx], y0, y_size, returnflag);
	double xidxp1INT = interpLinearNondirect(y, z[xidx+1], y0, y_size, returnflag);

	//return xidxINT + (xidxp1INT-xidxINT)/dx*(x0-x[xidx]);
	return xidxINT + (xidxp1INT-xidxINT)/(x[xidx+1]-x[xidx])*(x0-x[xidx]);
}


//**********************************************************************
double interpTriLinearDirect(double * x, double * y, double * z, double *** f, double x0, double y0, double z0, long x_size, long y_size, long z_size, bool returnflag /*= false*/, double default_return_value /* = 0*/)
{
	//long size = x->size();
	if (x_size==1 && y_size==1) {cout<<"interpTriLinearDirect warning: table size = 1"<<endl; return f[0][0][0];}
	double dx = x[1]-x[0]; // increment in x

	// assume not close to edges for now...
	// find x's integer index
	long xidx = floor((x0-x[0])/dx);

	if (xidx < 0 || xidx >= x_size-1)
	{
		if (!returnflag)	//i.e., if returnflag is false, exit
		{
			cerr << "interpTriLinearDirect(): index out of range!  Aborting!" << endl
				<< "interpTriLinearDirect(): x_size = " << x_size << ", x0 = " << x0 << ", xidx = " << xidx << endl;
			exit(1);
		}
		else return (default_return_value);
	}

	double xidxINT = interpBiLinearDirect(y, z, f[xidx], y0, z0, y_size, z_size, returnflag);
	double xidxp1INT = interpBiLinearDirect(y, z, f[xidx+1], y0, z0, y_size, z_size, returnflag);

	return xidxINT + (xidxp1INT-xidxINT)/dx*(x0-x[xidx]);
}


//**********************************************************************
double interpTriLinearNondirect(double * x, double * y, double * z, double *** f, double x0, double y0, double z0, long x_size, long y_size, long z_size, bool returnflag /*= false*/, double default_return_value /* = 0*/)
{
	//long size = x->size();
	if (x_size==1 && y_size==1 && z_size==1) {cout<<"interpTriLinearNondirect warning: table size = 1"<<endl; return f[0][0][0];}

	// assume not close to edges for now...
	// find x's integer index
	//long xidx = floor((x0-x[0])/dx);
	long xidx = binarySearch(x, x_size, x0, true);

	if (xidx < 0 || xidx >= x_size-1)
	{
		if (!returnflag)	//i.e., if returnflag is false, exit
		{
			cerr << "interpTriLinearNonDirect(): index out of range!  Aborting!" << endl
				<< "interpTriLinearNonDirect(): x_size = " << x_size << ", x0 = " << x0 << ", xidx = " << xidx << endl;
			cerr << "made it here" << endl;
			long xidx2 = binarySearch(x, x_size, x0, false, true);
			exit(1);
		}
		else return (default_return_value);
	}

	double xidxINT = interpBiLinearNondirect(y, z, f[xidx], y0, z0, y_size, z_size, returnflag);
	double xidxp1INT = interpBiLinearNondirect(y, z, f[xidx+1], y0, z0, y_size, z_size, returnflag);

	return xidxINT + (xidxp1INT-xidxINT)/(x[xidx+1]-x[xidx])*(x0-x[xidx]);
}


//**********************************************************************
double interpQuadriLinearNondirect(double * x, double * y, double * z, double * t, double **** f, double x0, double y0, double z0, double t0,
									long x_size, long y_size, long z_size, long t_size, bool returnflag /*= false*/, double default_return_value /* = 0*/)
{
	//long size = x->size();
	if (x_size==1 && y_size==1 && z_size==1 && t_size==1) {cout<<"interpQuadriLinearNondirect warning: table size = 1"<<endl; return f[0][0][0][0];}

	// assume not close to edges for now...
	// find x's integer index
	//long xidx = floor((x0-x[0])/dx);
	long xidx = binarySearch(x, x_size, x0, true);

	if (xidx < 0 || xidx >= x_size-1)
	{
		if (!returnflag)	//i.e., if returnflag is false, exit
		{
			cerr << "interpQuadriLinearNondirect(): index out of range!  Aborting!" << endl
				<< "interpQuadriLinearNondirect(): x_size = " << x_size << ", x0 = " << x0 << ", xidx = " << xidx << endl;
			exit(1);
		}
		else return (default_return_value);
	}

	double xidxINT = interpTriLinearNondirect(y, z, t, f[xidx], y0, z0, t0, y_size, z_size, t_size, returnflag);
	double xidxp1INT = interpTriLinearNondirect(y, z, t, f[xidx+1], y0, z0, t0, y_size, z_size, t_size, returnflag);

	return xidxINT + (xidxp1INT-xidxINT)/(x[xidx+1]-x[xidx])*(x0-x[xidx]);
}


/////////////////////////////////////////////////////////////////////////////////////////

double interpNewtonDirect(double * x, double * y, double x0, long size)
{
	//size = N + 1 = total number of data points
	int N = size - 1;

	//evaluate coefficients
	double * coeffs = new double [N+1];
	double * uvec = new double [N+2];
	uvec[0] = 0.0;

	for (int i = 1; i<=N+1; i++)
		uvec[i] = y[i-1];

	coeffs[0] = uvec[1];
	for (int m = 1; m<=N; m++)
	{
		for (int i = 1; i<=N-m+1; i++)
			uvec[i] = (uvec[i+1] - uvec[i]) / (x[i+m-1] - x[i-1]);
		coeffs[m] = uvec[1];
	}

	//compute interpolant polynomial
	double P = coeffs[N];
	for(int i = N; i>=1; i--)
		P = P*(x0 - x[i-1]) + coeffs[i-1];
	/*for (int i = 0; i <= N; i++)
	{
		cout << "coeffs[" << i << "] = " << coeffs[i] << endl;
		cout << "uvec[" << i+1 << "] = " << uvec[i+1] << endl;
	}*/

	return P;
}

/////////////////////////////////////////////////////////////////////////////////////////
//**********************************************************************
double interpCubicDirect(double * x, double * y, double x0, long size, bool returnflag /*= false*/, double default_return_value /* = 0*/)
// Returns the interpreted value of y=y(x) at x=x0 using cubic polynomial interpolation method.
// -- x,y: the independent and dependent tables; x is assumed to be equal spaced and increasing
// -- x0: where the interpolation should be performed
{
  //long size = x->size();
  if (size==1) {cout<<"interpCubicDirect warning: table size = 1"; return y[0];}
  double dx = x[1]-x[0]; // increment in x

  // if close to left end:
  if (abs(x0-x[0])<dx*1e-30) return y[0];

  // find x's integer index
  long idx = floor((x0-x[0])/dx);

  //if (idx<0 || idx>=size-1)
  //{
  //  cout    << "interpCubicDirect: x0 out of bounds." << endl
  //          << "x ranges from " << x[0] << " to " << x[size-1] << ", "
  //          << "x0=" << x0 << ", " << "dx=" << dx << ", " << "idx=" << idx << endl;
  //  exit(1);
  //}
	if (idx < 0 || idx >= size-1)
	{
		if (!returnflag)	//i.e., if returnflag is false, exit
		{
			cerr << "interpCubicDirect(): index out of range!  Aborting!" << endl
				<< "interpCubicDirect(): size = " << size << ", x0 = " << x0 << ", " << "dx=" << dx << ", " << "idx=" << idx << endl;
			exit(1);
		}
		else
		{
			idx = (idx<0) ? 0 : size-2;	//uses linear extrapolation
			return y[idx] + (y[idx+1]-y[idx])/(x[idx+1]-x[idx])*(x0-x[idx]);
		}
		//else return default_return_value;
	}

  if (idx==0)
  {
    // use quadratic interpolation at left end
    double A0 = y[0], A1 = y[1], A2 = y[2], deltaX = x0 - x[0]; // deltaX is the increment of x0 compared to the closest lattice point
    return (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX - (3.0*A0-4.0*A1+A2)/(2.0*dx)*deltaX + A0;
  }
  else if (idx==size-2)
  {
    // use quadratic interpolation at right end
    double A0 = y[size-3], A1 = y[size-2], A2 = y[size-1], deltaX = x0 - (x[0] + (idx-1)*dx);
    return (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX - (3.0*A0-4.0*A1+A2)/(2.0*dx)*deltaX + A0;
  }
  else
  {
    // use cubic interpolation
    double A0 = y[idx-1], A1 = y[idx], A2 = y[idx+1], A3 = y[idx+2], deltaX = x0 - (x[0] + idx*dx);
    //cout << A0 << "  " << A1 << "  " << A2 << "  " << A3 << endl;
    return (-A0+3.0*A1-3.0*A2+A3)/(6.0*dx*dx*dx)*deltaX*deltaX*deltaX
            + (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX
            - (2.0*A0+3.0*A1-6.0*A2+A3)/(6.0*dx)*deltaX
            + A1;
  }

}







/////////////////////////////////////////////////////////////////////////////////////////
//**********************************************************************
//double interpCubicNondirect(double * x, double * y, double xx, long size)
double interpCubicNonDirect(double * x, double * y, double xi, long size, bool returnflag /*= false*/, double default_return_value /* = 0*/)
// Returns the interpreted value of y=y(x) at x=x0 using cubic polynomial interpolation method.
// -- x,y: the independent and dependent double x0ables; x is *NOT* assumed to be equal spaced but it has to be increasing
// -- xi: where the interpolation should be performed
{
//cout << "interpCubicNondirect: " << xi << "   " << size << endl;
  //long size = x->size();
  if (size==1) {cout<<"interpCubicNondirect warning: table size = 1"<<endl; return y[0];}

  // if close to left end:
  if (abs(xi-x[0])<(x[1]-x[0])*1e-30) return y[0];

  // find x's integer index
  //long idx = binarySearch(x, xi);
  long idx = binarySearch(x, size, xi, true);
//cout << "interpCubicNondirect: " << idx << endl;

  //if (idx<0 || idx>=size-1)
  //{
  //  cout    << "interpCubicNondirect: x0 out of bounds." << endl
   //         << "x ranges from " << x[0] << " to " << x[size-1] << ", "
  //          << "x0=" << x0 << ", " << "idx=" << idx << endl;
   // exit(1);
  //}
	if (idx < 0 || idx >= size-1)
	{
		if (!returnflag)	//i.e., if returnflag is false, exit
		{
			cerr << "interpCubicNonDirect(): index out of range!  Aborting!" << endl
				<< "interpCubicNonDirect(): size = " << size << ", x0 = " << xi << ", " << "idx=" << idx << endl;
			exit(1);
		}
		else
		{
			idx = (idx<0) ? 0 : size-2;	//uses linear extrapolation
		}
		//else return default_return_value;
	}

  if (idx==0)
  {
    // use linear interpolation at the left end
    return y[0] + (y[1]-y[0])/(x[1]-x[0])*(xi-x[0]);
  }
  else if (idx==size-2)
  {
    // use linear interpolation at the right end
    return y[size-2] + (y[size-1]-y[size-2] )/(x[size-1]-x[size-2] )*(xi-x[size-2]);
  }
  else
  {
    // use cubic interpolation
    long double y0 = y[idx-1], y1 = y[idx], y2 = y[idx+1], y3 = y[idx+2];
    long double y01=y0-y1, y02=y0-y2, y03=y0-y3, y12=y1-y2, y13=y1-y3, y23=y2-y3;
    long double x0 = x[idx-1], x1 = x[idx], x2 = x[idx+1], x3 = x[idx+2];
    long double x01=x0-x1, x02=x0-x2, x03=x0-x3, x12=x1-x2, x13=x1-x3, x23=x2-x3;
    long double x0s=x0*x0, x1s=x1*x1, x2s=x2*x2, x3s=x3*x3;
    long double denominator = x01*x02*x12*x03*x13*x23;
    long double C0, C1, C2, C3;
    C0 = (x0*x02*x2*x03*x23*x3*y1
          + x1*x1s*(x0*x03*x3*y2+x2s*(-x3*y0+x0*y3)+x2*(x3s*y0-x0s*y3))
          + x1*(x0s*x03*x3s*y2+x2*x2s*(-x3s*y0+x0s*y3)+x2s*(x3*x3s*y0-x0*x0s*y3))
          + x1s*(x0*x3*(-x0s+x3s)*y2+x2*x2s*(x3*y0-x0*y3)+x2*(-x3*x3s*y0+x0*x0s*y3))
          )/denominator;
    C1 = (x0s*x03*x3s*y12
          + x2*x2s*(x3s*y01+x0s*y13)
          + x1s*(x3*x3s*y02+x0*x0s*y23-x2*x2s*y03)
          + x2s*(-x3*x3s*y01-x0*x0s*y13)
          + x1*x1s*(-x3s*y02+x2s*y03-x0s*y23)
          )/denominator;
    C2 = (-x0*x3*(x0s-x3s)*y12
          + x2*(x3*x3s*y01+x0*x0s*y13)
          + x1*x1s*(x3*y02+x0*y23-x2*y03)
          + x2*x2s*(-x3*y01-x0*y13)
          + x1*(-x3*x3s*y02+x2*x2s*y03-x0*x0s*y23)
          )/denominator;
    C3 = (x0*x03*x3*y12
          + x2s*(x3*y01+x0*y13)
          + x1*(x3s*y02+x0s*y23-x2s*y03)
          + x2*(-x3s*y01-x0s*y13)
          + x1s*(-x3*y02+x2*y03-x0*y23)
          )/denominator;
//cout << "interpCubicNondirect(list1): " << y0 << "   " << y1 << "   " << y2 << "   " << y3 << "   "
//				<< y01 << "   " << y02 << "   " << y03 << "   " << y12 << "   " << y13 << "   " << y23 << endl;
//cout << "interpCubicNondirect(list2): " << x0 << "   " << x1 << "   " << x2 << "   " << x3 << "   "
//				<< x01 << "   " << x02 << "   " << x03 << "   " << x12 << "   " << x13 << "   " << x23 << endl;
//cout << "interpCubicNondirect(list3): " << x0s << "   " << x1s << "   " << x2s << "   " << x3s << "   " << denominator << endl;
//cout << "interpCubicNondirect(list4): " << C0 << "   " << C1 << "   " << C2 << "   " << C3 << endl;
    return C0 + C1*xi + C2*xi*xi + C3*xi*xi*xi;
  }

}


/////////////////////////////////////////////////////////////////////////////////////////
//**********************************************************************
double interpBiCubicDirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size, bool returnflag /*= false*/, double default_return_value /* = 0*/)
{
	//long size = x->size();
	if (x_size==1 && y_size==1) {cout<<"interpBiCubicDirect warning: table size = 1"<<endl; return z[0][0];}
	double dx = x[1]-x[0]; // increment in x
	double dy = y[1]-y[0]; // increment in y
	// find x's integer index
	long xidx = floor((x0-x[0])/dx);
	long yidx = floor((y0-y[0])/dy);
	
	// check for out-of-bounds points
	if (xidx<0 || xidx>=x_size-1 || yidx<0 || yidx>=y_size-1)
	{
		if (!returnflag)	//i.e., if returnflag is false, exit
		{
			cout << "interpBiCubicDirect: point out of bounds." << endl
				<< "x ranges from " << x[0] << " to " << x[x_size-1] << ", "
				<< "x0=" << x0 << ", " << "dx=" << dx << ", " << "xidx=" << xidx << endl
				<< "y ranges from " << y[0] << " to " << y[y_size-1] << ", "
				<< "y0=" << y0 << ", " << "dy=" << dy << ", " << "yidx=" << yidx << endl;
    			exit(1);
		}
		else return (default_return_value);
	}

  if (xidx==0)
  {
    // use quadratic interpolation at left end
    double A0 = interpCubicDirect(y, z[0], y0, y_size, returnflag);
	double A1 = interpCubicDirect(y, z[1], y0, y_size, returnflag);
	double A2 = interpCubicDirect(y, z[2], y0, y_size, returnflag);
	double deltaX = x0 - x[0]; // deltaX is the increment of x0 compared to the closest lattice point
    return (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX - (3.0*A0-4.0*A1+A2)/(2.0*dx)*deltaX + A0;
  }
  else if (xidx==x_size-2)
  {
    // use quadratic interpolation at right end
    double A0 = interpCubicDirect(y, z[x_size-3], y0, y_size, returnflag);
	double A1 = interpCubicDirect(y, z[x_size-2], y0, y_size, returnflag);
	double A2 = interpCubicDirect(y, z[x_size-1], y0, y_size, returnflag);
	double deltaX = x0 - (x[0] + (xidx-1)*dx);
    return (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX - (3.0*A0-4.0*A1+A2)/(2.0*dx)*deltaX + A0;
  }
  else
  {
    // use cubic interpolation
    double A0 = interpCubicDirect(y, z[xidx-1], y0, y_size, returnflag);
	double A1 = interpCubicDirect(y, z[xidx], y0, y_size, returnflag);
	double A2 = interpCubicDirect(y, z[xidx+1], y0, y_size, returnflag);
	double A3 = interpCubicDirect(y, z[xidx+2], y0, y_size, returnflag);
	double deltaX = x0 - (x[0] + xidx*dx);
    //cout << A0 << "  " << A1 << "  " << A2 << "  " << A3 << endl;
    return (-A0+3.0*A1-3.0*A2+A3)/(6.0*dx*dx*dx)*deltaX*deltaX*deltaX
            + (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX
            - (2.0*A0+3.0*A1-6.0*A2+A3)/(6.0*dx)*deltaX
            + A1;
  }
}


/////////////////////////////////////////////////////////////////////////////////////////
//**********************************************************************
double interpTriCubicDirect(double * x, double * y, double * z, double *** f, double x0, double y0, double z0,
							long x_size, long y_size, long z_size, bool returnflag /*= false*/, double default_return_value /* = 0*/)
{
	//long size = x->size();
	if (x_size==1 && y_size==1 && z_size==1) {cout << "interpTriCubicDirect warning: table size = 1" << endl; return f[0][0][0];}
	double dx = x[1]-x[0]; // increment in x
	double dy = y[1]-y[0]; // increment in y
	double dz = z[1]-z[0]; // increment in z
	// find x's integer index
	long xidx = floor((x0-x[0])/dx);
	long yidx = floor((y0-y[0])/dy);
	long zidx = floor((z0-z[0])/dz);
	
	// check for out-of-bounds points
	if (xidx<0 || xidx>=x_size-1 || yidx<0 || yidx>=y_size-1 || zidx<0 || zidx>=z_size-1)
	{
		if (!returnflag)	//i.e., if returnflag is false, exit
		{
			cout << "interpTriCubicDirect: point out of bounds." << endl
				<< "x ranges from " << x[0] << " to " << x[x_size-1] << ", "
				<< "x0=" << x0 << ", " << "dx=" << dx << ", " << "xidx=" << xidx << endl
				<< "y ranges from " << y[0] << " to " << y[y_size-1] << ", "
				<< "y0=" << y0 << ", " << "dy=" << dy << ", " << "yidx=" << yidx << endl
				<< "z ranges from " << z[0] << " to " << z[z_size-1] << ", "
				<< "z0=" << z0 << ", " << "dz=" << dz << ", " << "zidx=" << zidx << endl;
    			exit(1);
		}
		else return (default_return_value);
	}

  if (xidx==0)
  {
    // use quadratic interpolation at left end
    double A0 = interpBiCubicDirect(y, z, f[0], y0, z0, y_size, z_size, returnflag);
	double A1 = interpBiCubicDirect(y, z, f[1], y0, z0, y_size, z_size, returnflag);
	double A2 = interpBiCubicDirect(y, z, f[2], y0, z0, y_size, z_size, returnflag);
	double deltaX = x0 - x[0]; // deltaX is the increment of x0 compared to the closest lattice point
    return (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX - (3.0*A0-4.0*A1+A2)/(2.0*dx)*deltaX + A0;
  }
  else if (xidx==x_size-2)
  {
    // use quadratic interpolation at right end
    double A0 = interpBiCubicDirect(y, z, f[x_size-3], y0, z0, y_size, z_size, returnflag);
	double A1 = interpBiCubicDirect(y, z, f[x_size-2], y0, z0, y_size, z_size, returnflag);
	double A2 = interpBiCubicDirect(y, z, f[x_size-1], y0, z0, y_size, z_size, returnflag);
	double deltaX = x0 - (x[0] + (xidx-1)*dx);
    return (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX - (3.0*A0-4.0*A1+A2)/(2.0*dx)*deltaX + A0;
  }
  else
  {
    // use cubic interpolation
    double A0 = interpBiCubicDirect(y, z, f[xidx-1], y0, z0, y_size, z_size, returnflag);
	double A1 = interpBiCubicDirect(y, z, f[xidx], y0, z0, y_size, z_size, returnflag);
	double A2 = interpBiCubicDirect(y, z, f[xidx+1], y0, z0, y_size, z_size, returnflag);
	double A3 = interpBiCubicDirect(y, z, f[xidx+2], y0, z0, y_size, z_size, returnflag);
	double deltaX = x0 - (x[0] + xidx*dx);
    //cout << A0 << "  " << A1 << "  " << A2 << "  " << A3 << endl;
    return (-A0+3.0*A1-3.0*A2+A3)/(6.0*dx*dx*dx)*deltaX*deltaX*deltaX
            + (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX
            - (2.0*A0+3.0*A1-6.0*A2+A3)/(6.0*dx)*deltaX
            + A1;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////
//**********************************************************************
double interpBiCubicNonDirectALT(double * x, double * y, double ** f, double xi, double yi, long x_size, long y_size, bool returnflag /*= false*/, double default_return_value /* = 0*/)
{
  //long size = x->size();
  if (x_size==1 && y_size==1) {cout<<"interpBiCubicNonDirectALT warning: table size = 1"<<endl; return f[0][0];}

  // if close to left end:
  if (abs(xi-x[0])<(x[1]-x[0])*1e-30) return interpCubicNonDirect(y, f[0], yi, y_size, returnflag, default_return_value);

  // find x's integer index
  long idx = binarySearch(x, x_size, xi, true);

	if (idx < 0 || idx >= x_size-1)
	{
		if (!returnflag)	//i.e., if returnflag is false, exit
		{
			cerr << "interpBiCubicNonDirectALT(): index out of range!  Aborting!" << endl
				<< "interpBiCubicNonDirectALT(): x_size = " << x_size << ", x0 = " << xi << ", " << "idx=" << idx << endl;
			exit(1);
		}
		else return (default_return_value);
	}

  if (idx==0)
  {
    // use linear interpolation at the left end
	double f0 = interpCubicNonDirect(y, f[0], yi, y_size, returnflag, default_return_value);
	double f1 = interpCubicNonDirect(y, f[1], yi, y_size, returnflag, default_return_value);
    return f0 + (f1-f0)/(x[1]-x[0])*(xi-x[0]);
  }
  else if (idx==x_size-2)
  {
    // use linear interpolation at the right end
	double fnm2 = interpCubicNonDirect(y, f[x_size-2], yi, y_size, returnflag, default_return_value);
	double fnm1 = interpCubicNonDirect(y, f[x_size-1], yi, y_size, returnflag, default_return_value);
    return fnm2 + (fnm1-fnm2 )/(x[x_size-1]-x[x_size-2] )*(xi-x[x_size-2]);
  }
  else
  {
    // use cubic interpolation
    long double f0 = interpCubicNonDirect(y, f[idx-1], yi, y_size, returnflag, default_return_value);
    long double f1 = interpCubicNonDirect(y, f[idx], yi, y_size, returnflag, default_return_value);
    long double f2 = interpCubicNonDirect(y, f[idx+1], yi, y_size, returnflag, default_return_value);
    long double f3 = interpCubicNonDirect(y, f[idx+2], yi, y_size, returnflag, default_return_value);
    long double f01=f0-f1, f02=f0-f2, f03=f0-f3, f12=f1-f2, f13=f1-f3, f23=f2-f3;
    long double x0 = x[idx-1], x1 = x[idx], x2 = x[idx+1], x3 = x[idx+2];
    long double x01=x0-x1, x02=x0-x2, x03=x0-x3, x12=x1-x2, x13=x1-x3, x23=x2-x3;
    long double x0s=x0*x0, x1s=x1*x1, x2s=x2*x2, x3s=x3*x3;
    long double denominator = x01*x02*x12*x03*x13*x23;
    long double C0, C1, C2, C3;
    C0 = (x0*x02*x2*x03*x23*x3*f1
          + x1*x1s*(x0*x03*x3*f2+x2s*(-x3*f0+x0*f3)+x2*(x3s*f0-x0s*f3))
          + x1*(x0s*x03*x3s*f2+x2*x2s*(-x3s*f0+x0s*f3)+x2s*(x3*x3s*f0-x0*x0s*f3))
          + x1s*(x0*x3*(-x0s+x3s)*f2+x2*x2s*(x3*f0-x0*f3)+x2*(-x3*x3s*f0+x0*x0s*f3))
          )/denominator;
    C1 = (x0s*x03*x3s*f12
          + x2*x2s*(x3s*f01+x0s*f13)
          + x1s*(x3*x3s*f02+x0*x0s*f23-x2*x2s*f03)
          + x2s*(-x3*x3s*f01-x0*x0s*f13)
          + x1*x1s*(-x3s*f02+x2s*f03-x0s*f23)
          )/denominator;
    C2 = (-x0*x3*(x0s-x3s)*f12
          + x2*(x3*x3s*f01+x0*x0s*f13)
          + x1*x1s*(x3*f02+x0*f23-x2*f03)
          + x2*x2s*(-x3*f01-x0*f13)
          + x1*(-x3*x3s*f02+x2*x2s*f03-x0*x0s*f23)
          )/denominator;
    C3 = (x0*x03*x3*f12
          + x2s*(x3*f01+x0*f13)
          + x1*(x3s*f02+x0s*f23-x2s*f03)
          + x2*(-x3s*f01-x0s*f13)
          + x1s*(-x3*f02+x2*f03-x0*f23)
          )/denominator;
    return C0 + C1*xi + C2*xi*xi + C3*xi*xi*xi;
  }
}



/////////////////////////////////////////////////////////////////////////////////////////
//**********************************************************************
double interpTriCubicNonDirect(double * x, double * y, double * z, double *** t, double x0, double y0, double z0,
									long x_size, long y_size, long z_size, bool returnflag /*= false*/, double default_return_value /* = 0*/)
{
	//long size = x->size();
	if (x_size==1 && y_size==1) {cout<<"interpBiCubicNonDirectALT warning: table size = 1"<<endl; return t[0][0][0];}
	double dx = x[1]-x[0]; // increment in x
	double dy = y[1]-y[0]; // increment in y
	double dz = z[1]-z[0]; // increment in z
	// find x's integer index
	long xidx = floor((x0-x[0])/dx);
	long yidx = floor((y0-y[0])/dy);
	long zidx = floor((z0-z[0])/dz);
	
	// check for out-of-bounds points
	if (xidx<0 || xidx>=x_size-1 || yidx<0 || yidx>=y_size-1 || zidx<0 || zidx>=z_size-1)
	{
		if (!returnflag)	//i.e., if returnflag is false, exit
		{
			cout << "interpTriCubicNonDirect: point out of bounds." << endl
				<< "x ranges from " << x[0] << " to " << x[x_size-1] << ", "
				<< "x0=" << x0 << ", " << "dx=" << dx << ", " << "xidx=" << xidx << endl
				<< "y ranges from " << y[0] << " to " << y[y_size-1] << ", "
				<< "y0=" << y0 << ", " << "dy=" << dy << ", " << "yidx=" << yidx << endl
				<< "z ranges from " << z[0] << " to " << z[z_size-1] << ", "
				<< "z0=" << z0 << ", " << "dz=" << dz << ", " << "zidx=" << zidx << endl;
    			exit(1);
		}
		else return (default_return_value);
	}

  if (xidx==0)
  {
    // use quadratic interpolation at left end
    double A0 = interpBiCubicNonDirectALT(y, z, t[0], y0, z0, y_size, z_size, returnflag);
	double A1 = interpBiCubicNonDirectALT(y, z, t[1], y0, z0, y_size, z_size, returnflag);
	double A2 = interpBiCubicNonDirectALT(y, z, t[2], y0, z0, y_size, z_size, returnflag);
	double deltaX = x0 - x[0]; // deltaX is the increment of x0 compared to the closest lattice point
    return (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX - (3.0*A0-4.0*A1+A2)/(2.0*dx)*deltaX + A0;
  }
  else if (xidx==x_size-2)
  {
    // use quadratic interpolation at right end
    double A0 = interpBiCubicNonDirectALT(y, z, t[x_size-3], y0, z0, y_size, z_size, returnflag);
	double A1 = interpBiCubicNonDirectALT(y, z, t[x_size-2], y0, z0, y_size, z_size, returnflag);
	double A2 = interpBiCubicNonDirectALT(y, z, t[x_size-1], y0, z0, y_size, z_size, returnflag);
	double deltaX = x0 - (x[0] + (xidx-1)*dx);
    return (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX - (3.0*A0-4.0*A1+A2)/(2.0*dx)*deltaX + A0;
  }
  else
  {
    // use cubic interpolation
    double A0 = interpBiCubicNonDirectALT(y, z, t[xidx-1], y0, z0, y_size, z_size, returnflag);
	double A1 = interpBiCubicNonDirectALT(y, z, t[xidx], y0, z0, y_size, z_size, returnflag);
	double A2 = interpBiCubicNonDirectALT(y, z, t[xidx+1], y0, z0, y_size, z_size, returnflag);
	double A3 = interpBiCubicNonDirectALT(y, z, t[xidx+2], y0, z0, y_size, z_size, returnflag);
	double deltaX = x0 - (x[0] + xidx*dx);
    //cout << A0 << "  " << A1 << "  " << A2 << "  " << A3 << endl;
    return (-A0+3.0*A1-3.0*A2+A3)/(6.0*dx*dx*dx)*deltaX*deltaX*deltaX
            + (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX
            - (2.0*A0+3.0*A1-6.0*A2+A3)/(6.0*dx)*deltaX
            + A1;
  }
}


/////////////////////////////////////////////////////////////////////////////////////////
//**********************************************************************
double interpQuadriCubicNonDirect(double * x, double * y, double * z, double * t, double **** f, double x0, double y0, double z0, double t0,
									long x_size, long y_size, long z_size, long t_size, bool returnflag /*= false*/, double default_return_value /* = 0*/)
{
	//long size = x->size();
	if (x_size==1 && y_size==1 && z_size==1 && t_size==1) {cout<<"interpBiCubicNonDirectALT warning: table size = 1"<<endl; return f[0][0][0][0];}
	double dx = x[1]-x[0]; // increment in x
	double dy = y[1]-y[0]; // increment in y
	double dz = z[1]-z[0]; // increment in z
	double dt = t[1]-t[0]; // increment in t
	// find x's integer index
	long xidx = floor((x0-x[0])/dx);
	long yidx = floor((y0-y[0])/dy);
	long zidx = floor((z0-z[0])/dz);
	long tidx = floor((t0-t[0])/dt);
	
	// check for out-of-bounds points
	if (xidx<0 || xidx>=x_size-1 || yidx<0 || yidx>=y_size-1 || zidx<0 || zidx>=z_size-1|| tidx<0 || tidx>=t_size-1)
	{
		if (!returnflag)	//i.e., if returnflag is false, exit
		{
			cout << "interpQuadriCubicNonDirect: point out of bounds." << endl
				<< "x ranges from " << x[0] << " to " << x[x_size-1] << ", "
				<< "x0=" << x0 << ", " << "dx=" << dx << ", " << "xidx=" << xidx << endl
				<< "y ranges from " << y[0] << " to " << y[y_size-1] << ", "
				<< "y0=" << y0 << ", " << "dy=" << dy << ", " << "yidx=" << yidx << endl
				<< "z ranges from " << z[0] << " to " << z[z_size-1] << ", "
				<< "z0=" << z0 << ", " << "dz=" << dz << ", " << "zidx=" << zidx << endl
				<< "t ranges from " << t[0] << " to " << t[t_size-1] << ", "
				<< "t0=" << t0 << ", " << "dt=" << dt << ", " << "tidx=" << tidx << endl;
    			exit(1);
		}
		else return (default_return_value);
	}

  if (xidx==0)
  {
    // use quadratic interpolation at left end
    double A0 = interpTriCubicNonDirect(y, z, t, f[0], y0, z0, t0, y_size, z_size, t_size, returnflag);
	double A1 = interpTriCubicNonDirect(y, z, t, f[1], y0, z0, t0, y_size, z_size, t_size, returnflag);
	double A2 = interpTriCubicNonDirect(y, z, t, f[2], y0, z0, t0, y_size, z_size, t_size, returnflag);
	double deltaX = x0 - x[0]; // deltaX is the increment of x0 compared to the closest lattice point
    return (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX - (3.0*A0-4.0*A1+A2)/(2.0*dx)*deltaX + A0;
  }
  else if (xidx==x_size-2)
  {
    // use quadratic interpolation at right end
    double A0 = interpTriCubicNonDirect(y, z, t, f[x_size-3], y0, z0, t0, y_size, z_size, t_size, returnflag);
	double A1 = interpTriCubicNonDirect(y, z, t, f[x_size-2], y0, z0, t0, y_size, z_size, t_size, returnflag);
	double A2 = interpTriCubicNonDirect(y, z, t, f[x_size-1], y0, z0, t0, y_size, z_size, t_size, returnflag);
	double deltaX = x0 - (x[0] + (xidx-1)*dx);
    return (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX - (3.0*A0-4.0*A1+A2)/(2.0*dx)*deltaX + A0;
  }
  else
  {
    // use cubic interpolation
    double A0 = interpTriCubicNonDirect(y, z, t, f[xidx-1], y0, z0, t0, y_size, z_size, t_size, returnflag);
	double A1 = interpTriCubicNonDirect(y, z, t, f[xidx], y0, z0, t0, y_size, z_size, t_size, returnflag);
	double A2 = interpTriCubicNonDirect(y, z, t, f[xidx+1], y0, z0, t0, y_size, z_size, t_size, returnflag);
	double A3 = interpTriCubicNonDirect(y, z, t, f[xidx+2], y0, z0, t0, y_size, z_size, t_size, returnflag);
	double deltaX = x0 - (x[0] + xidx*dx);
    //cout << A0 << "  " << A1 << "  " << A2 << "  " << A3 << endl;
    return (-A0+3.0*A1-3.0*A2+A3)/(6.0*dx*dx*dx)*deltaX*deltaX*deltaX
            + (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX
            - (2.0*A0+3.0*A1-6.0*A2+A3)/(6.0*dx)*deltaX
            + A1;
  }
}


/////////////////////////////////////////////////////////////////////////////////////////

//**********************************************************************
double interpPolyDirect(double * x, double * y, double x0, long size)
{
	double result, error;
	polint(x, y, size, x0, &result, &error);
	return result;
}


/////////////////////////////////////////////////////////////////////////////////////////


//**********************************************************************
double interpBiPolyDirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size)
{
	double result, error;
	polin2(x, y, z, x_size, y_size, x0, y0, &result, &error);
	return result;
}



/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////CALL THESE ROUTINES BELOW////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


//**********************************************************************
double interpolate1D(double * x, double * y, double x0, long size, int kind, bool uniform_spacing, bool returnflag /*= false*/, double default_return_value /* = 0*/)
{
// kind == 0: linear interpolation
// kind == 1: cubic interpolation
// kind == 2: polynomial interpolation
	switch (kind)
	{
		case 0:
		{
			if (uniform_spacing)
				return interpLinearDirect(x, y, x0, size, returnflag, default_return_value);
			else
				return interpLinearNondirect(x, y, x0, size, returnflag, default_return_value);
			break;
		}
		case 1:
		{
			if (uniform_spacing)
				return interpCubicDirect(x, y, x0, size, returnflag, default_return_value);
			else
				return interpCubicNonDirect(x, y, x0, size, returnflag, default_return_value);
				//cerr << "Error (interpolate1D): cubic interpolation with non-uniform spacing not supported!" << endl;
			break;
		}
		case 2:
		{
			if (uniform_spacing)
				return interpPolyDirect(x, y, x0, size);
			else
				cerr << "Error (interpolate1D): polynomial interpolation with non-uniform spacing not supported!" << endl;
			break;
		}
		default:
		{
			cerr << "Error (interpolate1D): interpolation kind not supported!" << endl;
			exit(1);
			break;
		}
	}
	return default_return_value;
}


//**********************************************************************
double interpolate2D(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size, int kind, bool uniform_spacing, bool returnflag /*= false*/, double default_return_value /* = 0*/)
{
// kind == 0: linear interpolation
// kind == 1: cubic interpolation
// kind == 2: polynomial interpolation
	switch (kind)
	{
		case 0:
		{
//debugger(__LINE__, __FILE__);
			if (uniform_spacing)
				return interpBiLinearDirect(x, y, z, x0, y0, x_size, y_size, returnflag, default_return_value);
			else
				return interpBiLinearNondirect(x, y, z, x0, y0, x_size, y_size, returnflag, default_return_value);
			break;
		}
		case 1:
		{
			if (uniform_spacing)
				return interpBiCubicDirect(x, y, z, x0, y0, x_size, y_size, returnflag, default_return_value);
			else
				return interpBiCubicNonDirectALT(x, y, z, x0, y0, x_size, y_size, returnflag, default_return_value);
				//cerr << "Error (interpolate2D): cubic interpolation with non-uniform spacing not supported!" << endl;
			break;
		}
		case 2:
		{
			if (uniform_spacing)
				return interpBiPolyDirect(x, y, z, x0, y0, x_size, y_size);
			else
				cerr << "Error (interpolate2D): polynomial interpolation with non-uniform spacing not supported!" << endl;
			break;
		}
		default:
		{
			cerr << "Error (interpolate2D): interpolation kind not supported!" << endl;
			exit(1);
			break;
		}
	}
	return default_return_value;
}



//**********************************************************************
double interpolate3D(double * x, double * y, double * z, double *** f, double x0, double y0, double z0,
			long x_size, long y_size, long z_size, int kind, bool uniform_spacing, bool returnflag /*= false*/, double default_return_value /* = 0*/)
{
// kind == 0: linear interpolation
// kind == 1: cubic interpolation
// kind == 2: polynomial interpolation
	switch (kind)
	{
		case 0:
		{
			if (uniform_spacing)
				return interpTriLinearDirect(x, y, z, f, x0, y0, z0, x_size, y_size, z_size, returnflag, default_return_value);
			else
				return interpTriLinearNondirect(x, y, z, f, x0, y0, z0, x_size, y_size, z_size, returnflag, default_return_value);
			break;
		}
		case 1:
		{
			if (uniform_spacing)
				return interpTriCubicDirect(x, y, z, f, x0, y0, z0, x_size, y_size, z_size, returnflag, default_return_value);
			else
			{
				cerr << "Error (interpolate3D): cubic interpolation not supported!" << endl;
				exit(1);
			}
			break;
		}
		case 2:
		{
			cerr << "Error (interpolate3D): polynomial interpolation not supported!" << endl;
			exit(1);
			break;
		}
		default:
		{
			cerr << "Error (interpolate3D): interpolation kind not supported!" << endl;
			exit(1);
			break;
		}
	}
	return default_return_value;
}


// End of file
