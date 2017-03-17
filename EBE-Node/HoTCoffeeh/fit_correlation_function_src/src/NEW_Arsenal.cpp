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

namespace arsDefs
{
	const double pi = 3.1415926535897932384626433;
	const bool defaultReturnFlag = false;
	const bool defaultAssumePeriodic = false;
	const bool defaultUseExtrapolation = false;
	const double defaultReturnValue = 0.0;
	const double defaultPeriod = 2.*pi;
}

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

/////////////////////////////////////////////////////////////////////////////////////////

//**********************************************************************
double interpLinearDirect(double * x, double * y, double x0, long size, void * interpOpts_ptr)
// Returns the interpreted value of y=y(x) at x=x0 using linear interpolation method.
// -- x,y: the independent and dependent tables; x is assumed to be equal spaced and increasing
// -- x0: where the interpolation should be performed
{
	interpolationOptions interpOpts = * (struct interpolationOptions *) interpOpts_ptr;

	//long size = x->size();
	if (size==1) {cout<<"interpLinearDirect warning: table size = 1"<<endl; return y[0];}
	double dx = x[1]-x[0]; // increment in x

	// if close to left end:
	if (abs(x0-x[0])<dx*1e-30) return y[0];

	// find x's integer index
	long idx = floor((x0-x[0])/dx);

	if (idx<0 || idx>=size-1)
	{
		if (!interpOpts.returnFlag)	//i.e., if returnflag is false, exit
		{
			cout    << "interpLinearDirect: x0 out of bounds." << endl
				<< "x ranges from " << x[0] << " to " << x[size-1] << ", "
				<< "x0=" << x0 << ", " << "dx=" << dx << ", " << "idx=" << idx << endl;
			exit(1);
		}
		else if (interpOpts.useExtrapolation)
		{
			idx = (idx<0) ? 0 : size-2;	//uses extrapolation
		}
		else return default_return_value;
	}

  return y[idx] + (y[idx+1]-y[idx])/dx*(x0-x[idx]);
}

//**********************************************************************
double interpLinearNondirect(double * x, double * y, double x0, long size, void * interpOpts_ptr)
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
double interpBiLinearDirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size, void * interpOpts_ptr)
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
double interpBiLinearNondirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size, void * interpOpts_ptr)
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
double interpTriLinearDirect(double * x, double * y, double * z, double *** f, double x0, double y0, double z0,
								long x_size, long y_size, long z_size, void * interpOpts_ptr)
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
double interpTriLinearNondirect(double * x, double * y, double * z, double *** f, double x0, double y0, double z0,
									long x_size, long y_size, long z_size, void * interpOpts_ptr)
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


/////////////////////////////////////////////////////////////////////////////////////////
//**********************************************************************
double interpCubicDirect(double * x, double * y, double x0, long size, void * interpOpts_ptr)
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
		else return (default_return_value);
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
double interpCubicNonDirect(double * x, double * y, double xi, long size, void * interpOpts_ptr)
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
		else return (default_return_value);
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
double interpBiCubicDirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size, void * interpOpts_ptr)
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
							long x_size, long y_size, long z_size, void * interpOpts_ptr)
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
double interpBiCubicNonDirect(double * x, double * y, double ** f, double xi, double yi, long x_size, long y_size, void * interpOpts_ptr)
{
  //long size = x->size();
  if (x_size==1 && y_size==1) {cout<<"interpBiCubicNonDirect warning: table size = 1"<<endl; return f[0][0];}

  // if close to left end:
  if (abs(xi-x[0])<(x[1]-x[0])*1e-30) return interpCubicNonDirect(y, f[0], yi, y_size, returnflag, default_return_value);

  // find x's integer index
  long idx = binarySearch(x, x_size, xi, true);

	if (idx < 0 || idx >= x_size-1)
	{
		if (!returnflag)	//i.e., if returnflag is false, exit
		{
			cerr << "interpBiCubicNonDirect(): index out of range!  Aborting!" << endl
				<< "interpBiCubicNonDirect(): x_size = " << x_size << ", x0 = " << xi << ", " << "idx=" << idx << endl;
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
									long x_size, long y_size, long z_size, void * interpOpts_ptr)
{
	//long size = x->size();
	if (x_size==1 && y_size==1) {cout<<"interpBiCubicNonDirect warning: table size = 1"<<endl; return t[0][0][0];}
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
    double A0 = interpBiCubicNonDirect(y, z, t[0], y0, z0, y_size, z_size, returnflag);
	double A1 = interpBiCubicNonDirect(y, z, t[1], y0, z0, y_size, z_size, returnflag);
	double A2 = interpBiCubicNonDirect(y, z, t[2], y0, z0, y_size, z_size, returnflag);
	double deltaX = x0 - x[0]; // deltaX is the increment of x0 compared to the closest lattice point
    return (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX - (3.0*A0-4.0*A1+A2)/(2.0*dx)*deltaX + A0;
  }
  else if (xidx==x_size-2)
  {
    // use quadratic interpolation at right end
    double A0 = interpBiCubicNonDirect(y, z, t[x_size-3], y0, z0, y_size, z_size, returnflag);
	double A1 = interpBiCubicNonDirect(y, z, t[x_size-2], y0, z0, y_size, z_size, returnflag);
	double A2 = interpBiCubicNonDirect(y, z, t[x_size-1], y0, z0, y_size, z_size, returnflag);
	double deltaX = x0 - (x[0] + (xidx-1)*dx);
    return (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX - (3.0*A0-4.0*A1+A2)/(2.0*dx)*deltaX + A0;
  }
  else
  {
    // use cubic interpolation
    double A0 = interpBiCubicNonDirect(y, z, t[xidx-1], y0, z0, y_size, z_size, returnflag);
	double A1 = interpBiCubicNonDirect(y, z, t[xidx], y0, z0, y_size, z_size, returnflag);
	double A2 = interpBiCubicNonDirect(y, z, t[xidx+1], y0, z0, y_size, z_size, returnflag);
	double A3 = interpBiCubicNonDirect(y, z, t[xidx+2], y0, z0, y_size, z_size, returnflag);
	double deltaX = x0 - (x[0] + xidx*dx);
    //cout << A0 << "  " << A1 << "  " << A2 << "  " << A3 << endl;
    return (-A0+3.0*A1-3.0*A2+A3)/(6.0*dx*dx*dx)*deltaX*deltaX*deltaX
            + (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX
            - (2.0*A0+3.0*A1-6.0*A2+A3)/(6.0*dx)*deltaX
            + A1;
  }
}


/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////CALL THESE ROUTINES BELOW////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


//**********************************************************************
double interpolate1D(double * x, double * y, double x0, long size, int kind, bool uniform_spacing, void * interpOpts_ptr)
{
// kind == 0: linear interpolation
// kind == 1: cubic interpolation
	interpolationOptions * options;
	
{
	bool returnFlag = arsDefs::defaultReturnFlag;				//whether to return a value or not
	double returnValue = arsDefs::defaultReturnValue;			//default value to return
	bool assumePeriodic = arsDefs::defaultAssumePeriodic;		//assumes range is periodic
	double period = arsDefs::defaultPeriod;						//default period to assume
	bool useExtrapolation arsDefs::defaultUseExtrapolation;		//assumes range is periodic
};


	switch (kind)
	{
		case 0:
		{
			if (uniform_spacing)
				return interpLinearDirect(x, y, x0, size, interpOpts_ptr);
			else
				return interpLinearNondirect(x, y, x0, size, interpOpts_ptr);
			break;
		}
		case 1:
		{
			if (uniform_spacing)
				return interpCubicDirect(x, y, x0, size, interpOpts_ptr);
			else
				return interpCubicNonDirect(x, y, x0, size, interpOpts_ptr);
				//cerr << "Error (interpolate1D): cubic interpolation with non-uniform spacing not supported!" << endl;
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
double interpolate2D(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size, int kind, bool uniform_spacing, void * interpOpts_ptr)
{
// kind == 0: linear interpolation
// kind == 1: cubic interpolation
	switch (kind)
	{
		case 0:
		{
			if (uniform_spacing)
				return interpBiLinearDirect(x, y, z, x0, y0, x_size, y_size, interpOpts_ptr);
			else
				return interpBiLinearNondirect(x, y, z, x0, y0, x_size, y_size, interpOpts_ptr);
			break;
		}
		case 1:
		{
			if (uniform_spacing)
				return interpBiCubicDirect(x, y, z, x0, y0, x_size, y_size, interpOpts_ptr);
			else
				return interpBiCubicNonDirect(x, y, z, x0, y0, x_size, y_size, interpOpts_ptr);
				//cerr << "Error (interpolate2D): cubic interpolation with non-uniform spacing not supported!" << endl;
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
			long x_size, long y_size, long z_size, int kind, bool uniform_spacing, void * interpOpts_ptr)
{
// kind == 0: linear interpolation
// kind == 1: cubic interpolation
	switch (kind)
	{
		case 0:
		{
			if (uniform_spacing)
				return interpTriLinearDirect(x, y, z, f, x0, y0, z0, x_size, y_size, z_size, interpOpts_ptr);
			else
				return interpTriLinearNondirect(x, y, z, f, x0, y0, z0, x_size, y_size, z_size, interpOpts_ptr);
			break;
		}
		case 1:
		{
			if (uniform_spacing)
				return interpTriCubicDirect(x, y, z, f, x0, y0, z0, x_size, y_size, z_size, interpOpts_ptr);
			else
			{
				cerr << "Error (interpolate3D): cubic interpolation not supported!" << endl;
				exit(1);
			}
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
