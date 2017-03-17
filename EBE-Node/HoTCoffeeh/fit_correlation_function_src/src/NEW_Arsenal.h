#ifndef ARSENAL_H
#define ARSENAL_H

using namespace std;

struct interpolationOptions
{
	bool returnFlag;			//whether to return a value or not
	double returnValue;			//default value to return
	bool assumePeriodic;		//assumes range is periodic
	double period;				//default period to assume
	bool useExtrapolation;		//assumes range is periodic
};

unsigned long int random_seed();
int sgn(double val);
void logspace(double * x, double a, double b, int n);
void linspace(double * x, double a, double b, int n);
void linspace(vector<double> & x, double a, double b);
void stratify_npts(double a, double b, int n1, int npts, double * x);
void scalepoints(double * x, double a, double b, double scale, int n);

//miscellaneous functions needed for interpolation routines
long binarySearch(double * A, int length, double value, bool skip_out_of_range = true, bool verbose = false);

//individual interpolation routines
double interpLinearDirect(double * x, double * y, double x0, long size, void * interpOpts_ptr);
double interpLinearNondirect(double * x, double * y, double x0, long size, void * interpOpts_ptr);

double interpBiLinearDirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size, void * interpOpts_ptr);
double interpBiLinearNondirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size, void * interpOpts_ptr);

double interpTriLinearDirect(double * x, double * y, double * z, double *** f, double x0, double y0, double z0, long x_size, long y_size, long z_size, void * interpOpts_ptr);
double interpTriLinearNondirect(double * x, double * y, double * z, double *** f, double x0, double y0, double z0, long x_size, long y_size, long z_size, void * interpOpts_ptr);

double interpCubicDirect(double * x, double * y, double x0, long size, void * interpOpts_ptr);
double interpCubicNonDirect(double * x, double * y, double x0, long size, void * interpOpts_ptr);

double interpBiCubicDirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size, void * interpOpts_ptr);
double interpBiCubicNonDirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size, void * interpOpts_ptr);

double interpTriCubicDirect(double * x, double * y, double * z, double *** f, double x0, double y0, double z0,
							long x_size, long y_size, long z_size, void * interpOpts_ptr);
double interpTriCubicNonDirect(double * x, double * y, double * z, double *** t, double x0, double y0, double z0,
									long x_size, long y_size, long z_size, void * interpOpts_ptr);

//main interpolation routines
double interpolate1D(double * x, double * y, double x0, long size, int kind, bool uniform_spacing, void * interpOpts_ptr);
double interpolate2D(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size, int kind, bool uniform_spacing, void * interpOpts_ptr);
double interpolate3D(double * x, double * y, double * z, double *** f, double x0, double y0, double z0,
			long x_size, long y_size, long z_size, int kind, bool uniform_spacing, void * interpOpts_ptr);

#endif
