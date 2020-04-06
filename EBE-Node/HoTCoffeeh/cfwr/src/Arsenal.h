#ifndef ARSENAL_H
#define ARSENAL_H

#include "stdlib.h"
#include <vector>
#include <string>

using namespace std;

unsigned long int random_seed();
int sgn(double val);
void logspace(double * x, double a, double b, int n);
void linspace(double * x, double a, double b, int n);
void linspace(vector<double> & x, double a, double b);
void stratify_npts(double a, double b, int n1, int npts, double * x);
void stratify_npts( double min, double a, double b, int n1, int npts, double * x );
void scalepoints(double * x, double a, double b, double scale, int n);

//miscellaneous functions needed for interpolation routines
long binarySearch(double * A, int length, double value, bool skip_out_of_range = true, bool verbose = false);

//ParameterReader routines
double stringToDouble(string);
vector<double> stringToDoubles(string);
string toLower(string str);
string trim(string str);
vector< vector<double>* >* readBlockData(istream &stream_in);
void releaseBlockData(vector< vector<double>* >* data);
void print_progressbar(double percentage, int length=50, string symbol="#");
void display_logo(int which=1);

//individual interpolation routines
double interpLinearDirect(double * x, double * y, double x0, long size,
				bool returnflag = false, double default_return_value = 0.0);
double interpLinearNondirect(double * x, double * y, double x0, long size,
				bool returnflag = false, double default_return_value = 0.0);
double interpBiLinearDirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size,
				bool returnflag = false, double default_return_value = 0.0);
double interpBiLinearNondirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size,
				bool returnflag = false, double default_return_value = 0.0);
double interpTriLinearDirect(double * x, double * y, double * z, double *** f, double x0, double y0, double z0, long x_size, long y_size, long z_size,
				bool returnflag = false, double default_return_value = 0.0);
double interpTriLinearNondirect(double * x, double * y, double * z, double *** f, double x0, double y0, double z0, long x_size, long y_size, long z_size,
				bool returnflag = false, double default_return_value = 0.0);
double interpCubicDirect(double * x, double * y, double x0, long size,
				bool returnflag = false, double default_return_value = 0.0);
//double interpCubicNondirect(double * x, double * y, double xx, long size);
double interpCubicNonDirect(double * x, double * y, double x0, long size,
				bool returnflag = false, double default_return_value = 0.0);
double interpBiCubicDirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size,
				bool returnflag = false, double default_return_value = 0.0);
double interpBiCubicNonDirectALT(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size,
				bool returnflag = false, double default_return_value = 0.0);

double interpTriCubicDirect(double * x, double * y, double * z, double *** f, double x0, double y0, double z0,
							long x_size, long y_size, long z_size, bool returnflag /*= false*/, double default_return_value /* = 0*/);

double interpTriCubicNonDirect(double * x, double * y, double * z, double *** t, double x0, double y0, double z0,
									long x_size, long y_size, long z_size, bool returnflag=false, double default_return_value=0);
double interpQuadriCubicNonDirect(double * x, double * y, double * z, double * t, double **** f, double x0, double y0, double z0, double t0,
									long x_size, long y_size, long z_size, long t_size, bool returnflag=false, double default_return_value=0);
double interpQuadriLinearNondirect(double * x, double * y, double * z, double * t, double **** f, double x0, double y0, double z0, double t0,
									long x_size, long y_size, long z_size, long t_size, bool returnflag=false, double default_return_value=0);

//main interpolation routines
double interpolate1D(double * x, double * y, double x0, long size, int kind, bool uniform_spacing,
			bool returnflag = false, double default_return_value = 0.0);
double interpolate2D(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size, int kind, bool uniform_spacing,
			bool returnflag = false, double default_return_value = 0.0);
double interpolate3D(double * x, double * y, double * z, double *** f, double x0, double y0, double z0,
			long x_size, long y_size, long z_size, int kind, bool uniform_spacing,
			bool returnflag = false, double default_return_value = 0.0);

#endif
