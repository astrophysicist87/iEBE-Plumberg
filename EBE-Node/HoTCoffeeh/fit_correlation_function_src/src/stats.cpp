#ifndef STATS_H
#define STATS_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

#include "stats.h"

double get_median(double * array, int length)
{
	vector<double> tmp ( array, array+length );
	sort( tmp.begin(), tmp.end() );
	return ( 0.5 * ( tmp[(length-1)/2] + tmp[length/2] ) );
}

double get_median_absolute_deviation(double * array, int length, double * median)
{
	*median = get_median( array, length );

	double * dev = new double [length];
	for (int i = 0; i < length; ++i)
		dev[i] = abs( array[i] - *median );

	double result = get_median( dev, length );

	delete [] dev;

	return ( result );
}

void find_outliers_Hampel(double * array, int length, bool * results, double * med, double Hampel_factor)
{
	//double factor - use of this factor defines Hampel method (Hampel 1985)

	double mad = get_median_absolute_deviation(array, length, med);

	//assumes results has already been allocated elsewhere
	if (abs(*med) <= 1.e-4 || mad <= 1.e-4)	//if all results are very close to zero or the spread is basically zero, don't return outliers
		for (int i = 0; i < length; ++i)
			results[i] = false;
	else
		for (int i = 0; i < length; ++i)
			results[i] = ( abs( array[i] - *med ) > Hampel_factor * mad );	//is true if array[i] is an outlier (by this detection algorithm)

	return;
}

int place_in_range(int idx, int minidx, int maxidx)
{
	int n = maxidx - minidx + 1;
	while (idx > maxidx)
		idx -= n;
	while (idx < minidx)
		idx += n;
	return (idx);
}

//applies Hampel filter only in window of fixed width
void find_outliers_window_Hampel(double * array, int length, bool * results, double * meds, double Hampel_factor, int width)
{
	double * window_array = new double [width];

	//window must be odd
	int hw = (width-1)/2;
	int min0 = place_in_range( 0 - hw, 0, length - 1 );
	int max0 = place_in_range( 0 + hw, 0, length - 1 );

	//set initial set of indices
	vector<int> inds;
	for (int im = min0; im < length; ++im)
		inds.push_back(im);
	for (int im = 0; im <= max0; ++im)
		inds.push_back(im);

	for (int i = 0; i < length; ++i)
	{
cout << "For element i = " << i << endl;
		for (int iw = 0; iw < width; ++iw)
		{
			window_array[iw] = array[inds[iw]];
cout << "\t" << iw << "   " << inds[iw] << "   " << array[inds[iw]] << endl;
		}

		double window_mad = get_median_absolute_deviation(window_array, width, &meds[i]);

		results[i] = ( abs( array[i] - meds[i] ) > Hampel_factor * window_mad );	//is true if array[i] is an outlier (by this detection algorithm)

cout << "Working with factor, median and m.a.d.: " << Hampel_factor << "   " << meds[i] << "   " << window_mad << " --> " << results[i] << endl;

		//update inds
		inds.erase(inds.begin());
		inds.push_back( place_in_range( i + hw + 1, 0, length - 1 ) );
cout << endl;
	}

	delete [] window_array;

	return;
}


#endif

//End of file
