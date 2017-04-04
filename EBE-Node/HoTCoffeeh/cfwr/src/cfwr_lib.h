#ifndef CFWR_LIB_H
#define CFWR_LIB_H

#include<iostream>
#include<cmath>
#include<algorithm>
#include<vector>
#include<stdio.h>

#include "cfwr.h"

using namespace std;

inline int CorrelationFunction::indexer(const int ipt, const int ipphi, const int ipY, const int iqt, const int iqx, const int iqy, const int iqz, const int itrig)
{
	return (
		( ( ( ( ( ( ipt * n_pphi_pts + ipphi ) * n_pY_pts + ipY ) * qtnpts + iqt ) * qxnpts + iqx ) * qynpts + iqy ) * qznpts + iqz ) * 2 + itrig
	);
}

inline int CorrelationFunction::FM_indexer(const int ipY, const int iqt, const int iqx, const int iqy, const int iqz)
{
	return (
		( ( ( ipY * qtnpts + iqt ) * qxnpts + iqx ) * qynpts + iqy ) * qznpts + iqz
	);
}

/*inline int CorrelationFunction::arb_indexer(const vector<int> indices, const vector<int> sizes)
{
	int n = sizes.size();
	int flat = indices[0];
	for (int i = 0; i < n; ++i)
		flat = sizes[i] * flat + indices[i+1];
	return (flat);
}*/

inline void CorrelationFunction::set_to_zero(double * array, size_t arraylength)
{
	for (size_t arrayidx=0; arrayidx<arraylength; ++arrayidx) array[arrayidx] = 0.0;
}

inline double CorrelationFunction::lin_int(double x_m_x1, double one_by_x2_m_x1, double f1, double f2)
{
	return ( f1 + (f2 - f1) * x_m_x1 * one_by_x2_m_x1 );
}

//dots 4-vectors together
inline double CorrelationFunction::dot_four_vectors(double * a, double * b)
{
	double sum = a[0]*b[0];
	for (size_t i=1; i<4; ++i) sum -= a[i]*b[i];
	return (sum);
}

inline void CorrelationFunction::addElementToQueue(priority_queue<pair<double, size_t> >& p, pair<double, size_t> elem, size_t max_size)
{
	if ( ( max_size <= p.size() ) && ( elem >= p.top() ) )
		return; // nothing to do.
	p.push(elem);
	if( max_size < p.size() )
		p.pop();
	return;
}

//*********************************************************************
//  Function to return the i'th best-fit parameter
inline double CorrelationFunction::get_fit_results(int i, gsl_multifit_fdfsolver * solver_ptr)
{
	return gsl_vector_get (solver_ptr->x, i);
}

//*********************************************************************
//  Function to retrieve the square root of the diagonal elements of
//   the covariance matrix.
inline double CorrelationFunction::get_fit_err (int i, gsl_matrix * covariance_ptr)
{
	return sqrt (gsl_matrix_get (covariance_ptr, i, i));
}

/************************************************************************/

//End of file

#endif
