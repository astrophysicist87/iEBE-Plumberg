#ifndef CFWR_LIB_H
#define CFWR_LIB_H

#include<iostream>
#include<cmath>
#include<algorithm>
#include<vector>
#include<stdio.h>

#include "cfwr.h"

using namespace std;

inline int CorrelationFunction::indexer(const int ipT, const int ipphi, const int iqt, const int iqx, const int iqy, const int iqz, const int itrig)
{
	return (
		( ( ( ( ( ipT * n_pphi_pts + ipphi ) * qtnpts + iqt ) * qxnpts + iqx ) * qynpts + iqy ) * qznpts + iqz ) * ntrig + itrig
	);
}

inline int CorrelationFunction::indexer(const int ipT, const int ipphi, const int ipY, const int iqt, const int iqx, const int iqy, const int iqz, const int itrig)
{
	return (
		( ( ( ( ( ( ipT * n_pphi_pts + ipphi ) * n_pY_pts + ipY )  * qtnpts + iqt ) * qxnpts + iqx ) * qynpts + iqy ) * qznpts + iqz ) * ntrig + itrig
	);
}

inline int CorrelationFunction::fixQTQZ_indexer(const int ipT, const int ipphi, const int ipY, const int iqx, const int iqy, const int itrig)
{
	return (
		( ( ( ( ipT * n_pphi_pts + ipphi ) * n_pY_pts + ipY )  * qxnpts + iqx ) * qynpts + iqy ) * ntrig + itrig
	);
}

inline int CorrelationFunction::FM_indexer(const int ipY, const int iqt, const int iqx, const int iqy, const int iqz)
{
	return (
		( ( ( ipY * qtnpts + iqt ) * qxnpts + iqx ) * qynpts + iqy ) * qznpts + iqz
	);
}

inline int CorrelationFunction::HDF_indexer(const int ir, const int iqt, const int iqz)
{
	return (
		( ir * qtnpts + iqt ) * qznpts + iqz
	);
}

inline int CorrelationFunction::mom_indexer(const int ipT, const int ipphi, const int ipY)
{
	return (
		( ipT * n_pphi_pts + ipphi ) * n_pY_pts + ipY
	);
}

inline int CorrelationFunction::indexer2(const int ipT, const int ipphi, const int ipY, const int iqt, const int iqx, const int iqy, const int iqz)
{
	return (
		( ( ( ( ( ipT * n_pphi_pts + ipphi ) * n_pY_pts + ipY ) * qtnpts + iqt ) * qxnpts + iqx ) * qynpts + iqy ) * qznpts + iqz
	);
}

inline int CorrelationFunction::indexer3(const int ipT, const int ipphi, const int isurf)
{
	return (
		( isurf * n_pT_pts + ipT ) * n_pphi_pts + ipphi
	);
}

inline int CorrelationFunction::indexer4(const int ipT, const int ipphi, const int iqx, const int iqy)
{
	return (
		//( ( ipT * n_pphi_pts + ipphi ) * qxnpts + iqx ) * qynpts + iqy
		( (  iqx * qynpts + iqy ) * n_pT_pts + ipT ) * n_pphi_pts + ipphi
	);
}

inline int CorrelationFunction::NB2_indexer(const int iv, const int izeta)
{
	return (
		iv * n_zeta_pts + izeta
	);
}

inline int CorrelationFunction::NB3_indexer(const int is, const int iv, const int izeta)
{
	return (
		( is * n_v_pts + iv ) * n_zeta_pts + izeta
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

//inline double CorrelationFunction::lin_int(double x_m_x1, double one_by_x2_m_x1, double f1, double f2)
inline double lin_int(double x_m_x1, double one_by_x2_m_x1, double f1, double f2)
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

	//while( max_size < p.size() )
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
