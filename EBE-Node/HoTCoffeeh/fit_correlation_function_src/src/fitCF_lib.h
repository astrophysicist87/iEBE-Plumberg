#ifndef FITCF_LIB_H
#define FITCF_LIB_H

#include<iostream>
#include<cmath>
#include<algorithm>
#include<vector>
#include<stdio.h>

#include "fitCF.h"

using namespace std;

inline int FitCF::indexer(const int ipt, const int ipphi, const int iqt, const int iqx, const int iqy, const int iqz, const int itrig)
{
	return (
		( ( ( ( ( ipt * n_interp_pphi_pts + ipphi ) * qtnpts + iqt ) * qxnpts + iqx ) * qynpts + iqy ) * qznpts + iqz ) * 2 + itrig
	);
}

inline int FitCF::indexer2(const int iKT, const int iKPHI, const int iqt, const int iqx, const int iqy, const int iqz, const int itrig)
{
	return (
		( ( ( ( ( iKT * n_localp_phi + iKPHI ) * qtnpts + iqt ) * qxnpts + iqx ) * qynpts + iqy ) * qznpts + iqz ) * 2 + itrig
	);
}

inline double FitCF::lin_int(double x_m_x1, double one_by_x2_m_x1, double f1, double f2)
{
	return ( f1 + (f2 - f1) * x_m_x1 * one_by_x2_m_x1 );
}

//*********************************************************************
//  Function to return the i'th best-fit parameter
inline double FitCF::get_fit_results(int i, gsl_multifit_fdfsolver * solver_ptr)
{
	return gsl_vector_get (solver_ptr->x, i);
}

//*********************************************************************
//  Function to retrieve the square root of the diagonal elements of
//   the covariance matrix.
inline double FitCF::get_fit_err (int i, gsl_matrix * covariance_ptr)
{
	return sqrt (gsl_matrix_get (covariance_ptr, i, i));
}

/************************************************************************/

//End of file

#endif
