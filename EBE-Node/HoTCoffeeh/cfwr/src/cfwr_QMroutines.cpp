#include<iostream>
#include<cmath>
#include<algorithm>
#include<vector>
#include<stdio.h>

#include <gsl/gsl_matrix.h>

#include "cfwr.h"

using namespace std;

void CorrelationFunction::Get_QM_HBTradii()
{
	*global_out_stream_ptr << "--> Getting HBT radii by q-moments method" << endl;

	Allocate_fleshed_out_CF();

	for (int ipt = 0; ipt < n_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	{
		*global_out_stream_ptr << "   --> Doing pT = " << SP_pT[ipt] << ", pphi = " << SP_pphi[ipphi] << "..." << endl;
		
		double sample_scale = 3.0;
		Flesh_out_CF(ipt, ipphi, sample_scale);

		//finally, do fits, depending on what kind you want to do
		Get_q_moments(fleshed_out_CF, ipt, ipphi);
	}

	Delete_fleshed_out_CF();

	return;
}

void CorrelationFunction::Get_q_moments(double *** current_C_slice, int ipt, int ipphi)
{
	//coordinates: o-s-l
	double invR2ij[3][3];
	double norm = 0.0;

	double ckp = cos_SP_pphi[ipphi], skp = sin_SP_pphi[ipphi];

	for (int iqx = 0; iqx < new_nqpts; ++iqx)
	for (int iqy = 0; iqy < new_nqpts; ++iqy)
	for (int iqz = 0; iqz < new_nqpts; ++iqz)
	{
		double C_at_q = current_C_slice[iqx][iqy][iqz];
		double qx0 = qx_fleshed_out_pts[iqx];
		double qy0 = qy_fleshed_out_pts[iqy];
		double qz0 = qz_fleshed_out_pts[iqz];
		cout << "QM CHECK: " << SP_pT[ipt] << "   " << SP_pphi[ipphi] << "   " << qx0 << "   " << qy0 << "   " << qz0 << "   " << C_at_q << endl;

		double factor1 = 1.;
		double factor2 = 1.;
		double factor3 = 1.;
		if (iqx == 0 || iqx == new_nqpts-1) factor1 = 0.5;
		if (iqy == 0 || iqy == new_nqpts-1) factor2 = 0.5;
		if (iqz == 0 || iqz == new_nqpts-1) factor3 = 0.5;

		double q_osl_point[3] = { qx0 * ckp + qy0 * skp, -qx0 * skp + qy0 * ckp, qz0 };

		for (int iq1 = 0; iq1 < 3; ++iq1)
		for (int iq2 = 0; iq2 < 3; ++iq2)
			invR2ij[iq1][iq2] += 2.0 * factor1 * factor2 * factor3 * ( C_at_q - 1.0 ) * q_osl_point[iq1] * q_osl_point[iq2];

		norm += factor1 * factor2 * factor3 * ( C_at_q - 1.0 );
	}

	//normalize appropriately
	for (int iq1 = 0; iq1 < 3; ++iq1)
	for (int iq2 = 0; iq2 < 3; ++iq2)
		invR2ij[iq1][iq2] /= norm;

	//now use GSL routines to invert matrix
	gsl_matrix * mat = gsl_matrix_alloc (3, 3);
	gsl_matrix * R2ij = gsl_matrix_alloc (3, 3);
  
	for (int iq1 = 0; iq1 < 3; ++iq1)
	for (int iq2 = 0; iq2 < 3; ++iq2)
		gsl_matrix_set (mat, iq1, iq2, invR2ij[iq1][iq2]);

	gsl_permutation * p = gsl_permutation_alloc (3);

	int s, status = 0;

	status += gsl_linalg_LU_decomp (mat, p, &s);

	status += gsl_linalg_LU_invert (mat, p, R2ij);

	R2_out_QM[ipt][ipphi] = gsl_matrix_get (R2ij, 0, 0) * hbarC * hbarC;
	R2_side_QM[ipt][ipphi] = gsl_matrix_get (R2ij, 1, 1) * hbarC * hbarC;
	R2_long_QM[ipt][ipphi] = gsl_matrix_get (R2ij, 2, 2) * hbarC * hbarC;
	R2_outside_QM[ipt][ipphi] = gsl_matrix_get (R2ij, 0, 1) * hbarC * hbarC;
	R2_sidelong_QM[ipt][ipphi] = gsl_matrix_get (R2ij, 1, 2) * hbarC * hbarC;
	R2_outlong_QM[ipt][ipphi] = gsl_matrix_get (R2ij, 0, 2) * hbarC * hbarC;

	status += gsl_linalg_LU_decomp (R2ij, p, &s);

	double det = gsl_linalg_LU_det (R2ij, s);

	lambda_QM[ipt][ipphi] = sqrt(abs(det)/(M_PI*M_PI*M_PI)) * norm;

	gsl_permutation_free (p);
	gsl_matrix_free (mat);
	gsl_matrix_free (R2ij);

	return;
}

//End of file
