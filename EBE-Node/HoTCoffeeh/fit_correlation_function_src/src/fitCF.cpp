#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<set>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<queue>
#include<map>
#include<cstdlib>
#include<numeric>

#include "fitCF.h"
#include "fitCF_lib.h"
#include "Arsenal.h"
#include "Stopwatch.h"
#include "CPStopwatch.h"
#include "gauss_quadrature.h"

using namespace std;

// only need to calculated interpolation grid of spacetime moments for each resonance, NOT each decay channel!
bool recycle_previous_moments = false;
bool recycle_similar_moments = false;
int reso_particle_id_of_moments_to_recycle = -1;
string reso_name_of_moments_to_recycle = "NULL";
string current_decay_channel_string = "NULL";

template < typename T >
void check_for_NaNs(string variable_name, const T variable_value, ofstream& localout)
{
	if (isnan(variable_value))
		localout << "ERROR: " << variable_name << " = " << variable_value << endl;
	return;
}

double FitCF::place_in_range(double phi, double min, double max)
{
	while (phi < min || phi > max)
	{
		if (phi < min) phi += twopi;
		else phi -= twopi;
	}

	return (phi);
}

//**************************************************************
//**************************************************************

void FitCF::R2_Fourier_transform(int iKT, double plane_psi, int mode)
{
	const int interpMode = 1;
	//int mode: 0 - GF, 1 - QM
	for(int Morder = 0; Morder < n_order; ++Morder)
	{
		double cos_mK_phi[n_localp_phi], sin_mK_phi[n_localp_phi];

		for(int iKphi = 0; iKphi < n_localp_phi; ++iKphi)
		{
			cos_mK_phi[iKphi] = cos(Morder*(K_phi[iKphi] - plane_psi));
			sin_mK_phi[iKphi] = sin(Morder*(K_phi[iKphi] - plane_psi));
		}

		double temp_sum_side_cos = 0.0, temp_sum_side_sin = 0.0;
		double temp_sum_out_cos = 0.0, temp_sum_out_sin = 0.0;
		double temp_sum_outside_cos = 0.0, temp_sum_outside_sin = 0.0;
		double temp_sum_long_cos = 0.0, temp_sum_long_sin = 0.0;
		double temp_sum_sidelong_cos = 0.0, temp_sum_sidelong_sin = 0.0;
		double temp_sum_outlong_cos = 0.0, temp_sum_outlong_sin = 0.0;

		for(int iKphi = 0; iKphi < n_localp_phi; ++iKphi)
		{
			//double point[2] = {K_T[iKT], K_phi[iKphi]};
			/*temp_sum_side_cos += (*approx_R2s).eval(point)*cos_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_side_sin += (*approx_R2s).eval(point)*sin_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_out_cos += (*approx_R2o).eval(point)*cos_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_out_sin += (*approx_R2o).eval(point)*sin_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_outside_cos += (*approx_R2os).eval(point)*cos_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_outside_sin += (*approx_R2os).eval(point)*sin_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_long_cos += (*approx_R2l).eval(point)*cos_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_long_sin += (*approx_R2l).eval(point)*sin_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_sidelong_cos += (*approx_R2sl).eval(point)*cos_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_sidelong_sin += (*approx_R2sl).eval(point)*sin_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_outlong_cos += (*approx_R2ol).eval(point)*cos_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_outlong_sin += (*approx_R2ol).eval(point)*sin_mK_phi[iKphi]*K_phi_weight[iKphi];*/
			double local_R2s = interpolate2D(SPinterp_pT, SPinterp_pphi, R2_side_GF, K_T[iKT], K_phi[iKphi], n_interp_pT_pts, n_interp_pphi_pts, interpMode, false, true);
			double local_R2o = interpolate2D(SPinterp_pT, SPinterp_pphi, R2_out_GF, K_T[iKT], K_phi[iKphi], n_interp_pT_pts, n_interp_pphi_pts, interpMode, false, true);
			double local_R2os = interpolate2D(SPinterp_pT, SPinterp_pphi, R2_outside_GF, K_T[iKT], K_phi[iKphi], n_interp_pT_pts, n_interp_pphi_pts, interpMode, false, true);
			double local_R2l = interpolate2D(SPinterp_pT, SPinterp_pphi, R2_long_GF, K_T[iKT], K_phi[iKphi], n_interp_pT_pts, n_interp_pphi_pts, interpMode, false, true);
			double local_R2sl = interpolate2D(SPinterp_pT, SPinterp_pphi, R2_sidelong_GF, K_T[iKT], K_phi[iKphi], n_interp_pT_pts, n_interp_pphi_pts, interpMode, false, true);
			double local_Rol = interpolate2D(SPinterp_pT, SPinterp_pphi, R2_outlong_GF, K_T[iKT], K_phi[iKphi], n_interp_pT_pts, n_interp_pphi_pts, interpMode, false, true);
			temp_sum_side_cos += local_R2s*cos_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_side_sin += local_R2s*sin_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_out_cos += local_R2o*cos_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_out_sin += local_R2o*sin_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_outside_cos += local_R2os*cos_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_outside_sin += local_R2os*sin_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_long_cos += local_R2l*cos_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_long_sin += local_R2l*sin_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_sidelong_cos += local_R2sl*cos_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_sidelong_sin += local_R2sl*sin_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_outlong_cos += local_Rol*cos_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_outlong_sin += local_Rol*sin_mK_phi[iKphi]*K_phi_weight[iKphi];
		}

		if (mode == 0)
		{
			R2_side_GF_C[iKT][Morder] = temp_sum_side_cos/(2.*M_PI);
			R2_side_GF_S[iKT][Morder] = temp_sum_side_sin/(2.*M_PI);
			R2_out_GF_C[iKT][Morder] = temp_sum_out_cos/(2.*M_PI);
			R2_out_GF_S[iKT][Morder] = temp_sum_out_sin/(2.*M_PI);
			R2_outside_GF_C[iKT][Morder] = temp_sum_outside_cos/(2.*M_PI);
			R2_outside_GF_S[iKT][Morder] = temp_sum_outside_sin/(2.*M_PI);
			R2_long_GF_C[iKT][Morder] = temp_sum_long_cos/(2.*M_PI);
			R2_long_GF_S[iKT][Morder] = temp_sum_long_sin/(2.*M_PI);
			R2_sidelong_GF_C[iKT][Morder] = temp_sum_sidelong_cos/(2.*M_PI);
			R2_sidelong_GF_S[iKT][Morder] = temp_sum_sidelong_sin/(2.*M_PI);
			R2_outlong_GF_C[iKT][Morder] = temp_sum_outlong_cos/(2.*M_PI);
			R2_outlong_GF_S[iKT][Morder] = temp_sum_outlong_sin/(2.*M_PI);
		}
	}

	return;
}

//performs extrapolation of running_sum of FO integrals to unity (1) by polynomial fit
double FitCF::gsl_polynomial_fit(const vector<double> &data_x, const vector<double> &data_y, const int order, double & chisq, bool verbose /* == false*/)
{
	const int n = data_x.size();
	double * in_data_x = new double [n];
	double * in_data_y = new double [n];
	gsl_vector *y, *c;
	gsl_matrix *X, *cov;
	y = gsl_vector_alloc (n);
	c = gsl_vector_alloc (order+1);
	X   = gsl_matrix_alloc (n, order+1);
	cov = gsl_matrix_alloc (order+1, order+1);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < order+1; j++)
		{
			in_data_x[i] = data_x[i];
			gsl_matrix_set (X, i, j, pow(in_data_x[i],j));
		}
		in_data_y[i] = data_y[i];
		gsl_vector_set (y, i, in_data_y[i]);
	}

	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, order+1);
	gsl_multifit_linear (X, y, c, cov, &chisq, work);
	gsl_multifit_linear_free (work);

	vector<double> vc;
	for (int i = 0; i < order+1; i++)
	{
		vc.push_back(gsl_vector_get(c,i));
		if (verbose) cerr << "In gsl_polynomial_fit(): vc[" << i << "] = " << vc[i] << endl;
	}

	gsl_vector_free (y);
	gsl_vector_free (c);
	gsl_matrix_free (X);
	gsl_matrix_free (cov);

	delete [] in_data_x;
	delete [] in_data_y;

	return ( accumulate(vc.begin(), vc.end(), 0.0) );
}



//End of file
