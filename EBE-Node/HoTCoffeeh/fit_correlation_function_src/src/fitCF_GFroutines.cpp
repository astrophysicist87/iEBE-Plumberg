#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>

#include "fitCF.h"
#include "fitCF_lib.h"
#include "CPStopwatch.h"
#include "Stopwatch.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"
#include "chebyshev.h"
#include "stats.h"

using namespace std;

bool runningOutsideCurrentFolderMode = true;

void FitCF::Get_GF_HBTradii()
{
	*global_out_stream_ptr << "--> Getting HBT radii by Gaussian fit method" << endl;

	if (FLESH_OUT_CF)
		Allocate_fleshed_out_CF();

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		*global_out_stream_ptr << "   --> Doing pT = " << SPinterp_pT[ipt] << ", pphi = " << SPinterp_pphi[ipphi] << "..." << endl;
		
		//determine whether to use fleshed out / projected CFvals
		double *** CF_for_fitting = avgCorrelation_function[ipt][ipphi];
		if (FLESH_OUT_CF)
		{
			Flesh_out_CF(ipt, ipphi);
			CF_for_fitting = fleshed_out_CF;
		}

		find_minimum_chisq_correlationfunction_full( CF_for_fitting, ipt, ipphi, FLESH_OUT_CF );
	}

	if (FLESH_OUT_CF)
		Delete_fleshed_out_CF();

	return;
}

//to save time, don't bother making grid large enough to interpolate to OSL
//this function leaves open the option of just interpolating over the qt-direction
void FitCF::Cal_correlationfunction()
{
	*global_out_stream_ptr << "Calculating the correlation function..." << endl;

	// Can't interpolate if there's only one point in qt-direction!
	if (qtnpts == 1)
		return;

	// chooses the qo, qs, ql (or qx, qy, ql) points at which to evaluate correlation function,
	// and allocates the array to hold correlation function values
	Set_correlation_function_q_pts();
	Allocate_CFvals();

	string local_name = all_particles[target_particle_id].name;
	replace_parentheses(local_name);

	string averageStem = ( nEvents > 1 ) ? "_evavg" : "";
	//string thermalOrResonanceStem = ( THERMAL_ONLY ) ? "_thermal" : "_resonance";

	ostringstream filename_stream_correlation_function_evavg;
	filename_stream_correlation_function_evavg << global_path
				<< "/correlation_function_" << local_name << averageStem << no_df_stem << ".dat";

	//now check if file exists
	bool avg_file_exists = fexists(filename_stream_correlation_function_evavg.str().c_str());

	if (avg_file_exists)
	{
		Read_in_correlationfunction_evavg(filename_stream_correlation_function_evavg.str());
	}
	else
	{
		////////////////////////////////////////////////////////////
		//reading in events...
		for (int iEvent = 0; iEvent < nEvents; ++iEvent)
		{
			*global_out_stream_ptr << "Reading in event = " << chosen_events[iEvent] << endl;

			//construct needed filename
			string resultsDirectory = "";
			if (currentfolderindex == -1 or runningOutsideCurrentFolderMode)
				resultsDirectory = "/results-" + patch::to_string(chosen_events[iEvent]);
			ostringstream filename_stream_correlation_function;
			filename_stream_correlation_function << global_path << resultsDirectory
						<< "/results/correlfunct3D_" << local_name << ".dat";

			//now check if file exists
			bool file_exists = fexists(filename_stream_correlation_function.str().c_str());
			if (file_exists)
				*global_out_stream_ptr << "    --> Note: file exists." << endl;
			else
			{
				*global_out_stream_ptr << "    --> Note: file does not exist." << endl;
				exit(8);
			}

			//now read the file in and average appropriately
			//keeps running sum for simplicity
			Read_in_correlationfunction(filename_stream_correlation_function.str());
		}

		Average_correlation_function();	//if nEvents=1, does nothing
	}
	////////////////////////////////////////////////////////////////////////

	return;
}



void FitCF::Average_correlation_function()
{
	*global_out_stream_ptr << "Ensemble-averaging EbE correlation functions..." << endl;

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		avgSpectra[target_particle_id][ipt][ipphi] /= double(nEvents);

		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int iqz = 0; iqz < qznpts; ++iqz)
		{
			avgThermalCFvals[ipt][ipphi][iqx][iqy][iqz] /= double(nEvents);
			avgCrosstermCFvals[ipt][ipphi][iqx][iqy][iqz] /= double(nEvents);
			avgResonancesCFvals[ipt][ipphi][iqx][iqy][iqz] /= double(nEvents);
			avgCorrelation_function_Numerator[ipt][ipphi][iqx][iqy][iqz] /= double(nEvents);
			avgCorrelation_function_Denominator[ipt][ipphi][iqx][iqy][iqz] /= double(nEvents);

			avgThermalCFvals[ipt][ipphi][iqx][iqy][iqz]
				/= avgCorrelation_function_Denominator[ipt][ipphi][iqx][iqy][iqz];
			avgCrosstermCFvals[ipt][ipphi][iqx][iqy][iqz]
				/= avgCorrelation_function_Denominator[ipt][ipphi][iqx][iqy][iqz];
			avgResonancesCFvals[ipt][ipphi][iqx][iqy][iqz]
				/= avgCorrelation_function_Denominator[ipt][ipphi][iqx][iqy][iqz];

			avgCorrelation_function[ipt][ipphi][iqx][iqy][iqz]
				= 1.0 + avgCorrelation_function_Numerator[ipt][ipphi][iqx][iqy][iqz]
					/ avgCorrelation_function_Denominator[ipt][ipphi][iqx][iqy][iqz];
		}
	}

	*global_out_stream_ptr << "Finished ensemble-averaging EbE correlation functions." << endl;

	return;
}

//******************************************************************
// Routines for refining correlation function grid via interpolation
//******************************************************************

void FitCF::Flesh_out_CF(int ipt, int ipphi, double sample_scale /*==1.0*/)
{
	//declare needed quantities here
	double qxmin = 0.9999*sample_scale*qx_pts[0], qxmax = 0.9999*sample_scale*qx_pts[qxnpts-1];
	double qymin = 0.9999*sample_scale*qy_pts[0], qymax = 0.9999*sample_scale*qy_pts[qynpts-1];
	double qzmin = 0.9999*sample_scale*qz_pts[0], qzmax = 0.9999*sample_scale*qz_pts[qznpts-1];

	double new_Del_qx = (qxmax - qxmin)/(double(new_nqpts-1)+1.e-100);
	double new_Del_qy = (qymax - qymin)/(double(new_nqpts-1)+1.e-100);
	double new_Del_qz = (qzmax - qzmin)/(double(new_nqpts-1)+1.e-100);

	double *** current_C_slice = CFvals[ipt][ipphi];
	
	for (int iqx = 0; iqx < new_nqpts; ++iqx)
	for (int iqy = 0; iqy < new_nqpts; ++iqy)
	for (int iqz = 0; iqz < new_nqpts; ++iqz)
	{
		double qx0 = qxmin + double(iqx) * new_Del_qx;
		double qy0 = qymin + double(iqy) * new_Del_qy;
		double qz0 = qzmin + double(iqz) * new_Del_qz;

		qx_fleshed_out_pts[iqx] = qx0;
		qy_fleshed_out_pts[iqy] = qy0;
		qz_fleshed_out_pts[iqz] = qz0;

		fleshed_out_thermal[iqx][iqy][iqz] = interpolate_CF(avgThermalCFvals[ipt][ipphi], qx0, qy0, qz0, ipt, 0);
		fleshed_out_crossterm[iqx][iqy][iqz] = interpolate_CF(avgCrosstermCFvals[ipt][ipphi], qx0, qy0, qz0, ipt, 1);
		fleshed_out_resonances[iqx][iqy][iqz] = interpolate_CF(avgResonancesCFvals[ipt][ipphi], qx0, qy0, qz0, ipt, 2);
		fleshed_out_CF[iqx][iqy][iqz] = 1.0 + fleshed_out_thermal[iqx][iqy][iqz] + fleshed_out_crossterm[iqx][iqy][iqz] + fleshed_out_resonances[iqx][iqy][iqz];
	}

	return;
}

double FitCF::interpolate_CF(double *** current_C_slice, double qx0, double qy0, double qz0, int ipt, int thermal_or_resonances)
{
	//int thermal_or_resonances - 	0: interpolate thermal CF (assuming basically Gaussian)
	//								1: interpolate resonance contributions (linear near origin, exponential further out)

	int iqx0_loc = binarySearch(qx_pts, qxnpts, qx0, true, true);
	int iqy0_loc = binarySearch(qy_pts, qynpts, qy0, true, true);
	int iqz0_loc = binarySearch(qz_pts, qznpts, qz0, true, true);

	if (iqx0_loc == -1 || iqy0_loc == -1 || iqz0_loc == -1)
	{
		cerr << "Interpolation failed: exiting!" << endl;
		exit(1);
	}

	double fx0y0z0 = current_C_slice[iqx0_loc][iqy0_loc][iqz0_loc];
	double fx0y0z1 = current_C_slice[iqx0_loc][iqy0_loc][iqz0_loc+1];
	double fx0y1z0 = current_C_slice[iqx0_loc][iqy0_loc+1][iqz0_loc];
	double fx0y1z1 = current_C_slice[iqx0_loc][iqy0_loc+1][iqz0_loc+1];
	double fx1y0z0 = current_C_slice[iqx0_loc+1][iqy0_loc][iqz0_loc];
	double fx1y0z1 = current_C_slice[iqx0_loc+1][iqy0_loc][iqz0_loc+1];
	double fx1y1z0 = current_C_slice[iqx0_loc+1][iqy0_loc+1][iqz0_loc];
	double fx1y1z1 = current_C_slice[iqx0_loc+1][iqy0_loc+1][iqz0_loc+1];

	double fx0 = current_C_slice[iqx0_loc][0][0];
	double fx1 = current_C_slice[iqx0_loc+1][0][0];

	bool use_log_current_C_slice = true;
	double *** log_current_C_slice = new double ** [qxnpts];
	//double *** sgn_current_C_slice = new double ** [qxnpts];
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	{
		log_current_C_slice[iqx] = new double * [qynpts];
		//sgn_current_C_slice[iqx] = new double * [qynpts];
		for (int iqy = 0; iqy < qynpts; ++iqy)
		{
			log_current_C_slice[iqx][iqy] = new double [qznpts];
			//sgn_current_C_slice[iqx][iqy] = new double [qznpts];
			for (int iqz = 0; iqz < qznpts; ++iqz)
			{
				log_current_C_slice[iqx][iqy][iqz] = 0.0;
				//sgn_current_C_slice[iqx][iqy][iqz] = 0.0;
			}
		}
	}


	if (use_log_current_C_slice)
	{
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int iqz = 0; iqz < qznpts; ++iqz)
		{
			log_current_C_slice[iqx][iqy][iqz] = log(current_C_slice[iqx][iqy][iqz]+1.0);	//make it positive to avoid NaNs
			//sgn_current_C_slice[iqx][iqy][iqz] = sgn(current_C_slice[iqx][iqy][iqz]);
		}
	}


	//interpolate here...
	double qx_0 = qx_pts[iqx0_loc];
	double qx_1 = qx_pts[iqx0_loc+1];
	double qy_0 = qy_pts[iqy0_loc];
	double qy_1 = qy_pts[iqy0_loc+1];
	double qz_0 = qz_pts[iqz0_loc];
	double qz_1 = qz_pts[iqz0_loc+1];

	double fxiyizi = 0.0;
	if (thermal_or_resonances == 0)
	{
		double sqx02 = sgn(qx0)*qx0*qx0;
		double sqy02 = sgn(qy0)*qy0*qy0;
		double sqz02 = sgn(qz0)*qz0*qz0;
		double sqx_02 = sgn(qx_0)*qx_0*qx_0;
		double sqy_02 = sgn(qy_0)*qy_0*qy_0;
		double sqz_02 = sgn(qz_0)*qz_0*qz_0;
		double sqx_12 = sgn(qx_1)*qx_1*qx_1;
		double sqy_12 = sgn(qy_1)*qy_1*qy_1;
		double sqz_12 = sgn(qz_1)*qz_1*qz_1;

		//interpolate over qz-points first
		///////////////////////
		//interpolate over each pair of qz-points
		double fx0y0zi = interpolate_qi(sqz02, sqz_02, sqz_12, fx0y0z0, fx0y0z1, false);
		double fx0y1zi = interpolate_qi(sqz02, sqz_02, sqz_12, fx0y1z0, fx0y1z1, false);    
		double fx1y0zi = interpolate_qi(sqz02, sqz_02, sqz_12, fx1y0z0, fx1y0z1, false);
		double fx1y1zi = interpolate_qi(sqz02, sqz_02, sqz_12, fx1y1z0, fx1y1z1, false);
		///////////////////////
	
		//interpolate over qy-points next
		double fx0yizi = interpolate_qi(sqy02, sqy_02, sqy_12, fx0y0zi, fx0y1zi, false);
		double fx1yizi = interpolate_qi(sqy02, sqy_02, sqy_12, fx1y0zi, fx1y1zi, false);
	
		//finally, interpolate over qx-points
		fxiyizi = interpolate_qi(sqx02, sqx_02, sqx_12, fx0yizi, fx1yizi, false);
		//fxiyizi = interpolate_qi(sqx02, sqx_02, sqx_12, fx0, fx1, false);
	}
	else
	{
		if (use_log_current_C_slice)
		{
			//separate into "thirds" in each direction == 27 "cubes" in q-space
			//use linear interpolation in middle "cube"; cubic everywhere else (seems to work best)
			bool in_central_cube = ( abs(qx0) <= qx_pts[qxnpts-1]/3.0 )
									and ( abs(qy0) <= qy_pts[qynpts-1]/3.0 )
									and ( abs(qz0) <= qz_pts[qznpts-1]/3.0 );
			fxiyizi = exp(interpolate3D(qx_pts, qy_pts, qz_pts, log_current_C_slice, qx0, qy0, qz0, qxnpts, qynpts, qznpts, 1 - int(in_central_cube), true)) - 1.0;
		}
		else
		{
			fxiyizi = interpolate3D(qx_pts, qy_pts, qz_pts, current_C_slice, qx0, qy0, qz0, qxnpts, qynpts, qznpts, 1, true);
		}
	}

	for (int iqx = 0; iqx < qxnpts; ++iqx)
	{
		for (int iqy = 0; iqy < qynpts; ++iqy)
			delete [] log_current_C_slice[iqx][iqy];
		delete [] log_current_C_slice[iqx];
	}
	delete [] log_current_C_slice;


	return (fxiyizi);
}

double FitCF::interpolate_qi(double q0, double qi0, double qi1, double f1, double f2, bool use_linear)
{
	double if1 = f1;
	double if2 = f2;
	bool use_log = (!use_linear) && (f1 > 0.0) && (f2 > 0.0);
	if (use_log)
	{
		if1 = log(f1);
		if2 = log(f2);
	}
	double tmp_result = lin_int(q0 - qi0, 1./(qi1-qi0), if1, if2);
	if (use_log)
		tmp_result = exp(tmp_result);

	return (tmp_result);
}

//**************************************************************
// Gaussian fit routine below
//**************************************************************
void FitCF::find_minimum_chisq_correlationfunction_full(
		double *** Correl_3D, int ipt, int ipphi, bool fleshing_out_CF /*== true*/)
{
	double * q1pts = qx_pts;
	double * q2pts = qy_pts;
	double * q3pts = qz_pts;
	if (fleshing_out_CF)
	{
		q1npts = new_nqxpts;
		q2npts = new_nqypts;
		q3npts = new_nqzpts;
		q1pts = qx_fleshed_out_pts;
		q2pts = qy_fleshed_out_pts;
		q3pts = qz_fleshed_out_pts;
	}
	const size_t data_length = q1npts*q2npts*q3npts;  // # of points

    double lambda, R_o, R_s, R_l, R_os;
    int dim = 5;
    int s_gsl;

    double *V = new double [dim];
    double *qweight = new double [dim];
    double **T = new double* [dim];
    for(int i = 0; i < dim; i++)
    {
        V[i] = 0.0;
        T[i] = new double [dim];
        for(int j = 0; j < dim; j++)
            T[i][j] = 0.0;
    }

    gsl_matrix * T_gsl = gsl_matrix_alloc (dim, dim);
    gsl_matrix * T_inverse_gsl = gsl_matrix_alloc (dim, dim);
    gsl_permutation * perm = gsl_permutation_alloc (dim);

	//double ckp = cos_SP_pphi[ipphi], skp = sin_SP_pphi[ipphi];
	double ckp = cos(SPinterp_pphi[ipphi]), skp = sin(SPinterp_pphi[ipphi]);
	double CF_err = 1.e-3;
	for (int i = 0; i < q1npts; i++)
	for (int j = 0; j < q2npts; j++)
	for (int k = 0; k < q3npts; k++)
    {
        double q_out_local = q1pts[i] * ckp + q2pts[j] * skp;
        double q_side_local = -q1pts[i] * skp + q2pts[j] * ckp;
        double q_long_local = q3pts[k];
        double correl_local = Correl_3D[i][j][k]-1;
		//*global_out_stream_ptr << "\t\t" << q1pts[i] << "   " << q2pts[i] << "   " << q3pts[i] << "   " << q_out_local << "   " << q_side_local << "   " << q_long_local << "   " << correl_local << endl;
        if(correl_local < 1e-15) continue;
		//if (i==(q1npts-1)/2 && j==(q2npts-1)/2 && k==(q3npts-1)/2)
		//	Correlfun3D_data.sigma[idx] = 1.e10;	//ignore central point
        double sigma_k_prime = CF_err/correl_local;
            
        double inv_sigma_k_prime_sq = 1./(sigma_k_prime*sigma_k_prime);
        double log_correl_over_sigma_sq = log(correl_local)*inv_sigma_k_prime_sq;

        qweight[0] = - 1.0;
        qweight[1] = q_out_local*q_out_local;
        qweight[2] = q_side_local*q_side_local;
        qweight[3] = q_long_local*q_long_local;
        qweight[4] = q_out_local*q_side_local;

        for(int ij = 0; ij < dim; ij++)
        {
            V[ij] += qweight[ij]*log_correl_over_sigma_sq;
            T[0][ij] += qweight[ij]*inv_sigma_k_prime_sq;
        }

        for(int ij = 1; ij < dim; ij++)
            T[ij][0] = T[0][ij];
            

        for(int ij = 1; ij < dim; ij++)
        {
            for(int lm = 1; lm < dim; lm++)
                T[ij][lm] += -qweight[ij]*qweight[lm]*inv_sigma_k_prime_sq;
        }
    }
    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
            gsl_matrix_set(T_gsl, i, j, T[i][j]);

    // Make LU decomposition of matrix T_gsl
    gsl_linalg_LU_decomp (T_gsl, perm, &s_gsl);
    // Invert the matrix m
    gsl_linalg_LU_invert (T_gsl, perm, T_inverse_gsl);

    double **T_inverse = new double* [dim];
    for(int i = 0; i < dim; i++)
    {
        T_inverse[i] = new double [dim];
        for(int j = 0; j < dim; j++)
            T_inverse[i][j] = gsl_matrix_get(T_inverse_gsl, i, j);
    }
    double *results = new double [dim];
    for(int i = 0; i < dim; i++)
    {
        results[i] = 0.0;
        for(int j = 0; j < dim; j++)
            results[i] += T_inverse[i][j]*V[j];
    }

	lambda_Correl[ipt][ipphi] = exp(results[0]);
	lambda_Correl_err[ipt][ipphi] = 0.0;
	R2_out_GF[ipt][ipphi] = results[1]*hbarC*hbarC;
	R2_side_GF[ipt][ipphi] = results[2]*hbarC*hbarC;
	R2_long_GF[ipt][ipphi] = results[3]*hbarC*hbarC;
	R2_outside_GF[ipt][ipphi] = results[4]*hbarC*hbarC;
	R2_out_err[ipt][ipphi] = 0.0;
	R2_side_err[ipt][ipphi] = 0.0;
	R2_long_err[ipt][ipphi] = 0.0;
	R2_outside_err[ipt][ipphi] = 0.0;


    double chi_sq = 0.0;
	for (int i = 0; i < q1npts; i++)
	for (int j = 0; j < q2npts; j++)
	for (int k = 0; k < q3npts; k++)
    {
        double q_out_local = q1pts[i] * ckp + q2pts[j] * skp;
        double q_side_local = -q1pts[i] * skp + q2pts[j] * ckp;
        double q_long_local = q3pts[k];
        double correl_local = Correl_3D[i][j][k]-1;
        if(correl_local < 1e-15) continue;
        double sigma_k_prime = CF_err/correl_local;

        chi_sq += pow((log(correl_local) - results[0] 
                       + results[1]*q_out_local*q_out_local 
                       + results[2]*q_side_local*q_side_local
                       + results[3]*q_long_local*q_long_local
                       + results[4]*q_out_local*q_side_local), 2)
                  /sigma_k_prime/sigma_k_prime;
    }
    //cout << "chi_sq/d.o.f = " << chi_sq/(qnpts - dim) << endl;
    //chi_sq_per_dof = chi_sq/(qnpts - dim);

    // clean up
    gsl_matrix_free (T_gsl);
    gsl_matrix_free (T_inverse_gsl);
    gsl_permutation_free (perm);

    delete [] qweight;
    delete [] V;
    for(int i = 0; i < dim; i++)
    {
        delete [] T[i];
        delete [] T_inverse[i];
    }
    delete [] T;
    delete [] T_inverse;
    delete [] results;
}


//End of file
