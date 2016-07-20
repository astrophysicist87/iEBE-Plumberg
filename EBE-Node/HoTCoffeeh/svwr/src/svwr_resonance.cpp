#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>
#include<time.h>

#include<gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>            // gsl random number generators
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting

#include "svwr.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"

using namespace std;

const double PTCHANGE = 1.0;	//GeV
const bool USE_PTCHANGE = true;

double SourceVariances::get_Q()
{
	double smin = (m2+m3)*(m2+m3);
	double smax = (Mres-mass)*(Mres-mass);
	double sum = 0.;
	
	for (int is = 0; is < n_s_pts; ++is)
	{
		double sp = s_pts[is];
		double f1 = (Mres+mass)*(Mres+mass) - sp;
		double f2 = smax - sp;
		double f3 = smin - sp;
		double f4 = (m2-m3)*(m2-m3) - sp;
		sum += s_wts[is]*sqrt(f1*f2*f3*f4)/(sp+1.e-15);
	}

	return sum;
}

double SourceVariances::g(double s)
{
	double gs_pstar_loc = sqrt( ((Mres+mass)*(Mres+mass) - s)*((Mres-mass)*(Mres-mass) - s) )/(2.0*Mres);
	double g_res = br/(4.*M_PI*gs_pstar_loc);
	if (n_body == 3 || n_body == 4)		//both set up to work the same way
	{
		double pre_f = (Mres*br)/(2.*M_PI*s);
		double num = sqrt( (s - (m2+m3)*(m2+m3)) * (s - (m2-m3)*(m2-m3)) );
		double den = Qfunc;
		g_res = pre_f * num / den;
	}

	return g_res;
}

void SourceVariances::get_rapidity_dependence(double * rap_indep_vector, double * rap_dep_vector, double rap_val)
{
	//assumes currently calculating all 15 SVs
	//rap_indep_vector - SV vector evaluated at momentum rapidity y_p = 0
	//rap_dep_vector - SV vector evaluated at momentum rapidity y_p = rap_val
	//rap_val - value of rapidity to evaluate at
	
	//for SVs which don't change, just copy them over
	//[{1}_r]_{r-->\pi}
	rap_dep_vector[0] = rap_indep_vector[0];									//1
	//if (INCLUDE_SOURCE_VARIANCES)
	//{
		double ch_rap_val = cosh(rap_val);
		double sh_rap_val = sinh(rap_val);
	
		//[{x}_r]_{r-->\pi}
		rap_dep_vector[1] = rap_indep_vector[1];
		//[{x2}_r]_{r-->\pi}
		rap_dep_vector[2] = rap_indep_vector[2];
		//[{y}_r]_{r-->\pi}
		rap_dep_vector[3] = rap_indep_vector[3];
		//[{y2}_r]_{r-->\pi}
		rap_dep_vector[4] = rap_indep_vector[4];
		//[{z}_r]_{r-->\pi}
		rap_dep_vector[5] = ch_rap_val * rap_indep_vector[5]
					+ sh_rap_val * rap_indep_vector[7];
		//[{z2}_r]_{r-->\pi}
		rap_dep_vector[6] = ch_rap_val * ch_rap_val * rap_indep_vector[6]
					+ 2. * ch_rap_val * sh_rap_val * rap_indep_vector[14]
					+ sh_rap_val * sh_rap_val * rap_indep_vector[8];
		//[{t}_r]_{r-->\pi}
		rap_dep_vector[7] = ch_rap_val * rap_indep_vector[7]
					+ sh_rap_val * rap_indep_vector[5];
		//[{t2}_r]_{r-->\pi}
		rap_dep_vector[8] = ch_rap_val * ch_rap_val * rap_indep_vector[8]
					+ 2. * ch_rap_val * sh_rap_val * rap_indep_vector[14]
					+ sh_rap_val * sh_rap_val * rap_indep_vector[6];
		//[{xy}_r]_{r-->\pi}
		rap_dep_vector[9] = rap_indep_vector[9];
		//[{xz}_r]_{r-->\pi}
		rap_dep_vector[10] = ch_rap_val * rap_indep_vector[10]
					+ sh_rap_val * rap_indep_vector[11];
		//[{xt}_r]_{r-->\pi}
		rap_dep_vector[11] = ch_rap_val * rap_indep_vector[11]
					+ sh_rap_val * rap_indep_vector[10];
		//[{yz}_r]_{r-->\pi}
		rap_dep_vector[12] = ch_rap_val * rap_indep_vector[12]
					+ sh_rap_val * rap_indep_vector[13];
		//[{yt}_r]_{r-->\pi}
		rap_dep_vector[13] = ch_rap_val * rap_indep_vector[13]
					+ sh_rap_val * rap_indep_vector[12];
		//[{zt}_r]_{r-->\pi}
		rap_dep_vector[14] = (ch_rap_val * ch_rap_val + sh_rap_val * sh_rap_val) * rap_indep_vector[14]
					+ sh_rap_val * ch_rap_val * (rap_indep_vector[6] + rap_indep_vector[8]);
	//}
	
	return;
}

void SourceVariances::combine_sourcevariances(double * output, double * input, double * alpha_vec)
{
	//0, 1, 2, 3 --> t, x, y, z
	//[{1}_r]_{r-->\pi}
	output[0] += input[0];
	//if (INCLUDE_SOURCE_VARIANCES)
	//{
		double ax = alpha_vec[1], ay = alpha_vec[2], az = alpha_vec[3], at = alpha_vec[0];
		//[{x}_r]_{r-->\pi}
		output[1] += input[1] + ax*input[0];
		//[{x2}_r]_{r-->\pi}
		output[2] += input[2] + 2.*ax*input[1] + 2.*ax*ax*input[0];
		//[{y}_r]_{r-->\pi}
		output[3] += input[3] + ay*input[0];
		//[{y2}_r]_{r-->\pi}
		output[4] += input[4] + 2.*ay*input[3] + 2.*ay*ay*input[0];
		//[{z}_r]_{r-->\pi}
		output[5] += input[5] + az*input[0];
		//[{z2}_r]_{r-->\pi}
		output[6] += input[6] + 2.*az*input[5] + 2.*az*az*input[0];
		//[{t}_r]_{r-->\pi}
		output[7] += input[7] + at*input[0];
		//[{t2}_r]_{r-->\pi}
		output[8] += input[8] + 2.*at*input[7] + 2.*at*at*input[0];
		//[{xy}_r]_{r-->\pi}
		output[9] += input[9] + ax*input[3] + ay*input[1] + 2.*ax*ay*input[0];
		//[{xz}_r]_{r-->\pi}
		output[10] += input[10] + ax*input[5] + az*input[1] + 2.*ax*az*input[0];
		//[{xt}_r]_{r-->\pi}
		output[11] += input[11] + ax*input[7] + at*input[1] + 2.*ax*at*input[0];
		//[{yz}_r]_{r-->\pi}
		output[12] += input[12] + ay*input[5] + az*input[3] + 2.*ay*az*input[0];
		//[{yt}_r]_{r-->\pi}
		output[13] += input[13] + ay*input[7] + at*input[3] + 2.*ay*at*input[0];
		//[{zt}_r]_{r-->\pi}
		output[14] += input[14] + az*input[7] + at*input[5] + 2.*az*at*input[0];
	//}
	
	return;
}

void SourceVariances::Do_resonance_integrals(int parent_resonance_particle_id, int daughter_particle_id, int decay_channel)
{
	time_t rawtime;
  	struct tm * timeinfo;
	double * ssum_vec = new double [n_weighting_functions];
	double * vsum_vec = new double [n_weighting_functions];
	double * zetasum_vec = new double [n_weighting_functions];
	double * Csum_vec = new double [n_weighting_functions];
	double * rap_indep_y_of_r = new double [n_weighting_functions];
	double * y_of_r = new double [n_weighting_functions];
	double * alphavec = new double [4];
	set_to_zero(ssum_vec, n_weighting_functions);
	set_to_zero(vsum_vec, n_weighting_functions);
	set_to_zero(zetasum_vec, n_weighting_functions);
	set_to_zero(Csum_vec, n_weighting_functions);
	set_to_zero(rap_indep_y_of_r, n_weighting_functions);
	set_to_zero(y_of_r, n_weighting_functions);

	res_sign_info = sign_of_dN_dypTdpTdphi_moments[parent_resonance_particle_id];
	res_log_info = ln_dN_dypTdpTdphi_moments[parent_resonance_particle_id];
	res_moments_info = dN_dypTdpTdphi_moments[parent_resonance_particle_id];

	n_body = current_reso_nbody;

	if (n_body == 2)
	{
		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
//cerr << "Entering " << n_body << "   " << ipt << "   " << ipphi << endl;
			double local_pT = SPinterp_pT[ipt];
			double local_pphi = SPinterp_pphi[ipphi];
			set_to_zero(ssum_vec, n_weighting_functions);
			set_to_zero(vsum_vec, n_weighting_functions);
			set_to_zero(zetasum_vec, n_weighting_functions);
			set_to_zero(Csum_vec, n_weighting_functions);
			set_to_zero(rap_indep_y_of_r, n_weighting_functions);
			set_to_zero(y_of_r, n_weighting_functions);
			Load_decay_channel_info(decay_channel, local_pT, local_pphi);	// set decay channel information

			//then g(s) is delta-function, skip s-integration entirely
			//double s_loc = m2*m2;
			for (int iv = 0; iv < n_v_pts; ++iv)
			{
				double zetasum = 0.0;
				time (&rawtime);
				timeinfo = localtime (&rawtime);
				set_to_zero(zetasum_vec, n_weighting_functions);
				for (int izeta = 0; izeta < n_zeta_pts; ++izeta)
				{
					alphavec = VEC_n2_alpha[iv][izeta];
					set_to_zero(Csum_vec, n_weighting_functions);
					double PKT = VEC_n2_PT[iv][izeta];
					double PKY = VEC_n2_P_Y[iv];
					double PKphi = VEC_n2_PPhi_tilde[iv][izeta];
					for (int tempidx = 1; tempidx <= 2; ++tempidx)
					{
						if (tempidx != 1)
						{
							//Phi only changes sign, does NOT get shifted by pi!
							PKphi = VEC_n2_PPhi_tildeFLIP[iv][izeta];		//also takes Pp --> Pm
							alphavec = VEC_n2_alpha_m[iv][izeta];
						}
						Edndp3(PKT, PKphi, rap_indep_y_of_r);
						get_rapidity_dependence(rap_indep_y_of_r, y_of_r, PKY);
						combine_sourcevariances(Csum_vec, y_of_r, alphavec);
					}																					// end of tempidx sum
					for (int iweight = 0; iweight < n_weighting_functions; ++iweight)
						zetasum_vec[iweight] += VEC_n2_zeta_factor[iv][izeta]*Csum_vec[iweight];
				}																						// end of zeta sum
				for (int iweight = 0; iweight < n_weighting_functions; ++iweight)
					vsum_vec[iweight] += VEC_n2_v_factor[iv]*zetasum_vec[iweight];
			}																							// end of v sum
			for (int iweight = 0; iweight < n_weighting_functions; ++iweight)
				ssum_vec[iweight] += Mres*VEC_n2_s_factor*vsum_vec[iweight];

			//update all gridpoints for daughter moments
			for (int iweight = 0; iweight < n_weighting_functions; ++iweight)
			{
				dN_dypTdpTdphi_moments[daughter_particle_id][iweight][ipt][ipphi] += ssum_vec[iweight];
				double temp = dN_dypTdpTdphi_moments[daughter_particle_id][iweight][ipt][ipphi];
				ln_dN_dypTdpTdphi_moments[daughter_particle_id][iweight][ipt][ipphi] = log(abs(temp) + 1.e-100);
				sign_of_dN_dypTdpTdphi_moments[daughter_particle_id][iweight][ipt][ipphi] = sgn(temp);
			}
	
			if (isnan(dN_dypTdpTdphi_moments[daughter_particle_id][0][ipt][ipphi]))
			{
				*global_out_stream_ptr << "ERROR: NaNs encountered!" << endl
										<< "dN_dypTdpTdphi_moments[daughter_particle_id][0][" << ipt << "][" << ipphi << "] = "
										<< setw(8) << setprecision(15) << dN_dypTdpTdphi_moments[daughter_particle_id][0][ipt][ipphi] << endl
										<< "  --> pt = " << local_pT << std::endl
										<< "  --> pphi = " << local_pphi << std::endl
										<< "daughter_particle_id = " << daughter_particle_id << endl
										<< "parent_resonance_particle_id = " << parent_resonance_particle_id << endl
										<< "  --> Qfunc = " << Qfunc << endl
										<< "  --> n_body = " << n_body << endl
										<< "  --> gRES = " << gRES << endl
										<< "  --> Mres = " << Mres << endl
										<< "  --> mass = " << mass << endl
										<< "  --> Gamma = " << Gamma << endl
										<< "  --> br = " << br << endl
										<< "  --> m2 = " << m2 << endl
										<< "  --> m3 = " << m3 << endl << endl;
				exit(1);
			}
		}																								// end of pT, pphi loops
	}																									// end of nbody == 2
	else
	{
		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
//cerr << "Entering " << n_body << "   " << ipt << "   " << ipphi << endl;
			double local_pT = SPinterp_pT[ipt];
			double local_pphi = SPinterp_pphi[ipphi];
			set_to_zero(ssum_vec, n_weighting_functions);
			set_to_zero(vsum_vec, n_weighting_functions);
			set_to_zero(zetasum_vec, n_weighting_functions);
			set_to_zero(Csum_vec, n_weighting_functions);
			set_to_zero(rap_indep_y_of_r, n_weighting_functions);
			set_to_zero(y_of_r, n_weighting_functions);
			Load_decay_channel_info(decay_channel, local_pT, local_pphi);	// set decay channel information

			for (int is = 0; is < n_s_pts; ++is)
			{
				double vsum = 0.0;
 		  		set_to_zero(vsum_vec, n_weighting_functions);
				for (int iv = 0; iv < n_v_pts; ++iv)
				{
					set_to_zero(zetasum_vec, n_weighting_functions);
					for (int izeta = 0; izeta < n_zeta_pts; ++izeta)
					{
						alphavec = VEC_alpha[is][iv][izeta];
						set_to_zero(Csum_vec, n_weighting_functions);
						double PKT = VEC_PT[is][iv][izeta];
						double PKY = VEC_P_Y[is][iv];
						double PKphi = VEC_PPhi_tilde[is][iv][izeta];
						for (int tempidx = 1; tempidx <= 2; ++tempidx)
						{
							if (tempidx != 1)
							{
								PKphi = VEC_PPhi_tildeFLIP[is][iv][izeta];		//also takes Pp --> Pm
								alphavec = VEC_alpha_m[is][iv][izeta];
							}
							Edndp3(PKT, PKphi, rap_indep_y_of_r);
							get_rapidity_dependence(rap_indep_y_of_r, y_of_r, PKY);
							//now compute appropriate linear combinations
							combine_sourcevariances(Csum_vec, y_of_r, alphavec);
						}																					// end of tempidx sum
						for (int iweight = 0; iweight < n_weighting_functions; ++iweight)
							zetasum_vec[iweight] += VEC_zeta_factor[is][iv][izeta]*Csum_vec[iweight];
					}																						// end of zeta sum
					for (int iweight = 0; iweight < n_weighting_functions; ++iweight)
						vsum_vec[iweight] += VEC_v_factor[is][iv]*zetasum_vec[iweight];
				}																							// end of v sum
				for (int iweight = 0; iweight < n_weighting_functions; ++iweight)
					ssum_vec[iweight] += Mres*VEC_s_factor[is]*vsum_vec[iweight];
			}																								// end of s sum
			//update all gridpoints for daughter moments
			for (int iweight = 0; iweight < n_weighting_functions; ++iweight)
			{
				dN_dypTdpTdphi_moments[daughter_particle_id][iweight][ipt][ipphi] += ssum_vec[iweight];
				double temp = dN_dypTdpTdphi_moments[daughter_particle_id][iweight][ipt][ipphi];
				ln_dN_dypTdpTdphi_moments[daughter_particle_id][iweight][ipt][ipphi] = log(abs(temp) + 1.e-100);
				sign_of_dN_dypTdpTdphi_moments[daughter_particle_id][iweight][ipt][ipphi] = sgn(temp);
			}
	
			if (isnan(dN_dypTdpTdphi_moments[daughter_particle_id][0][ipt][ipphi]))
			{
				*global_out_stream_ptr << "ERROR: NaNs encountered!" << endl
										<< "dN_dypTdpTdphi_moments[daughter_particle_id][0][" << ipt << "][" << ipphi << "] = "
										<< setw(8) << setprecision(15) << dN_dypTdpTdphi_moments[daughter_particle_id][0][ipt][ipphi] << endl
										<< "  --> pt = " << local_pT << std::endl
										<< "  --> pphi = " << local_pphi << std::endl
										<< "daughter_particle_id = " << daughter_particle_id << endl
										<< "parent_resonance_particle_id = " << parent_resonance_particle_id << endl
										<< "  --> Qfunc = " << Qfunc << endl
										<< "  --> n_body = " << n_body << endl
										<< "  --> gRES = " << gRES << endl
										<< "  --> Mres = " << Mres << endl
										<< "  --> mass = " << mass << endl
										<< "  --> Gamma = " << Gamma << endl
										<< "  --> br = " << br << endl
										<< "  --> m2 = " << m2 << endl
										<< "  --> m3 = " << m3 << endl << endl;
				exit(1);
			}
		}																									// end of pT, pphi loops
	}																										// end of nbody == 3



	//clean up
	delete [] ssum_vec;
	delete [] vsum_vec;
	delete [] zetasum_vec;
	delete [] Csum_vec;
	delete [] rap_indep_y_of_r;
	delete [] y_of_r;

	return;
}

inline void SourceVariances::set_to_zero(double * array, size_t arraylength)
{
	for (size_t arrayidx=0; arrayidx<arraylength; ++arrayidx) array[arrayidx] = 0.0;
}

inline double SourceVariances::lin_int(double x_m_x1, double one_by_x2_m_x1, double f1, double f2)
{
	return ( f1 + (f2 - f1) * x_m_x1 * one_by_x2_m_x1 );
}

void SourceVariances::Edndp3(double ptr, double phir, double * results)
{//move loop over n_weighting_functions inside this function
	double phi0, phi1;
	double f1, f2;

	int npphi_max = n_interp_pphi_pts - 1;
	int npT_max = n_interp_pT_pts - 1;

	// locate pT interval
	int npt = 1;
	while ((ptr > SPinterp_pT[npt]) &&
			(npt < npT_max)) ++npt;
	double pT0 = SPinterp_pT[npt-1];
	double pT1 = SPinterp_pT[npt];

	// locate pphi interval
	int nphi = 1, nphim1 = 0;
	if(phir < SPinterp_pphi[0])			//if angle is less than minimum angle grid point
	{
		phi0 = SPinterp_pphi[npphi_max] - 2. * M_PI;
		phi1 = SPinterp_pphi[0];
		nphi = 0;
		nphim1 = npphi_max;
	}
	else if(phir > SPinterp_pphi[npphi_max])	//if angle is greater than maximum angle grid point
	{
		phi0 = SPinterp_pphi[npphi_max];
		phi1 = SPinterp_pphi[0] + 2. * M_PI;
		nphi = 0;
		nphim1 = npphi_max;
	}
	else						//if angle is within grid range
	{
		while ((phir > SPinterp_pphi[nphi]) &&
				(nphi < npphi_max)) ++nphi;
		nphim1 = nphi - 1;
		phi0 = SPinterp_pphi[nphim1];
		phi1 = SPinterp_pphi[nphi];
	}

	if (pT0==pT1 || phi0==phi1)
	{
		cerr << "ERROR in Edndp3(): pT and/or pphi values equal!" << endl;
		exit(1);
	}

	double one_by_pTdiff = 1./(pT1 - pT0), one_by_pphidiff = 1./(phi1 - phi0);

	for (int wfi = 0; wfi < n_weighting_functions; ++wfi)
	{
		double ** temp_res_sign_info = res_sign_info[wfi];
		double ** temp_res_log_info = res_log_info[wfi];
		double ** temp_res_moments_info = res_moments_info[wfi];
		
		// interpolate over pT values first
		if(ptr > PTCHANGE && USE_PTCHANGE)				// if pT interpolation point is larger than PTCHANGE (currently 1.0 GeV)
		{
			double sign_of_f11 = temp_res_sign_info[npt-1][nphim1];
			double sign_of_f12 = temp_res_sign_info[npt-1][nphi];
			double sign_of_f21 = temp_res_sign_info[npt][nphim1];
			double sign_of_f22 = temp_res_sign_info[npt][nphi];
	
			//*******************************************************************************************************************
			// set f1 first
			//*******************************************************************************************************************
			// if using extrapolation and spectra at pT1 has larger magnitude than at pT0, just return zero
			if (ptr > pT1 && ( temp_res_log_info[npt][nphim1] > temp_res_log_info[npt-1][nphim1] || sign_of_f11 * sign_of_f21 < 0 ) )
				f1 = 0.0;
			else if (sign_of_f11 * sign_of_f21 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
				f1 = sign_of_f11 * exp( lin_int(ptr-pT0, one_by_pTdiff, temp_res_log_info[npt-1][nphim1], temp_res_log_info[npt][nphim1]) );
			else					// otherwise, just interpolate original vals
				f1 = lin_int(ptr-pT0, one_by_pTdiff, temp_res_moments_info[npt-1][nphim1], temp_res_moments_info[npt][nphim1]);

			//*******************************************************************************************************************
			// set f2 next
			//*******************************************************************************************************************
			if (ptr > pT1 && ( temp_res_log_info[npt][nphi] > temp_res_log_info[npt-1][nphi] || sign_of_f12 * sign_of_f22 < 0 ) )
				f2 = 0.0;
			else if (sign_of_f12 * sign_of_f22 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
				f2 = sign_of_f12 * exp( lin_int(ptr-pT0, one_by_pTdiff, temp_res_log_info[npt-1][nphi], temp_res_log_info[npt][nphi]) );
			else					// otherwise, just interpolate original vals
				f2 = lin_int(ptr-pT0, one_by_pTdiff, temp_res_moments_info[npt-1][nphi], temp_res_moments_info[npt][nphi]);
			//*******************************************************************************************************************
		}
		else						// if pT is smaller than PTCHANGE, just use linear interpolation, no matter what
		{
			f1 = lin_int(ptr-pT0, one_by_pTdiff, temp_res_moments_info[npt-1][nphim1], temp_res_moments_info[npt][nphim1]);
			f2 = lin_int(ptr-pT0, one_by_pTdiff, temp_res_moments_info[npt-1][nphi], temp_res_moments_info[npt][nphi]);

		}
	
		// now, interpolate f1 and f2 over the pphi direction
		results[wfi] = lin_int(phir-phi0, one_by_pphidiff, f1, f2);
	
		if ( isnan( results[wfi] ) )
		{
			*global_out_stream_ptr << "ERROR in Edndp3(double, double, double*): problems encountered!" << endl
				<< "results[" << wfi << "] = " << setw(8) << setprecision(15) << results[wfi] << endl
				<< "  --> ptr = " << ptr << endl
				<< "  --> pt0 = " << pT0 << endl
				<< "  --> pt1 = " << pT1 << endl
				<< "  --> phir = " << phir << endl
				<< "  --> phi0 = " << phi0 << endl
				<< "  --> phi1 = " << phi1 << endl
				<< "  --> f11 = " << temp_res_moments_info[npt-1][nphim1] << endl
				<< "  --> f12 = " << temp_res_moments_info[npt-1][nphi] << endl
				<< "  --> f21 = " << temp_res_moments_info[npt][nphim1] << endl
				<< "  --> f22 = " << temp_res_moments_info[npt][nphi] << endl
				<< "  --> f1 = " << f1 << endl
				<< "  --> f2 = " << f2 << endl;
			exit(1);
		}
		/*if (current_level_of_output > 0) cout << "Edndp3(): results[" << wfi << "] = "
			<< setw(8) << setprecision(15) << results[wfi] << endl
			<< "  --> ptr = " << ptr << endl
			<< "  --> pt0 = " << pT0 << endl
			<< "  --> pt1 = " << pT1 << endl
			<< "  --> phir = " << phir << endl
			<< "  --> phi0 = " << phi0 << endl
			<< "  --> phi1 = " << phi1 << endl
			<< "  --> f11 = " << temp_res_moments_info[npt-1][nphim1] << endl
			<< "  --> f12 = " << temp_res_moments_info[npt-1][nphi] << endl
			<< "  --> f21 = " << temp_res_moments_info[npt][nphim1] << endl
			<< "  --> f22 = " << temp_res_moments_info[npt][nphi] << endl
			<< "  --> f1 = " << f1 << endl
			<< "  --> f2 = " << f2 << endl;*/
	}

	return;
}

//End of file
