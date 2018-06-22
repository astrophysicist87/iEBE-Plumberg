#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>
#include<time.h>
#include<algorithm>
#include <set>

#include "fitCF.h"
#include "fitCF_lib.h"
#include "Arsenal.h"
#include "chebyshev.h"
#include "gauss_quadrature.h"

using namespace std;

template <typename T> int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

FitCF::FitCF(particle_info* all_particles_in, int Nparticle_in, int particle_idx, vector<int> chosen_events_in, ofstream& myout,
					const int n_interp_pT_pts_in, const int n_interp_pphi_pts_in, const int qtnpts_in, const int qxnpts_in, const int qynpts_in, const int qznpts_in)
{
	//set header info
	n_interp_pT_pts = n_interp_pT_pts_in;
	n_interp_pphi_pts = n_interp_pphi_pts_in;
	qtnpts = qtnpts_in;
	qxnpts = qxnpts_in;
	qynpts = qynpts_in;
	qznpts = qznpts_in;
	init_qt = -0.5*double(qtnpts-1)*delta_qt;
	init_qx = -0.5*double(qxnpts-1)*delta_qx;
	init_qy = -0.5*double(qynpts-1)*delta_qy;
	init_qz = -0.5*double(qznpts-1)*delta_qz;

	//set ofstream for output file
	global_out_stream_ptr = &myout;
	
	//particle information (both final-state particle used in HBT and all decay decay_channels)
	target_particle_id = particle_idx;
	current_total_resonance_percentage = 0.0;

	//set events
	nEvents = (int)chosen_events_in.size();
	for (int iEvent = 0; iEvent < nEvents; ++iEvent)
		chosen_events.push_back( chosen_events_in[iEvent] );

	//actually copy over array of particle info to avoid memory leak issues
	*global_out_stream_ptr << "Filling fitCF.all_particles array..." << endl;

	Nparticle = Nparticle_in;
	all_particles = new particle_info [Nparticle];
	//all_particles = all_particles_in;
	for (int iPart = 0; iPart < Nparticle; ++iPart)
		all_particles[iPart] = all_particles_in[iPart];

	particle = &all_particles[particle_idx];
	chosen_particle = &all_particles[particle_idx];
	particle_name = particle->name;
	particle_mass = particle->mass;
	particle_sign = particle->sign;
	particle_gspin = particle->gspin;
	particle_monval = particle->monval;

	*global_out_stream_ptr << "Finished filling fitCF.all_particles array." << endl;

	thermal_pions_only = false;
	read_in_all_dN_dypTdpTdphi = false;
	output_all_dN_dypTdpTdphi = true;
	currentfolderindex = -1;
	current_level_of_output = 0;
	qspace_cs_slice_length = qtnpts*qxnpts*qynpts*qznpts*2;		//factor of 2 for sin or cos

	gsl_set_error_handler_off();

	//set arrays containing q points
	Set_q_points();

	//sort by proximity to origin (where CF is largest) and do those points first
	Set_sorted_q_pts_list();

	
	//default: use delta_f in calculations
	use_delta_f = true;
	no_df_stem = "";

	thermal_target_dN_dypTdpTdphi_moments = new double [n_interp_pT_pts * n_interp_pphi_pts * qtnpts * qxnpts * qynpts * qznpts * ntrig];
	full_target_dN_dypTdpTdphi_moments = new double [n_interp_pT_pts * n_interp_pphi_pts * qtnpts * qxnpts * qynpts * qznpts * ntrig];

	int qidx = 0;
	qlist = new double * [qtnpts*qxnpts*qynpts*qznpts];
	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	{
		qlist[qidx] = new double [4];
		qlist[qidx][0] = qt_pts[iqt];
		qlist[qidx][1] = qx_pts[iqx];
		qlist[qidx][2] = qy_pts[iqy];
		qlist[qidx][3] = qz_pts[iqz];
		qidx++;
	}

	// initialize spectra and correlation function arrays
	spectra = new double ** [Nparticle];
	avgSpectra = new double ** [Nparticle];
	abs_spectra = new double ** [Nparticle];
	thermal_spectra = new double ** [Nparticle];
	log_spectra = new double ** [Nparticle];
	sign_spectra = new double ** [Nparticle];
	for (int ir = 0; ir < Nparticle; ++ir)
	{
		spectra[ir] = new double * [n_interp_pT_pts];
		avgSpectra[ir] = new double * [n_interp_pT_pts];
		abs_spectra[ir] = new double * [n_interp_pT_pts];
		thermal_spectra[ir] = new double * [n_interp_pT_pts];
		log_spectra[ir] = new double * [n_interp_pT_pts];
		sign_spectra[ir] = new double * [n_interp_pT_pts];
		for (int ipT = 0; ipT < n_interp_pT_pts; ++ipT)
		{
			spectra[ir][ipT] = new double [n_interp_pphi_pts];
			avgSpectra[ir][ipT] = new double [n_interp_pphi_pts];
			abs_spectra[ir][ipT] = new double [n_interp_pphi_pts];
			thermal_spectra[ir][ipT] = new double [n_interp_pphi_pts];
			log_spectra[ir][ipT] = new double [n_interp_pphi_pts];
			sign_spectra[ir][ipT] = new double [n_interp_pphi_pts];
			for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
			{
				spectra[ir][ipT][ipphi] = 0.0;
				avgSpectra[ir][ipT][ipphi] = 0.0;
				abs_spectra[ir][ipT][ipphi] = 0.0;
				thermal_spectra[ir][ipT][ipphi] = 0.0;
				log_spectra[ir][ipT][ipphi] = 0.0;
				sign_spectra[ir][ipT][ipphi] = 0.0;
			}
		}
	}

	//set pT and pphi points
	SPinterp_pT = new double [n_interp_pT_pts];
	SPinterp_pT_wts = new double [n_interp_pT_pts];
	SPinterp_pphi = new double [n_interp_pphi_pts];
	SPinterp_pphi_wts = new double [n_interp_pphi_pts];
	sin_SPinterp_pphi = new double [n_interp_pphi_pts];
	cos_SPinterp_pphi = new double [n_interp_pphi_pts];
	if (USE_OLD_INTERP)
	{
		gauss_quadrature(n_interp_pT_pts, 5, 0.0, 0.0, 0.0, 13.0, SPinterp_pT, SPinterp_pT_wts);
		gauss_quadrature(n_interp_pphi_pts, 1, 0.0, 0.0, interp_pphi_min, interp_pphi_max, SPinterp_pphi, SPinterp_pphi_wts);
	}
	else
	{
		for(int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		{
			double del = 0.5 * (interp_pT_max - interp_pT_min);
			double cen = 0.5 * (interp_pT_max + interp_pT_min);
			SPinterp_pT[ipt] = cen - del * cos( M_PI*(2.*(ipt+1.) - 1.) / (2.*n_interp_pT_pts) );
			//double tmp = - cos( M_PI*(2.*(ipt+1.) - 1.) / (2.*n_interp_pT_pts) );
			//SPinterp_pT[ipt] = ( (1.-sin(M_PI/n_interp_pT_pts))/(1.+sin(M_PI/n_interp_pT_pts)) ) * ( 1. + tmp ) / ( 1. - tmp );
cout << "SPinterp_pT[" << ipt << "] = " << SPinterp_pT[ipt] << endl;
		}
		for(int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			double del = 0.5 * (interp_pphi_max - interp_pphi_min);
			double cen = 0.5 * (interp_pphi_max + interp_pphi_min);
			SPinterp_pphi[ipphi] = cen - del * cos( M_PI*(2.*(ipphi+1.) - 1.) / (2.*n_interp_pphi_pts) );
cout << "SPinterp_pphi[" << ipphi << "] = " << SPinterp_pT[ipphi] << endl;
		}
	}
	for(int ipphi=0; ipphi<n_interp_pphi_pts; ipphi++)
	{
		sin_SPinterp_pphi[ipphi] = sin(SPinterp_pphi[ipphi]);
		cos_SPinterp_pphi[ipphi] = cos(SPinterp_pphi[ipphi]);
	}

	//set p0 and pz points
	SPinterp_p0 = new double * [n_interp_pT_pts];
	SPinterp_pz = new double * [n_interp_pT_pts];
	for(int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		SPinterp_p0[ipt] = new double [eta_s_npts];
		SPinterp_pz[ipt] = new double [eta_s_npts];
	}

	//pair momentum
	K_T = new double [n_localp_T];
	double dK_T = (localp_T_max - localp_T_min)/(n_localp_T - 1 + 1e-100);
	for (int i = 0; i < n_localp_T; ++i)
		K_T[i] = localp_T_min + i*dK_T;
	//K_y = p_y;
	K_y = 0.;
	ch_K_y = cosh(K_y);
	sh_K_y = sinh(K_y);
	beta_l = sh_K_y/ch_K_y;
	K_phi = new double [n_localp_phi];
	K_phi_weight = new double [n_localp_phi];
	gauss_quadrature(n_localp_phi, 1, 0.0, 0.0, localp_phi_min, localp_phi_max, K_phi, K_phi_weight);

	//spatial rapidity grid
	eta_s = new double [eta_s_npts];
	eta_s_weight = new double [eta_s_npts];
	gauss_quadrature(eta_s_npts, 1, 0.0, 0.0, eta_s_i, eta_s_f, eta_s, eta_s_weight);
	ch_eta_s = new double [eta_s_npts];
	sh_eta_s = new double [eta_s_npts];
	for (int ieta = 0; ieta < eta_s_npts; ieta++)
	{
		ch_eta_s[ieta] = cosh(eta_s[ieta]);
		sh_eta_s[ieta] = sinh(eta_s[ieta]);
	}

//debugger(__LINE__, __FILE__);

	//set HBT radii
	R2_side_GF = new double * [n_interp_pT_pts];
	R2_out_GF = new double * [n_interp_pT_pts];
	R2_long_GF = new double * [n_interp_pT_pts];
	R2_outside_GF = new double * [n_interp_pT_pts];
	R2_sidelong_GF = new double * [n_interp_pT_pts];
	R2_outlong_GF = new double * [n_interp_pT_pts];

	R2_side_GF_C = new double * [n_localp_T];
	R2_out_GF_C = new double * [n_localp_T];
	R2_long_GF_C = new double * [n_localp_T];
	R2_outside_GF_C = new double * [n_localp_T];
	R2_sidelong_GF_C = new double * [n_localp_T];
	R2_outlong_GF_C = new double * [n_localp_T];

	R2_side_GF_S = new double * [n_localp_T];
	R2_out_GF_S = new double * [n_localp_T];
	R2_long_GF_S = new double * [n_localp_T];
	R2_outside_GF_S = new double * [n_localp_T];
	R2_sidelong_GF_S = new double * [n_localp_T];
	R2_outlong_GF_S = new double * [n_localp_T];

	R2_side_err = new double * [n_interp_pT_pts];
	R2_out_err = new double * [n_interp_pT_pts];
	R2_long_err = new double * [n_interp_pT_pts];
	R2_outside_err = new double * [n_interp_pT_pts];
	R2_sidelong_err = new double * [n_interp_pT_pts];
	R2_outlong_err = new double * [n_interp_pT_pts];

	lambda_Correl = new double * [n_interp_pT_pts];
	lambda_Correl_err = new double * [n_interp_pT_pts];

	lambda_QM = new double * [n_interp_pT_pts];

	for(int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		R2_side_GF[ipt] = new double [n_interp_pphi_pts];
		R2_out_GF[ipt] = new double [n_interp_pphi_pts];
		R2_outside_GF[ipt] = new double [n_interp_pphi_pts];
		R2_long_GF[ipt] = new double [n_interp_pphi_pts];
		R2_sidelong_GF[ipt] = new double [n_interp_pphi_pts];
		R2_outlong_GF[ipt] = new double [n_interp_pphi_pts];

		R2_side_err[ipt] = new double [n_interp_pphi_pts];
		R2_out_err[ipt] = new double [n_interp_pphi_pts];
		R2_long_err[ipt] = new double [n_interp_pphi_pts];
		R2_outside_err[ipt] = new double [n_interp_pphi_pts];
		R2_sidelong_err[ipt] = new double [n_interp_pphi_pts];
		R2_outlong_err[ipt] = new double [n_interp_pphi_pts];

		lambda_Correl[ipt] = new double [n_interp_pphi_pts];
		lambda_Correl_err[ipt] = new double [n_interp_pphi_pts];

		lambda_QM[ipt] = new double [n_interp_pphi_pts];

		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			R2_side_GF[ipt][ipphi] = 0.;
			R2_out_GF[ipt][ipphi] = 0.;
			R2_long_GF[ipt][ipphi] = 0.;
			R2_outside_GF[ipt][ipphi] = 0.;
			R2_sidelong_GF[ipt][ipphi] = 0.;
			R2_outlong_GF[ipt][ipphi] = 0.;

			R2_side_err[ipt][ipphi] = 0.;
			R2_out_err[ipt][ipphi] = 0.;
			R2_long_err[ipt][ipphi] = 0.;
			R2_outside_err[ipt][ipphi] = 0.;
			R2_sidelong_err[ipt][ipphi] = 0.;
			R2_outlong_err[ipt][ipphi] = 0.;

			lambda_Correl[ipt][ipphi] = 0.0;
			lambda_Correl_err[ipt][ipphi] = 0.0;
		}
	}

	for (int iKT = 0; iKT < n_localp_T; ++iKT)
	{
		R2_side_GF_C[iKT] = new double [n_order];
		R2_out_GF_C[iKT] = new double [n_order];
		R2_outside_GF_C[iKT] = new double [n_order];
		R2_long_GF_C[iKT] = new double [n_order];
		R2_sidelong_GF_C[iKT] = new double [n_order];
		R2_outlong_GF_C[iKT] = new double [n_order];

		R2_side_GF_S[iKT] = new double [n_order];
		R2_out_GF_S[iKT] = new double [n_order];
		R2_outside_GF_S[iKT] = new double [n_order];
		R2_long_GF_S[iKT] = new double [n_order];
		R2_sidelong_GF_S[iKT] = new double [n_order];
		R2_outlong_GF_S[iKT] = new double [n_order];

		for (int in = 0; in < n_order; ++in)
		{
			R2_side_GF_C[iKT][in] = 0.;
			R2_out_GF_C[iKT][in] = 0.;
			R2_long_GF_C[iKT][in] = 0.;
			R2_outside_GF_C[iKT][in] = 0.;
			R2_sidelong_GF_C[iKT][in] = 0.;
			R2_outlong_GF_C[iKT][in] = 0.;

			R2_side_GF_S[iKT][in] = 0.;
			R2_out_GF_S[iKT][in] = 0.;
			R2_long_GF_S[iKT][in] = 0.;
			R2_outside_GF_S[iKT][in] = 0.;
			R2_sidelong_GF_S[iKT][in] = 0.;
			R2_outlong_GF_S[iKT][in] = 0.;
		}
	}

   return;
}


void FitCF::Fill_out_pts(double * pointsarray, int numpoints, double max_val, int spacing_type)
{
	// spacing_type:		0 - uniform spacing
	//						1 - Chebyshev-node based spacing
	if (numpoints == 1)
		pointsarray[0] = 0.0;
	else
	{
		// if I want q-points equally spaced...
		if (spacing_type == 0)
		{
			for (int iqd = 0; iqd < numpoints; ++iqd)
				pointsarray[iqd] = -max_val + (double)iqd * 2.*max_val / double(numpoints - 1+1e-100);
		}
		// else, use Chebyshev nodes instead...
		else if (spacing_type == 1)
		{
			//double local_scale = max_val / cos(M_PI / (2.*qtnpts));
			double local_scale = -max_val;
			for (int iqd = 0; iqd < numpoints; ++iqd)
				pointsarray[iqd] = local_scale * cos( M_PI*(2.*(iqd+1.) - 1.) / (2.*numpoints) );
		}
	}
	return;
}

FitCF::~FitCF()
{
   delete [] K_T;
   delete [] K_phi;
   delete [] K_phi_weight;
   delete [] eta_s;
   delete [] eta_s_weight;

	for(int ipt=0; ipt<n_interp_pT_pts; ipt++)
	{
		delete [] lambda_Correl[ipt];
		delete [] R2_side_GF[ipt];
		delete [] R2_out_GF[ipt];
		delete [] R2_long_GF[ipt];
		delete [] R2_outside_GF[ipt];
		delete [] R2_sidelong_GF[ipt];
		delete [] R2_outlong_GF[ipt];

		delete [] lambda_Correl_err[ipt];
		delete [] R2_side_err[ipt];
		delete [] R2_out_err[ipt];
		delete [] R2_long_err[ipt];
		delete [] R2_outside_err[ipt];
		delete [] R2_sidelong_err[ipt];
		delete [] R2_outlong_err[ipt];
	}

	delete [] R2_side_GF;
	delete [] R2_out_GF;
	delete [] R2_long_GF;
	delete [] R2_outside_GF;
	delete [] R2_sidelong_GF;
	delete [] R2_outlong_GF;

	delete [] lambda_Correl_err;
	delete [] R2_side_err;
	delete [] R2_out_err;
	delete [] R2_long_err;
	delete [] R2_outside_err;
	delete [] R2_sidelong_err;
	delete [] R2_outlong_err;

	for (int ir = 0; ir < Nparticle; ++ir)
	{
		for (int ipT = 0; ipT < n_interp_pT_pts; ++ipT)
			delete [] spectra[ir][ipT];
		delete [] spectra[ir];
	}
	delete [] spectra;

   return;
}

void FitCF::Allocate_CFvals()
{
	CFvals = new double **** [n_interp_pT_pts];
	thermalCFvals = new double **** [n_interp_pT_pts];
	crosstermCFvals = new double **** [n_interp_pT_pts];
	resonancesCFvals = new double **** [n_interp_pT_pts];
	for (int ipT = 0; ipT < n_interp_pT_pts; ++ipT)
	{
		CFvals[ipT] = new double *** [n_interp_pphi_pts];
		thermalCFvals[ipT] = new double *** [n_interp_pphi_pts];
		crosstermCFvals[ipT] = new double *** [n_interp_pphi_pts];
		resonancesCFvals[ipT] = new double *** [n_interp_pphi_pts];
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			CFvals[ipT][ipphi] = new double ** [qxnpts];
			thermalCFvals[ipT][ipphi] = new double ** [qxnpts];
			crosstermCFvals[ipT][ipphi] = new double ** [qxnpts];
			resonancesCFvals[ipT][ipphi] = new double ** [qxnpts];
			for (int iqx = 0; iqx < qxnpts; ++iqx)
			{
				CFvals[ipT][ipphi][iqx] = new double * [qynpts];
				thermalCFvals[ipT][ipphi][iqx] = new double * [qynpts];
				crosstermCFvals[ipT][ipphi][iqx] = new double * [qynpts];
				resonancesCFvals[ipT][ipphi][iqx] = new double * [qynpts];
				for (int iqy = 0; iqy < qynpts; ++iqy)
				{
					CFvals[ipT][ipphi][iqx][iqy] = new double [qznpts];
					thermalCFvals[ipT][ipphi][iqx][iqy] = new double [qznpts];
					crosstermCFvals[ipT][ipphi][iqx][iqy] = new double [qznpts];
					resonancesCFvals[ipT][ipphi][iqx][iqy] = new double [qznpts];
					for (int iqz = 0; iqz < qznpts; ++iqz)
					{
						CFvals[ipT][ipphi][iqx][iqy][iqz] = 0.0;
						thermalCFvals[ipT][ipphi][iqx][iqy][iqz] = 0.0;
						crosstermCFvals[ipT][ipphi][iqx][iqy][iqz] = 0.0;
						resonancesCFvals[ipT][ipphi][iqx][iqy][iqz] = 0.0;
					}
				}
			}
		}
	}
	avgCorrelation_function = new double **** [n_interp_pT_pts];
	avgCorrelation_function_Numerator = new double **** [n_interp_pT_pts];
	avgCorrelation_function_Denominator = new double **** [n_interp_pT_pts];
	avgThermalCFvals = new double **** [n_interp_pT_pts];
	avgCrosstermCFvals = new double **** [n_interp_pT_pts];
	avgResonancesCFvals = new double **** [n_interp_pT_pts];
	for (int ipT = 0; ipT < n_interp_pT_pts; ++ipT)
	{
		avgCorrelation_function[ipT] = new double *** [n_interp_pphi_pts];
		avgCorrelation_function_Numerator[ipT] = new double *** [n_interp_pphi_pts];
		avgCorrelation_function_Denominator[ipT] = new double *** [n_interp_pphi_pts];
		avgThermalCFvals[ipT] = new double *** [n_interp_pphi_pts];
		avgCrosstermCFvals[ipT] = new double *** [n_interp_pphi_pts];
		avgResonancesCFvals[ipT] = new double *** [n_interp_pphi_pts];
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			avgCorrelation_function[ipT][ipphi] = new double ** [qxnpts];
			avgCorrelation_function_Numerator[ipT][ipphi] = new double ** [qxnpts];
			avgCorrelation_function_Denominator[ipT][ipphi] = new double ** [qxnpts];
			avgThermalCFvals[ipT][ipphi] = new double ** [qxnpts];
			avgCrosstermCFvals[ipT][ipphi] = new double ** [qxnpts];
			avgResonancesCFvals[ipT][ipphi] = new double ** [qxnpts];
			for (int iqx = 0; iqx < qxnpts; ++iqx)
			{
				avgCorrelation_function[ipT][ipphi][iqx] = new double * [qynpts];
				avgCorrelation_function_Numerator[ipT][ipphi][iqx] = new double * [qynpts];
				avgCorrelation_function_Denominator[ipT][ipphi][iqx] = new double * [qynpts];
				avgThermalCFvals[ipT][ipphi][iqx] = new double * [qynpts];
				avgCrosstermCFvals[ipT][ipphi][iqx] = new double * [qynpts];
				avgResonancesCFvals[ipT][ipphi][iqx] = new double * [qynpts];
				for (int iqy = 0; iqy < qynpts; ++iqy)
				{
					avgCorrelation_function[ipT][ipphi][iqx][iqy] = new double [qznpts];
					avgCorrelation_function_Numerator[ipT][ipphi][iqx][iqy] = new double [qznpts];
					avgCorrelation_function_Denominator[ipT][ipphi][iqx][iqy] = new double [qznpts];
					avgThermalCFvals[ipT][ipphi][iqx][iqy] = new double [qznpts];
					avgCrosstermCFvals[ipT][ipphi][iqx][iqy] = new double [qznpts];
					avgResonancesCFvals[ipT][ipphi][iqx][iqy] = new double [qznpts];
					for (int iqz = 0; iqz < qznpts; ++iqz)
					{
						avgCorrelation_function[ipT][ipphi][iqx][iqy][iqz] = 0.0;
						avgCorrelation_function_Numerator[ipT][ipphi][iqx][iqy][iqz] = 0.0;
						avgCorrelation_function_Denominator[ipT][ipphi][iqx][iqy][iqz] = 0.0;
						avgThermalCFvals[ipT][ipphi][iqx][iqy][iqz] = 0.0;
						avgCrosstermCFvals[ipT][ipphi][iqx][iqy][iqz] = 0.0;
						avgResonancesCFvals[ipT][ipphi][iqx][iqy][iqz] = 0.0;
					}
				}
			}
		}
	}

	return;
}

void FitCF::Delete_CFvals()
{
	for (int ipT = 0; ipT < n_interp_pT_pts; ++ipT)
	{
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			for (int iqx = 0; iqx < qxnpts; ++iqx)
			{
				for (int iqy = 0; iqy < qynpts; ++iqy)
				{
					delete [] CFvals[ipT][ipphi][iqx][iqy];
					delete [] thermalCFvals[ipT][ipphi][iqx][iqy];
					delete [] crosstermCFvals[ipT][ipphi][iqx][iqy];
					delete [] resonancesCFvals[ipT][ipphi][iqx][iqy];
				}
				delete [] CFvals[ipT][ipphi][iqx];
				delete [] thermalCFvals[ipT][ipphi][iqx];
				delete [] crosstermCFvals[ipT][ipphi][iqx];
				delete [] resonancesCFvals[ipT][ipphi][iqx];
			}
			delete [] CFvals[ipT][ipphi];
			delete [] thermalCFvals[ipT][ipphi];
			delete [] crosstermCFvals[ipT][ipphi];
			delete [] resonancesCFvals[ipT][ipphi];
		}
		delete [] CFvals[ipT];
		delete [] thermalCFvals[ipT];
		delete [] crosstermCFvals[ipT];
		delete [] resonancesCFvals[ipT];
	}
	delete [] CFvals;
	delete [] thermalCFvals;
	delete [] crosstermCFvals;
	delete [] resonancesCFvals;

	return;
}

void FitCF::Set_correlation_function_q_pts()
{
	q1npts = qxnpts;
	q2npts = qynpts;
	q3npts = qznpts;
	
	// initialize error matrix
	Correl_3D_err = new double ** [qxnpts];
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	{
		Correl_3D_err[iqx] = new double * [qynpts];
		for (int iqy = 0; iqy < qynpts; ++iqy)
		{
			Correl_3D_err[iqx][iqy] = new double [qznpts];
			for (int iqz = 0; iqz < qznpts; ++iqz)
				Correl_3D_err[iqx][iqy][iqz] = 1e-3;	//naive choice for now
		}
	}

	return;
}

// sets points in q-space for computing weighted spectra grid
void FitCF::Set_q_points()
{
	q_pts = new double [qnpts];
	for (int iq = 0; iq < qnpts; ++iq)
		q_pts[iq] = init_q + (double)iq * delta_q;
	q_axes = new double [3];

	qt_pts = new double [qtnpts];
	qx_pts = new double [qxnpts];
	qy_pts = new double [qynpts];
	qz_pts = new double [qznpts];

	double local_pT_max = SPinterp_pT[n_interp_pT_pts-1];	//max pT value

	//q0 maximized by maximizing pT, maximizing qo, and setting qs==ql==0
	double mpion = all_particles[target_particle_id].mass;
	double qxmax = abs(init_qx);
	double qymax = abs(init_qy);
	double qxymax = sqrt(qxmax*qxmax+qymax*qymax);
	double xi2 = mpion*mpion + interp_pT_max*interp_pT_max + 0.25*qxymax*qxymax;	//pretend that Kphi == 0, qx == qo and qs == ql == 0, to maximize qtmax
	double qtmax = sqrt(xi2 + interp_pT_max*qxymax) - sqrt(xi2 - interp_pT_max*qxymax) + 1.e-10;

	Fill_out_pts(qt_pts, qtnpts, qtmax, QT_POINTS_SPACING);
	Fill_out_pts(qx_pts, qxnpts, abs(init_qx), QX_POINTS_SPACING);
	Fill_out_pts(qy_pts, qynpts, abs(init_qy), QY_POINTS_SPACING);
	Fill_out_pts(qz_pts, qznpts, abs(init_qz), QZ_POINTS_SPACING);

	iqt0 = (qtnpts - 1) / 2;
	iqx0 = (qxnpts - 1) / 2;
	iqy0 = (qynpts - 1) / 2;
	iqz0 = (qznpts - 1) / 2;

	//cout << "Output iq*0 = " << iqt0 << "   " << iqx0 << "   " << iqy0 << "   " << iqz0 << endl;

	return;
}

inline int norm (vector<int> v) { int norm2 = 0; for (size_t iv = 0; iv < v.size(); ++iv) norm2+=v[iv]*v[iv]; return (norm2); }

inline bool qpt_comparator (vector<int> i, vector<int> j) { return (norm(i) < norm(j)); }

void FitCF::Set_sorted_q_pts_list()
{
	vector<int> tmp (4);

	//sort through all q-points with qt <= 0 by proximity to origin in q-space
	//allows to do largest CF values first
	for (int iqt = 0; iqt < (qtnpts / 2) + 1; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	{
		tmp[0] = iqt - iqt0;
		tmp[1] = iqx - iqx0;
		tmp[2] = iqy - iqy0;
		tmp[3] = iqz - iqz0;

		sorted_q_pts_list.push_back(tmp);
	}

	sort(sorted_q_pts_list.begin(), sorted_q_pts_list.end(), qpt_comparator);

	//cout << "Checking sorting of q-points: " << endl;
	for (size_t iq = 0; iq < sorted_q_pts_list.size(); ++iq)
	{
		sorted_q_pts_list[iq][0] += iqt0;
		sorted_q_pts_list[iq][1] += iqx0;
		sorted_q_pts_list[iq][2] += iqy0;
		sorted_q_pts_list[iq][3] += iqz0;

		//cout << "   --> iq = " << iq << ": ";
		//for (size_t iqmu = 0; iqmu < 4; ++iqmu)
		//	cout << sorted_q_pts_list[iq][iqmu] << "   ";
		//cout << endl;
	}

	return;
}

// returns points in q-space for computing weighted spectra grid corresponding to to given q and K choices
// weighted spectra grid thus needs to be interpolated at point returned in qgridpts
void FitCF::Get_q_points(double q1, double q2, double q3, double pT, double pphi, double * qgridpts)
{
	double mtarget = all_particles[target_particle_id].mass;
	double xi2 = mtarget*mtarget + pT*pT + 0.25*(q1*q1 + q2*q2 + q3*q3);
	double ckp = cos(pphi), skp = sin(pphi);

	double qo = ckp * q1 + skp * q2;
		
	// set qpts at which to interpolate spectra
	qgridpts[0] = sqrt(xi2 + qo*pT) - sqrt(xi2 - qo*pT);	//set qt component
	qgridpts[1] = q1;										//set qx component
	qgridpts[2] = q2;										//set qy component
	qgridpts[3] = q3;										//set qz component, since qz = ql

	return;
}

bool FitCF::fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile;
}

//print output to output filestream, one line at a time
void FitCF::Set_path(string localpath)
{
	global_path = localpath;

	return;
}

void FitCF::Set_runfolder(string localrunfolder)
{
	global_runfolder = localrunfolder;

	return;
}

void FitCF::Set_resultsfolder_stem(string usrdef_resultsfolder_stem)
{
	global_resultsfolder_stem = usrdef_resultsfolder_stem;

	return;
}

void FitCF::Allocate_fleshed_out_CF()
{
	fleshed_out_CF = new double ** [new_nqpts];
	fleshed_out_thermal = new double ** [new_nqpts];
	fleshed_out_crossterm = new double ** [new_nqpts];
	fleshed_out_resonances = new double ** [new_nqpts];
	for (int iqx = 0; iqx < new_nqpts; ++iqx)
	{
		fleshed_out_CF[iqx] = new double * [new_nqpts];
		fleshed_out_thermal[iqx] = new double * [new_nqpts];
		fleshed_out_crossterm[iqx] = new double * [new_nqpts];
		fleshed_out_resonances[iqx] = new double * [new_nqpts];
		for (int iqy = 0; iqy < new_nqpts; ++iqy)
		{
			fleshed_out_CF[iqx][iqy] = new double [new_nqpts];
			fleshed_out_thermal[iqx][iqy] = new double [new_nqpts];
			fleshed_out_crossterm[iqx][iqy] = new double [new_nqpts];
			fleshed_out_resonances[iqx][iqy] = new double [new_nqpts];
			for (int iqz = 0; iqz < new_nqpts; ++iqz)
			{
				fleshed_out_CF[iqx][iqy][iqz] = 0.0;
				fleshed_out_thermal[iqx][iqy][iqz] = 0.0;
				fleshed_out_crossterm[iqx][iqy][iqz] = 0.0;
				fleshed_out_resonances[iqx][iqy][iqz] = 0.0;
			}
		}
	}

	qx_fleshed_out_pts = new double [new_nqpts];
	qy_fleshed_out_pts = new double [new_nqpts];
	qz_fleshed_out_pts = new double [new_nqpts];

	return;
}

void FitCF::Set_use_delta_f(bool usrdef_usedeltaf)
{
	use_delta_f = usrdef_usedeltaf;
	if (!use_delta_f)
		no_df_stem = "_no_df";
	return;
}

void FitCF::Delete_fleshed_out_CF()
{
	for (int iqx = 0; iqx < new_nqpts; ++iqx)
	{
		for (int iqy = 0; iqy < new_nqpts; ++iqy)
		{
			delete [] fleshed_out_CF[iqx][iqy];
			delete [] fleshed_out_thermal[iqx][iqy];
			delete [] fleshed_out_crossterm[iqx][iqy];
			delete [] fleshed_out_resonances[iqx][iqy];
		}
		delete [] fleshed_out_CF[iqx];
		delete [] fleshed_out_thermal[iqx];
		delete [] fleshed_out_crossterm[iqx];
		delete [] fleshed_out_resonances[iqx];
	}
	delete [] fleshed_out_CF;
	delete [] fleshed_out_thermal;
	delete [] fleshed_out_crossterm;
	delete [] fleshed_out_resonances;

	delete [] qx_fleshed_out_pts;
	delete [] qy_fleshed_out_pts;
	delete [] qz_fleshed_out_pts;

	return;
}

//End of file
