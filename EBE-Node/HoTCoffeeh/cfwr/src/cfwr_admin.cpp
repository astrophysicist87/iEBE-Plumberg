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

#include "cfwr.h"
#include "cfwr_lib.h"
#include "Arsenal.h"
#include "chebyshev.h"
#include "gauss_quadrature.h"

using namespace std;

template <typename T> int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

CorrelationFunction::CorrelationFunction(particle_info* particle, particle_info* all_particles_in, int Nparticle_in,
					FO_surf* FOsurf_ptr, vector<int> chosen_resonances_in, int particle_idx, ofstream& myout,
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
	particle_name = particle->name;
	particle_mass = particle->mass;
	particle_sign = particle->sign;
	particle_gspin = particle->gspin;
	particle_id = particle_idx;
	target_particle_id = particle_id;
	particle_monval = particle->monval;
	S_prefactor = 1.0/(8.0*(M_PI*M_PI*M_PI))/hbarC/hbarC/hbarC;
	current_total_resonance_percentage = 0.0;
	all_particles = all_particles_in;
	for (int icr = 0; icr < (int)chosen_resonances_in.size(); icr++)
		chosen_resonances.push_back(chosen_resonances_in[icr]);
	thermal_pions_only = false;
	Nparticle = Nparticle_in;
	NchosenParticle = (int)chosen_resonances_in.size();
	read_in_all_dN_dypTdpTdphi = false;
	output_all_dN_dypTdpTdphi = true;
	currentfolderindex = -1;
	current_level_of_output = 0;
	//qspace_cs_slice_length = qnpts*qnpts*qnpts*qnpts*2;		//factor of 2 for sin or cos
	qspace_cs_slice_length = qtnpts*qxnpts*qynpts*qznpts*2;		//factor of 2 for sin or cos

	gsl_set_error_handler_off();

	number_of_percentage_markers = UDPMsize;

	//set arrays containing q points
	Set_q_points();

	//sort by proximity to origin (where CF is largest) and do those points first
	Set_sorted_q_pts_list();

	n_zeta_pts = zeta_npts;
	n_v_pts = v_npts;
	n_s_pts = s_npts;

	v_min = -1.;
	v_max = 1.;
	zeta_min = 0.;
	zeta_max = M_PI;
	
	//default: use delta_f in calculations
	use_delta_f = true;
	no_df_stem = "";
	current_FOsurf_ptr = FOsurf_ptr;

//****************************************************************************************************
//OLD CODE FOR READING IN SELECTED decay_channels...
	current_resonance_mass = 0.0;
	current_resonance_mu = 0.0;
	current_resonance_Gamma = 0.0;
	current_resonance_total_br = 0.0;
	current_resonance_decay_masses = new double [2];
	current_resonance_decay_masses[0] = 0.0;
	current_resonance_decay_masses[1] = 0.0;
	previous_resonance_particle_id = -1;
	previous_decay_channel_idx = -1;				//different for each decay channel
	previous_resonance_mass = 0.0;
	previous_resonance_Gamma = 0.0;
	previous_resonance_total_br = 0.0;
	if (chosen_resonances.size() == 0)
	{
		n_decay_channels = 1;
		n_resonance = 0;
		thermal_pions_only = true;
		if (VERBOSE > 0) *global_out_stream_ptr << "Thermal pion(+) only!" << endl;
		decay_channels = new decay_info [n_decay_channels];
		decay_channels[0].resonance_decay_masses = new double [Maxdecaypart];	// Maxdecaypart == 5
	}
	else
	{
		//n_decay_channels is actually total number of decay channels which can generate pions
		//from chosen decay_channels
		n_decay_channels = get_number_of_decay_channels(chosen_resonances, all_particles);
		n_resonance = (int)chosen_resonances.size();
		if (VERBOSE > 0) *global_out_stream_ptr << "Computed n_decay_channels = " << n_decay_channels << endl
							<< "Computed n_resonance = " << n_resonance << endl;
		decay_channels = new decay_info [n_decay_channels];
		int temp_idx = 0;
		for (int icr = 0; icr < n_resonance; icr++)
		{
			particle_info particle_temp = all_particles[chosen_resonances[icr]];
			if (VERBOSE > 0) *global_out_stream_ptr << "Loading resonance: " << particle_temp.name
					<< ", chosen_resonances[" << icr << "] = " << chosen_resonances[icr] << endl;
			for (int idecay = 0; idecay < particle_temp.decays; idecay++)
			{
				if (VERBOSE > 0) *global_out_stream_ptr << "Current temp_idx = " << temp_idx << endl;
				if (temp_idx == n_decay_channels)	//i.e., all contributing decay channels have been loaded
					break;
				decay_channels[temp_idx].resonance_name = particle_temp.name;		// set name of resonance

				//check if effective branching is too small for inclusion in source variances
				bool effective_br_is_too_small = false;
				if (particle_temp.decays_effective_branchratio[idecay] <= 1.e-12)
					effective_br_is_too_small = true;

				decay_channels[temp_idx].resonance_particle_id = chosen_resonances[icr];	// set index of resonance in all_particles
				decay_channels[temp_idx].resonance_idx = icr;					// set index of resonance in chosen_resonances
				decay_channels[temp_idx].resonance_decay_masses = new double [Maxdecaypart];	// Maxdecaypart == 5
				decay_channels[temp_idx].resonance_decay_monvals = new int [Maxdecaypart];	// Maxdecaypart == 5
				decay_channels[temp_idx].resonance_decay_Gammas = new double [Maxdecaypart];	// Maxdecaypart == 5



				//*** SETTING RESONANCE DECAY MASSES DIFFERENTLY FOR NEW ANALYZE SF
				for (int ii = 0; ii < Maxdecaypart; ii++)
				{
					decay_channels[temp_idx].resonance_decay_monvals[ii] = particle_temp.decays_part[idecay][ii];
					if (particle_temp.decays_part[idecay][ii] == 0)
					{
						decay_channels[temp_idx].resonance_decay_masses[ii] = 0.0;
						decay_channels[temp_idx].resonance_decay_Gammas[ii] = 0.0;

					}
					else
					{
						int tempID = lookup_particle_id_from_monval(all_particles, Nparticle, particle_temp.decays_part[idecay][ii]);
						decay_channels[temp_idx].resonance_decay_masses[ii] = all_particles[tempID].mass;
						decay_channels[temp_idx].resonance_decay_Gammas[ii] = all_particles[tempID].width;

					}
				}
				decay_channels[temp_idx].resonance_mu = particle_temp.mu;
				decay_channels[temp_idx].resonance_gspin = particle_temp.gspin;
				decay_channels[temp_idx].resonance_sign = particle_temp.sign;
				decay_channels[temp_idx].resonance_mass = particle_temp.mass;
				decay_channels[temp_idx].nbody = abs(particle_temp.decays_Npart[idecay]);
				decay_channels[temp_idx].resonance_Gamma = particle_temp.width;
				decay_channels[temp_idx].resonance_total_br = particle_temp.decays_effective_branchratio[idecay];
				decay_channels[temp_idx].resonance_direct_br = particle_temp.decays_branchratio[idecay];

				
				//check if particle lifetime is too long for inclusion in source variances
				bool lifetime_is_too_long = false;
				if (decay_channels[temp_idx].resonance_Gamma < hbarC / max_lifetime)
					lifetime_is_too_long = true;		//i.e., for lifetimes longer than 100 fm/c, skip decay channel

				if (VERBOSE > 0) *global_out_stream_ptr << "Resonance = " << decay_channels[temp_idx].resonance_name << ", decay channel " << idecay + 1
						<< ": mu=" << decay_channels[temp_idx].resonance_mu
						<< ", gs=" << decay_channels[temp_idx].resonance_gspin << ", sign=" << decay_channels[temp_idx].resonance_sign
						<< ", M=" << decay_channels[temp_idx].resonance_mass << ", nbody=" << decay_channels[temp_idx].nbody
						<< ", Gamma=" << decay_channels[temp_idx].resonance_Gamma << ", total br=" << decay_channels[temp_idx].resonance_total_br
						<< ", direct br=" << decay_channels[temp_idx].resonance_direct_br << endl;

				if (VERBOSE > 0) *global_out_stream_ptr << "Resonance = " << decay_channels[temp_idx].resonance_name << ": ";
				for (int decay_part_idx = 0; decay_part_idx < decay_channels[temp_idx].nbody; decay_part_idx++)
					if (VERBOSE > 0) *global_out_stream_ptr << "m[" << decay_part_idx << "] = "
						<< decay_channels[temp_idx].resonance_decay_masses[decay_part_idx] << "   "
						<< decay_channels[temp_idx].resonance_decay_monvals[decay_part_idx] << "   ";
				if (VERBOSE > 0) *global_out_stream_ptr << endl << endl;

				// if decay channel parent resonance is not too long-lived
				// and decay channel contains at least one target daughter particle,
				// include channel
				decay_channels[temp_idx].include_channel = !lifetime_is_too_long && !effective_br_is_too_small;

				temp_idx++;
			}
		}
	}

	//try flattening
	const int giant_flat_array_size = n_interp_pT_pts * n_interp_pphi_pts * qtnpts * qxnpts * qynpts * qznpts * ntrig;
	thermal_target_dN_dypTdpTdphi_moments = new double [giant_flat_array_size];
	current_dN_dypTdpTdphi_moments = new double [giant_flat_array_size];
	current_ln_dN_dypTdpTdphi_moments = new double [giant_flat_array_size];
	current_sign_of_dN_dypTdpTdphi_moments = new double [giant_flat_array_size];
	for (int i = 0; i < giant_flat_array_size; ++i)
	{
		thermal_target_dN_dypTdpTdphi_moments[i] = 0.0;
		current_dN_dypTdpTdphi_moments[i] = 0.0;
		current_ln_dN_dypTdpTdphi_moments[i] = 0.0;
		current_sign_of_dN_dypTdpTdphi_moments[i] = 0.0;
	}


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

	res_log_info = new double ** [n_interp_pT_pts];
	res_sign_info = new double ** [n_interp_pT_pts];
	res_moments_info = new double ** [n_interp_pT_pts];
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		res_log_info[ipt] = new double * [n_interp_pphi_pts];
		res_sign_info[ipt] = new double * [n_interp_pphi_pts];
		res_moments_info[ipt] = new double * [n_interp_pphi_pts];
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			res_log_info[ipt][ipphi] = new double [qspace_cs_slice_length];
			res_sign_info[ipt][ipphi] = new double [qspace_cs_slice_length];
			res_moments_info[ipt][ipphi] = new double [qspace_cs_slice_length];
		}
	}

	// initialize spectra and correlation function arrays
	spectra = new double ** [Nparticle];
	abs_spectra = new double ** [Nparticle];
	thermal_spectra = new double ** [Nparticle];
	log_spectra = new double ** [Nparticle];
	sign_spectra = new double ** [Nparticle];
	for (int ir = 0; ir < Nparticle; ++ir)
	{
		spectra[ir] = new double * [n_interp_pT_pts];
		abs_spectra[ir] = new double * [n_interp_pT_pts];
		thermal_spectra[ir] = new double * [n_interp_pT_pts];
		log_spectra[ir] = new double * [n_interp_pT_pts];
		sign_spectra[ir] = new double * [n_interp_pT_pts];
		for (int ipT = 0; ipT < n_interp_pT_pts; ++ipT)
		{
			spectra[ir][ipT] = new double [n_interp_pphi_pts];
			abs_spectra[ir][ipT] = new double [n_interp_pphi_pts];
			thermal_spectra[ir][ipT] = new double [n_interp_pphi_pts];
			log_spectra[ir][ipT] = new double [n_interp_pphi_pts];
			sign_spectra[ir][ipT] = new double [n_interp_pphi_pts];
			for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
			{
				spectra[ir][ipT][ipphi] = 0.0;
				abs_spectra[ir][ipT][ipphi] = 0.0;
				thermal_spectra[ir][ipT][ipphi] = 0.0;
				log_spectra[ir][ipT][ipphi] = 0.0;
				sign_spectra[ir][ipT][ipphi] = 0.0;
			}
		}
	}

	flat_spectra = new double [n_interp_pT_pts*n_interp_pphi_pts];
	for (int iii = 0; iii < n_interp_pT_pts*n_interp_pphi_pts; ++iii)
		flat_spectra[iii] = 0.0;

	tmp_moments_real = new double **** [qtnpts];
	tmp_moments_imag = new double **** [qtnpts];
	for (int iqt = 0; iqt < qtnpts; ++iqt)
	{
		tmp_moments_real[iqt] = new double *** [qxnpts];
		tmp_moments_imag[iqt] = new double *** [qxnpts];
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		{
			tmp_moments_real[iqt][iqx] = new double ** [qynpts];
			tmp_moments_imag[iqt][iqx] = new double ** [qynpts];
			for (int iqy = 0; iqy < qynpts; ++iqy)
			{
				tmp_moments_real[iqt][iqx][iqy] = new double * [qznpts];
				tmp_moments_imag[iqt][iqx][iqy] = new double * [qznpts];
				for (int iqz = 0; iqz < qznpts; ++iqz)
				{
					tmp_moments_real[iqt][iqx][iqy][iqz] = new double [n_interp_pT_pts*n_interp_pphi_pts];
					tmp_moments_imag[iqt][iqx][iqy][iqz] = new double [n_interp_pT_pts*n_interp_pphi_pts];
					for (int iii = 0; iii < n_interp_pT_pts*n_interp_pphi_pts; ++iii)
					{
						tmp_moments_real[iqt][iqx][iqy][iqz][iii] = 0.0;
						tmp_moments_imag[iqt][iqx][iqy][iqz][iii] = 0.0;	
					}
				}
			}
		}
	}

	// used for keeping track of how many FO cells are important for given pT, pphi
	// also set up q-space cutoffs array
	number_of_FOcells_above_cutoff_array = new int * [n_interp_pT_pts];
	current_q_space_cutoff = new double * [n_interp_pT_pts];
	correlator_minus_one_cutoff_norms = new int ** [n_interp_pT_pts];
	for (int ipT = 0; ipT < n_interp_pT_pts; ++ipT)
	{
		number_of_FOcells_above_cutoff_array[ipT] = new int [n_interp_pphi_pts];
		current_q_space_cutoff[ipT] = new double [n_interp_pphi_pts];
		correlator_minus_one_cutoff_norms[ipT] = new int * [n_interp_pphi_pts];
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			correlator_minus_one_cutoff_norms[ipT][ipphi] = new int [4];
			for (int ii = 0; ii < 4; ++ii)
				correlator_minus_one_cutoff_norms[ipT][ipphi][ii] = qtnpts*qtnpts + qxnpts*qxnpts + qynpts*qynpts + qznpts*qznpts;
				// i.e., something larger than maximum q-array ranges, only made smaller if correlator cutoff threshhold reached, making large q-values redundant
			number_of_FOcells_above_cutoff_array[ipT][ipphi] = 0;
			current_q_space_cutoff[ipT][ipphi] = 0.0;
		}
	}

	// set-up integration points for resonance integrals
	v_pts = new double [n_v_pts];
	v_wts = new double [n_v_pts];
	zeta_pts = new double [n_zeta_pts];
	zeta_wts = new double [n_zeta_pts];
	//initialize all gaussian points for resonance integrals
	//syntax: int gauss_quadrature(int order, int kind, double alpha, double beta, double a, double b, double x[], double w[])
	gauss_quadrature(n_zeta_pts, 1, 0.0, 0.0, zeta_min, zeta_max, zeta_pts, zeta_wts);
	gauss_quadrature(n_v_pts, 1, 0.0, 0.0, v_min, v_max, v_pts, v_wts);

	//set pT and pphi points
	SPinterp_pT = new double [n_interp_pT_pts];
	SPinterp_pT_wts = new double [n_interp_pT_pts];
	SPinterp_pphi = new double [n_interp_pphi_pts];
	SPinterp_pphi_wts = new double [n_interp_pphi_pts];
	sin_SPinterp_pphi = new double [n_interp_pphi_pts];
	cos_SPinterp_pphi = new double [n_interp_pphi_pts];
	gauss_quadrature(n_interp_pT_pts, 5, 0.0, 0.0, 0.0, 13.0, SPinterp_pT, SPinterp_pT_wts);
	gauss_quadrature(n_interp_pphi_pts, 1, 0.0, 0.0, interp_pphi_min, interp_pphi_max, SPinterp_pphi, SPinterp_pphi_wts);
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

	R2_side_QM = new double * [n_interp_pT_pts];
	R2_out_QM = new double * [n_interp_pT_pts];
	R2_long_QM = new double * [n_interp_pT_pts];
	R2_outside_QM = new double * [n_interp_pT_pts];
	R2_sidelong_QM = new double * [n_interp_pT_pts];
	R2_outlong_QM = new double * [n_interp_pT_pts];

	R2_side_QM_C = new double * [n_localp_T];
	R2_out_QM_C = new double * [n_localp_T];
	R2_long_QM_C = new double * [n_localp_T];
	R2_outside_QM_C = new double * [n_localp_T];
	R2_sidelong_QM_C = new double * [n_localp_T];
	R2_outlong_QM_C = new double * [n_localp_T];

	R2_side_QM_S = new double * [n_localp_T];
	R2_out_QM_S = new double * [n_localp_T];
	R2_long_QM_S = new double * [n_localp_T];
	R2_outside_QM_S = new double * [n_localp_T];
	R2_sidelong_QM_S = new double * [n_localp_T];
	R2_outlong_QM_S = new double * [n_localp_T];

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

		R2_side_QM[ipt] = new double [n_interp_pphi_pts];
		R2_out_QM[ipt] = new double [n_interp_pphi_pts];
		R2_outside_QM[ipt] = new double [n_interp_pphi_pts];
		R2_long_QM[ipt] = new double [n_interp_pphi_pts];
		R2_sidelong_QM[ipt] = new double [n_interp_pphi_pts];
		R2_outlong_QM[ipt] = new double [n_interp_pphi_pts];

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

		R2_side_QM_C[iKT] = new double [n_order];
		R2_out_QM_C[iKT] = new double [n_order];
		R2_outside_QM_C[iKT] = new double [n_order];
		R2_long_QM_C[iKT] = new double [n_order];
		R2_sidelong_QM_C[iKT] = new double [n_order];
		R2_outlong_QM_C[iKT] = new double [n_order];

		R2_side_QM_S[iKT] = new double [n_order];
		R2_out_QM_S[iKT] = new double [n_order];
		R2_outside_QM_S[iKT] = new double [n_order];
		R2_long_QM_S[iKT] = new double [n_order];
		R2_sidelong_QM_S[iKT] = new double [n_order];
		R2_outlong_QM_S[iKT] = new double [n_order];

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

			R2_side_QM_C[iKT][in] = 0.;
			R2_out_QM_C[iKT][in] = 0.;
			R2_long_QM_C[iKT][in] = 0.;
			R2_outside_QM_C[iKT][in] = 0.;
			R2_sidelong_QM_C[iKT][in] = 0.;
			R2_outlong_QM_C[iKT][in] = 0.;

			R2_side_QM_S[iKT][in] = 0.;
			R2_out_QM_S[iKT][in] = 0.;
			R2_long_QM_S[iKT][in] = 0.;
			R2_outside_QM_S[iKT][in] = 0.;
			R2_sidelong_QM_S[iKT][in] = 0.;
			R2_outlong_QM_S[iKT][in] = 0.;
		}
	}

   return;
}

void CorrelationFunction::Allocate_osc_arrays(int FOarray_length)
{
	osc0 = new double *** [qtnpts];
	osc1 = new double ** [qxnpts];
	osc2 = new double ** [qynpts];
	osc3 = new double *** [qznpts];

	//allocate qtpts
	for (int iqt = 0; iqt < qtnpts; ++iqt)
	{
		osc0[iqt] = new double ** [FOarray_length];
		for (int isurf = 0; isurf < FOarray_length; ++isurf)
		{
			osc0[iqt][isurf] = new double * [eta_s_npts];
			for (int ieta = 0; ieta < eta_s_npts; ++ieta)
				osc0[iqt][isurf][ieta] = new double [2];
		}
	}
	//allocate qxpts
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	{
		osc1[iqx] = new double * [FOarray_length];
		for (int isurf = 0; isurf < FOarray_length; ++isurf)
			osc1[iqx][isurf] = new double [2];
	}
	//allocate qypts
	for (int iqy = 0; iqy < qynpts; ++iqy)
	{
		osc2[iqy] = new double * [FOarray_length];
		for (int isurf = 0; isurf < FOarray_length; ++isurf)
			osc2[iqy][isurf] = new double [2];
	}
	//allocate qzpts
	for (int iqz = 0; iqz < qznpts; ++iqz)
	{
		osc3[iqz] = new double ** [FOarray_length];
		for (int isurf = 0; isurf < FOarray_length; ++isurf)
		{
			osc3[iqz][isurf] = new double * [eta_s_npts];
			for (int ieta = 0; ieta < eta_s_npts; ++ieta)
				osc3[iqz][isurf][ieta] = new double [2];
		}
	}

	return;
}

void CorrelationFunction::Delete_osc_arrays()
{
	for (int iqt = 0; iqt < qtnpts; ++iqt)
	{
		for (int isurf = 0; isurf < FO_length; ++isurf)
		{
			for (int ieta = 0; ieta < eta_s_npts; ++ieta)
				delete [] osc0[iqt][isurf][ieta];
			delete [] osc0[iqt][isurf];
		}
		delete [] osc0[iqt];
	}
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	{
		for (int isurf = 0; isurf < FO_length; ++isurf)
			delete [] osc1[iqx][isurf];
		delete [] osc1[iqx];
	}
	for (int iqy = 0; iqy < qynpts; ++iqy)
	{
		for (int isurf = 0; isurf < FO_length; ++isurf)
			delete [] osc2[iqy][isurf];
		delete [] osc2[iqy];
	}
	for (int iqz = 0; iqz < qznpts; ++iqz)
	{
		for (int isurf = 0; isurf < FO_length; ++isurf)
		{
			for (int ieta = 0; ieta < eta_s_npts; ++ieta)
				delete [] osc3[iqz][isurf][ieta];
			delete [] osc3[iqz][isurf];
		}
		delete [] osc3[iqz];
	}
}

void CorrelationFunction::Update_sourcefunction(particle_info* particle, int FOarray_length, int particle_idx)
{
	full_FO_length = FOarray_length * eta_s_npts;

	Allocate_osc_arrays(FOarray_length);

	for (int isurf = 0; isurf < FOarray_length; ++isurf)
	{
		FO_surf * surf = &current_FOsurf_ptr[isurf];

		double tau = surf->tau;
		double xpt = surf->xpt;
		double ypt = surf->ypt;

		for (int iqx = 0; iqx < qxnpts; ++iqx)
		{
			osc1[iqx][isurf][0] = cos(hbarCm1*qx_pts[iqx]*xpt);
			osc1[iqx][isurf][1] = sin(hbarCm1*qx_pts[iqx]*xpt);
		}

		for (int iqy = 0; iqy < qynpts; ++iqy)
		{
			osc2[iqy][isurf][0] = cos(hbarCm1*qy_pts[iqy]*ypt);
			osc2[iqy][isurf][1] = sin(hbarCm1*qy_pts[iqy]*ypt);
		}

		for (int ieta = 0; ieta < eta_s_npts; ++ieta)
		{
			double tpt = tau*ch_eta_s[ieta];
			double zpt = tau*sh_eta_s[ieta];

			for (int iqt = 0; iqt < qtnpts; ++iqt)
			{
				osc0[iqt][isurf][ieta][0] = cos(hbarCm1*qt_pts[iqt]*tpt);
				osc0[iqt][isurf][ieta][1] = sin(hbarCm1*qt_pts[iqt]*tpt);
			}

			for (int iqz = 0; iqz < qznpts; ++iqz)
			{
				osc3[iqz][isurf][ieta][0] = cos(hbarCm1*qz_pts[iqz]*zpt);
				osc3[iqz][isurf][ieta][1] = sin(hbarCm1*qz_pts[iqz]*zpt);
			}
		}
	}

	S_p_withweight_array = new double ** [n_interp_pT_pts];
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		S_p_withweight_array[ipt] = new double * [n_interp_pphi_pts];

   //particle information
   particle_name = particle->name;
   particle_mass = particle->mass;
   particle_sign = particle->sign;
   particle_gspin = particle->gspin;
   particle_id = particle_idx;

	*global_out_stream_ptr << "Inside Update_sourcefunction(...): using fraction_of_resonances = " << fraction_of_resonances << endl;

   FO_length = FOarray_length;

	// set the rest later
	most_important_FOcells = new size_t ** [n_interp_pT_pts];
	for(int ipt=0; ipt<n_interp_pT_pts; ipt++)
		most_important_FOcells[ipt] = new size_t * [n_interp_pphi_pts];

   return;
}

void CorrelationFunction::Fill_out_pts(double * pointsarray, int numpoints, double max_val, int spacing_type)
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

CorrelationFunction::~CorrelationFunction()
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

	Delete_osc_arrays();

   return;
}

void CorrelationFunction::Allocate_CFvals()
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

	return;
}

void CorrelationFunction::Delete_CFvals()
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


void CorrelationFunction::Delete_S_p_withweight_array()
{
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		delete [] S_p_withweight_array[ipt];
	delete [] S_p_withweight_array;

	return;
}

void CorrelationFunction::Allocate_decay_channel_info()
{
	if (VERBOSE > 2) *global_out_stream_ptr << "Reallocating memory for decay channel information..." << endl;
	VEC_n2_v_factor = new double [n_v_pts];
	VEC_n2_zeta_factor = new double * [n_v_pts];
	VEC_n2_P_Y = new double [n_v_pts];
	VEC_n2_MTbar = new double [n_v_pts];
	VEC_n2_DeltaMT = new double [n_v_pts];
	VEC_n2_MTp = new double [n_v_pts];
	VEC_n2_MTm = new double [n_v_pts];
	VEC_n2_MT = new double * [n_v_pts];
	VEC_n2_PPhi_tilde = new double * [n_v_pts];
	VEC_n2_PPhi_tildeFLIP = new double * [n_v_pts];
	VEC_n2_PT = new double * [n_v_pts];
	VEC_n2_Ppm = new double *** [n_v_pts];
	for (int iv = 0; iv < n_v_pts; ++iv)
	{
		VEC_n2_MT[iv] = new double [n_zeta_pts];
		VEC_n2_PPhi_tilde[iv] = new double [n_zeta_pts];
		VEC_n2_PPhi_tildeFLIP[iv] = new double [n_zeta_pts];
		VEC_n2_PT[iv] = new double [n_zeta_pts];
		VEC_n2_zeta_factor[iv] = new double [n_zeta_pts];
		VEC_n2_Ppm[iv] = new double ** [n_zeta_pts];
		for (int izeta = 0; izeta < n_zeta_pts; ++izeta)
		{
			VEC_n2_Ppm[iv][izeta] = new double * [2];		//two corresponds to +/-
			for (int i = 0; i < 2; ++i)
				VEC_n2_Ppm[iv][izeta][i] = new double [4];	//four corresponds to space-time components
		}
	}
	s_pts = new double [n_s_pts];
	s_wts = new double [n_s_pts];
	VEC_pstar = new double [n_s_pts];
	VEC_Estar = new double [n_s_pts];
	VEC_DeltaY = new double [n_s_pts];
	VEC_g_s = new double [n_s_pts];
	VEC_s_factor = new double [n_s_pts];
	VEC_v_factor = new double * [n_s_pts];
	VEC_zeta_factor = new double ** [n_s_pts];
	VEC_Yp = new double [n_s_pts];
	VEC_Ym = new double [n_s_pts];
	VEC_P_Y = new double * [n_s_pts];
	VEC_MTbar = new double * [n_s_pts];
	VEC_DeltaMT = new double * [n_s_pts];
	VEC_MTp = new double * [n_s_pts];
	VEC_MTm = new double * [n_s_pts];
	VEC_MT = new double ** [n_s_pts];
	VEC_PPhi_tilde = new double ** [n_s_pts];
	VEC_PPhi_tildeFLIP = new double ** [n_s_pts];
	VEC_PT = new double ** [n_s_pts];
	VEC_Ppm = new double **** [n_s_pts];
	for(int is = 0; is < n_s_pts; ++is)
	{
		VEC_v_factor[is] = new double [n_v_pts];
		VEC_zeta_factor[is] = new double * [n_v_pts];
		VEC_P_Y[is] = new double [n_v_pts];
		VEC_MTbar[is] = new double [n_v_pts];
		VEC_DeltaMT[is] = new double [n_v_pts];
		VEC_MTp[is] = new double [n_v_pts];
		VEC_MTm[is] = new double [n_v_pts];
		VEC_MT[is] = new double * [n_v_pts];
		VEC_PPhi_tilde[is] = new double * [n_v_pts];
		VEC_PPhi_tildeFLIP[is] = new double * [n_v_pts];
		VEC_PT[is] = new double * [n_v_pts];
		VEC_Ppm[is] = new double *** [n_v_pts];
		for(int iv = 0; iv < n_v_pts; ++iv)
		{
			VEC_MT[is][iv] = new double [n_zeta_pts];
			VEC_PPhi_tilde[is][iv] = new double [n_zeta_pts];
			VEC_PPhi_tildeFLIP[is][iv] = new double [n_zeta_pts];
			VEC_PT[is][iv] = new double [n_zeta_pts];
			VEC_zeta_factor[is][iv] = new double [n_zeta_pts];
			VEC_Ppm[is][iv] = new double ** [n_zeta_pts];
			for (int izeta = 0; izeta < n_zeta_pts; ++izeta)
			{
				VEC_Ppm[is][iv][izeta] = new double * [2];		//two corresponds to +/-
				for (int i = 0; i < 2; ++i)
					VEC_Ppm[is][iv][izeta][i] = new double [4];	//four corresponds to space-time components
			}
		}
	}
	if (VERBOSE > 2) *global_out_stream_ptr << "Reallocated memory for decay channel information." << endl;

	return;
}

void CorrelationFunction::Delete_decay_channel_info()
{
	if (VERBOSE > 2) *global_out_stream_ptr << "Deleting memory for decay channel information..." << endl;
	for(int iv = 0; iv < n_v_pts; ++iv)
	{
		delete [] VEC_n2_MT[iv];
		delete [] VEC_n2_PPhi_tilde[iv];
		delete [] VEC_n2_PPhi_tildeFLIP[iv];
		delete [] VEC_n2_PT[iv];
		delete [] VEC_n2_zeta_factor[iv];
		for (int izeta = 0; izeta < n_zeta_pts; ++izeta)
		{
			for (int i = 0; i < 2; ++i)
				delete [] VEC_n2_Ppm[iv][izeta][i];
			delete [] VEC_n2_Ppm[iv][izeta];
		}
		delete [] VEC_n2_Ppm[iv];
	}
	delete [] VEC_n2_v_factor;
	delete [] VEC_n2_zeta_factor;
	delete [] VEC_n2_P_Y;
	delete [] VEC_n2_MTbar ;
	delete [] VEC_n2_DeltaMT;
	delete [] VEC_n2_MTp;
	delete [] VEC_n2_MTm;
	delete [] VEC_n2_MT;
	delete [] VEC_n2_PPhi_tilde;
	delete [] VEC_n2_PPhi_tildeFLIP;
	delete [] VEC_n2_PT;
	delete [] VEC_n2_Ppm;

	for(int is = 0; is < n_s_pts; ++is)
	{
		for(int iv = 0; iv < n_v_pts; ++iv)
		{
			delete [] VEC_MT[is][iv];
			delete [] VEC_PPhi_tilde[is][iv];
			delete [] VEC_PPhi_tildeFLIP[is][iv];
			delete [] VEC_PT[is][iv];
			delete [] VEC_zeta_factor[is][iv];
			for (int izeta = 0; izeta < n_zeta_pts; ++izeta)
			{
				for (int i = 0; i < 2; ++i)
					delete [] VEC_Ppm[is][iv][izeta][i];
				delete [] VEC_Ppm[is][iv][izeta];
			}
			delete [] VEC_Ppm[is][iv];
		}
		delete [] VEC_v_factor[is];
		delete [] VEC_zeta_factor[is];
		delete [] VEC_P_Y[is];
		delete [] VEC_MTbar[is];
		delete [] VEC_DeltaMT[is];
		delete [] VEC_MTp[is];
		delete [] VEC_MTm[is];
		delete [] VEC_MT[is];
		delete [] VEC_PPhi_tilde[is];
		delete [] VEC_PPhi_tildeFLIP[is];
		delete [] VEC_PT[is];
		delete [] VEC_Ppm[is];
	}
	delete [] s_pts;
	delete [] s_wts;
	delete [] VEC_pstar;
	delete [] VEC_Estar;
	delete [] VEC_DeltaY;
	delete [] VEC_g_s;
	delete [] VEC_s_factor;
	delete [] VEC_v_factor;
	delete [] VEC_zeta_factor;
	delete [] VEC_Yp;
	delete [] VEC_Ym;
	delete [] VEC_P_Y;
	delete [] VEC_MTbar;
	delete [] VEC_DeltaMT;
	delete [] VEC_MTp;
	delete [] VEC_MTm;
	delete [] VEC_MT;
	delete [] VEC_PPhi_tilde;
	delete [] VEC_PPhi_tildeFLIP;
	delete [] VEC_PT;
	delete [] VEC_Ppm;
	if (VERBOSE > 2) *global_out_stream_ptr << "Deleted memory for decay channel information." << endl;

	return;
}

void CorrelationFunction::Set_correlation_function_q_pts()
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
void CorrelationFunction::Set_q_points()
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

	cout << "Output iq*0 = " << iqt0 << "   " << iqx0 << "   " << iqy0 << "   " << iqz0 << endl;

	return;
}

inline int norm (vector<int> v) { int norm2 = 0; for (size_t iv = 0; iv < v.size(); ++iv) norm2+=v[iv]*v[iv]; return (norm2); }

inline bool qpt_comparator (vector<int> i, vector<int> j) { return (norm(i) < norm(j)); }

void CorrelationFunction::Set_sorted_q_pts_list()
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

	cout << "Checking sorting of q-points: " << endl;
	for (size_t iq = 0; iq < sorted_q_pts_list.size(); ++iq)
	{
		sorted_q_pts_list[iq][0] += iqt0;
		sorted_q_pts_list[iq][1] += iqx0;
		sorted_q_pts_list[iq][2] += iqy0;
		sorted_q_pts_list[iq][3] += iqz0;

		cout << "   --> iq = " << iq << ": ";
		for (size_t iqmu = 0; iqmu < 4; ++iqmu)
			cout << sorted_q_pts_list[iq][iqmu] << "   ";
		cout << endl;
	}

	return;
}

// returns points in q-space for computing weighted spectra grid corresponding to to given q and K choices
// weighted spectra grid thus needs to be interpolated at point returned in qgridpts
void CorrelationFunction::Get_q_points(double q1, double q2, double q3, double pT, double pphi, double * qgridpts)
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

bool CorrelationFunction::fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile;
}

//print output to output filestream, one line at a time
void CorrelationFunction::Set_ofstream(ofstream& myout)
{
	global_out_stream_ptr = &myout;

	return;
}

//print output to output filestream, one line at a time
void CorrelationFunction::Set_path(string localpath)
{
	global_path = localpath;

	return;
}

void CorrelationFunction::Set_runfolder(string localrunfolder)
{
	global_runfolder = localrunfolder;

	return;
}

void CorrelationFunction::Set_resultsfolder_stem(string usrdef_resultsfolder_stem)
{
	global_resultsfolder_stem = usrdef_resultsfolder_stem;

	return;
}

void CorrelationFunction::Set_use_delta_f(bool usrdef_usedeltaf)
{
	use_delta_f = usrdef_usedeltaf;
	if (!use_delta_f)
		no_df_stem = "_no_df";
	return;
}

void CorrelationFunction::Set_particle_mass(double usrdef_particle_mass)
{
	particle_mass = usrdef_particle_mass;
	return;
}

void CorrelationFunction::Set_current_FOsurf_ptr(FO_surf* FOsurf_ptr)
{
	current_FOsurf_ptr = FOsurf_ptr;
	return;
}

void CorrelationFunction::Get_current_decay_string(int dc_idx, string * decay_string)
{
	// N.B. - dc_idx == 0 is thermal pion(+)s in calling loop, dc_idx > 0 gives resonance decays
	//      ==> need to use dc_idx - 1 here
	*decay_string = decay_channels[dc_idx - 1].resonance_name + " --->> ";
	int temp_monval, tempID;
	for (int decay_part_idx = 0; decay_part_idx < decay_channels[dc_idx - 1].nbody; decay_part_idx++)
	{
		temp_monval = decay_channels[dc_idx - 1].resonance_decay_monvals[decay_part_idx];
		//if (VERBOSE > 0) *global_out_stream_ptr << "Get_current_decay_string(): temp_monval = " << temp_monval << endl;
		if (temp_monval == 0)
			continue;
		else
		{
			tempID = lookup_particle_id_from_monval(all_particles, Nparticle, temp_monval);
			*decay_string += all_particles[tempID].name;
			if (decay_part_idx < decay_channels[dc_idx - 1].nbody - 1) *decay_string += " + ";
		}
	}
	return;
}

int CorrelationFunction::Set_daughter_list(int parent_pid)
{
	// reset list
	daughter_resonance_indices.clear();
	
//for (int i = 0; i < (int)chosen_resonances.size(); ++i)
//	cout << chosen_resonances[i] << "   " << all_particles[chosen_resonances[i]].name << endl;

	// then re-populate it
	particle_info parent = all_particles[parent_pid];
	if (parent.stable == 1 && parent.decays_Npart[0] == 1)
		return (0);									// no daughters to worry about if parent resonance is actually stable
	int number_of_decays = parent.decays;
	for (int k = 0; k < number_of_decays; k++)		// loop through decays for parent resonance
	{
		int nb = abs(parent.decays_Npart[k]);		// for each decay, nb is the number of daughter particles
		for (int l = 0; l < nb; l++)				// loop through each daughter particle
		{
			int pid = lookup_particle_id_from_monval(all_particles, Nparticle, parent.decays_part[k][l]);
//cout << "Searching for " << pid << "   (" << all_particles[pid].name << ", " << all_particles[pid].effective_branchratio << ")"<< endl;
			if ( all_particles[pid].effective_branchratio >= 1.e-12 || pid == target_particle_id )
				daughter_resonance_indices.insert(pid);		// using a <set> object will automatically remove duplicates and keep pid's in a fixed order
		}
	}

//cout << endl << "Ended up with n_daughter = " << daughter_resonance_indices.size() << " for " << all_particles[parent_pid].name << " (" << parent_pid << ")" << endl;
//int i = 0;
//for (set<int>::iterator it = daughter_resonance_indices.begin(); it != daughter_resonance_indices.end(); ++it)
//	cout << i++ << "   " << *it << "   " << all_particles[*it].name << "   " << all_particles[*it].effective_branchratio << endl;

	// return value is total number of daughters found
	return (daughter_resonance_indices.size());
}

int CorrelationFunction::lookup_resonance_idx_from_particle_id(int pid)
{
	// pid - particle index in all_particles array
	// looks up location in chosen_resonances of given value particle_id
	int result = -1;

	for (int ii = 0; ii < (int)chosen_resonances.size(); ii++)
	{
		if (chosen_resonances[ii] == pid)
		{
			result = ii;
			break;
		}
	}

	// if pid is not one of the chosen_resonances, is not the target daughter (pion(+)), is not stable and has a non-zero effective branching ratio
	if (result < 0 && pid != particle_id && all_particles[pid].stable == 0 && all_particles[pid].effective_branchratio >= 1.e-12)
	{
		*global_out_stream_ptr << " *** lookup_resonance_idx_from_particle_id(): Particle_id = " << pid
					<< " (" << all_particles[pid].name <<") not found in chosen_resonances!" << endl
					<< " *** br = " << all_particles[pid].effective_branchratio << endl;
	}
	return (result);
}

void CorrelationFunction::Allocate_resonance_running_sum_vectors()
{
    ssum_vec = new double [qspace_cs_slice_length];
    vsum_vec = new double [qspace_cs_slice_length];
    zetasum_vec = new double [qspace_cs_slice_length];
    Csum_vec = new double [qspace_cs_slice_length];
}

void CorrelationFunction::Delete_resonance_running_sum_vectors()
{
    delete [] ssum_vec;
    delete [] vsum_vec;
    delete [] zetasum_vec;
    delete [] Csum_vec;
}

void CorrelationFunction::Zero_resonance_running_sum_vector(double * vec)
{
	for (int tmp = 0; tmp < qspace_cs_slice_length; ++tmp)
		vec[tmp] = 0.0;
}

void CorrelationFunction::Setup_current_daughters_dN_dypTdpTdphi_moments(int n_daughter)
{
	const int giant_flat_array_size = n_interp_pT_pts * n_interp_pphi_pts * qtnpts * qxnpts * qynpts * qznpts * ntrig;
	current_daughters_dN_dypTdpTdphi_moments = new double * [n_daughter];
	current_daughters_ln_dN_dypTdpTdphi_moments = new double * [n_daughter];
	current_daughters_sign_of_dN_dypTdpTdphi_moments = new double * [n_daughter];
	for (int id = 0; id < n_daughter; ++id)
	{
		current_daughters_dN_dypTdpTdphi_moments[id] = new double [giant_flat_array_size];
		current_daughters_ln_dN_dypTdpTdphi_moments[id] = new double [giant_flat_array_size];
		current_daughters_sign_of_dN_dypTdpTdphi_moments[id] = new double [giant_flat_array_size];
		for (int i = 0; i < giant_flat_array_size; ++i)
		{
			current_daughters_dN_dypTdpTdphi_moments[id][i] = 0.0;
			current_daughters_ln_dN_dypTdpTdphi_moments[id][i] = 0.0;
			current_daughters_sign_of_dN_dypTdpTdphi_moments[id][i] = 0.0;
		}
	}
	return;
}

void CorrelationFunction::Cleanup_current_daughters_dN_dypTdpTdphi_moments(int n_daughter)
{
	for (int id = 0; id < n_daughter; ++id)
	{
		delete [] current_daughters_dN_dypTdpTdphi_moments[id];
		delete [] current_daughters_ln_dN_dypTdpTdphi_moments[id];
		delete [] current_daughters_sign_of_dN_dypTdpTdphi_moments[id];
	}
	delete [] current_daughters_dN_dypTdpTdphi_moments;
	delete [] current_daughters_ln_dN_dypTdpTdphi_moments;
	delete [] current_daughters_sign_of_dN_dypTdpTdphi_moments;

	return;
}

void CorrelationFunction::Allocate_fleshed_out_CF()
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

void CorrelationFunction::Delete_fleshed_out_CF()
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
