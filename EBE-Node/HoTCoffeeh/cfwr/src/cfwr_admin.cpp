#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <time.h>
#include <algorithm>
#include <set>

#include "../include/cfwr.h"
#include "../include/cfwr_lib.h"
#include "../include/Arsenal.h"
#include "../include/chebyshev.h"
#include "../include/Stopwatch.h"
#include "../include/gauss_quadrature.h"
#include "../include/bessel.h"

using namespace std;

template <typename T> int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

struct bessel_params { double beta; double gamma; };

CorrelationFunction::CorrelationFunction(
	ParameterReader * paraRdr_in,
	particle_info* particle,
	particle_info* all_particles_in,
	int Nparticle_in,
	vector<int> chosen_resonances_in,
	int particle_idx, ofstream& myout)
{
	paraRdr = paraRdr_in;

	USE_PLANE_PSI_ORDER
		= paraRdr->getVal("use_plane_psi_order");
	INCLUDE_DELTA_F
		= paraRdr->getVal("include_delta_f");
	GROUPING_PARTICLES
		= paraRdr->getVal("grouping_particles");
	PARTICLE_DIFF_TOLERANCE
		= paraRdr->getVal("particle_diff_tolerance");
	IGNORE_LONG_LIVED_RESONANCES
		= paraRdr->getVal("ignore_long_lived_resonances");
	FIT_WITH_PROJECTED_CFVALS
		= paraRdr->getVal("fit_with_projected_cfvals");
	FLESH_OUT_CF
		= paraRdr->getVal("flesh_out_cf");
	CALCULATE_CF_MODE
		= paraRdr->getVal("calculate_CF_mode");

	n_order 		= paraRdr->getVal("n_order");
	tol 			= paraRdr->getVal("tolerance");
	flagneg 		= paraRdr->getVal("flag_negative_S");
	max_lifetime	= paraRdr->getVal("max_lifetime");

	//Set header info
	//Define various grid sizes
	// - SP momentum points at which to evaluate correlation function
	n_pT_pts 		= paraRdr->getVal("CF_npT");
	n_pphi_pts 		= paraRdr->getVal("CF_npphi");
	n_pY_pts 		= paraRdr->getVal("CF_npY");
	// - pair momenta points at which to interpolate HBT results
	nKT 			= paraRdr->getVal("nKT");
	nKphi 			= paraRdr->getVal("nKphi");
	KT_min 			= paraRdr->getVal("KTmin");
	KT_max 			= paraRdr->getVal("KTmax");
	// - relative momentum points at which to evaluate
	//   correlation function
	qtnpts 			= paraRdr->getVal("qtnpts");
	qxnpts 			= paraRdr->getVal("qxnpts");
	qynpts 			= paraRdr->getVal("qynpts");
	qznpts 			= paraRdr->getVal("qznpts");
	// - step size in q directions
	delta_qt 		= 0.1;
	delta_qx 		= paraRdr->getVal("delta_qx");
	delta_qy 		= paraRdr->getVal("delta_qy");
	delta_qz 		= paraRdr->getVal("delta_qz");
	// - minimum value in each q direction
	init_qt 		= -0.5*double(qtnpts-1)*delta_qt;
	init_qx 		= -0.5*double(qxnpts-1)*delta_qx;
	init_qy 		= -0.5*double(qynpts-1)*delta_qy;
	init_qz 		= -0.5*double(qznpts-1)*delta_qz;

	// - used for rapid Chebyshev evaluation of complex Bessel functions;
	//   31 points works reasonably well
	n_alpha_points 	= 31;

	// - number of points to use when fleshing out correlation
	//   function in each direction
	new_nqxpts 		= ( qxnpts > 1 ) ? new_nqpts : 1;
	new_nqypts 		= ( qynpts > 1 ) ? new_nqpts : 1;
	new_nqzpts 		= ( qznpts > 1 ) ? new_nqpts : 1;

	//set ofstream for output file
	out 			= &myout;
	
	//particle information (both for the final-state particle
	//used in HBT correlation function and for all decay channels)
	particle_name 	= particle->name;
	particle_mass 	= particle->mass;
	particle_sign 	= particle->sign;
	particle_gspin 	= particle->gspin;
	particle_id 	= particle_idx;
	target_particle_id
					= particle_id;
	particle_monval = particle->monval;

	current_total_resonance_percentage
					= 0.0;

	// load all particle information
	all_particles 	= all_particles_in;

	// load list of chosen resonances to use in calculation
	for (int icr = 0; icr < (int)chosen_resonances_in.size(); icr++)
	{
		chosen_resonances.push_back(chosen_resonances_in[icr]);
		//cout << "chosen_resonances_in[" << icr << "] = "
		//		<< chosen_resonances_in[icr] << endl;
	}

	// total number of particles included in database
	Nparticle 		= Nparticle_in;

	// total number of resonances being used
	NchosenParticle = (int)chosen_resonances_in.size();

	// set some defaults
	thermal_pions_only 			= false;
	read_in_all_dN_dypTdpTdphi 	= false;
	output_all_dN_dypTdpTdphi 	= true;
	current_level_of_output 	= 0;

	//Length of q-vectors, counting complexity and even/odd:
	//factor of 4 for sin or cos, real or imaginary
	qspace_cs_slice_length 		= qxnpts*qynpts*ntrig;

	//Define index of the point where pY==0
	ipY0 = ( USE_RAPIDITY_SYMMETRY ) ? 0 : (n_pY_pts - 1) / 2;

	//gsl_set_error_handler_off();

	//set arrays containing q points
	Set_q_points();

	//Define limits of phase-space integrals
	v_min 		= -1.;
	v_max 		= 1.;
	zeta_min 	= 0.;
	zeta_max 	= M_PI;
	
	//Set whether viscous corrections included
	Set_use_delta_f();

	// Set up keeping track of current and previous resonances
	current_resonance_mass 				= 0.0;
	current_resonance_mu 				= 0.0;
	current_resonance_Gamma 			= 0.0;
	current_resonance_total_br 			= 0.0;
	current_resonance_decay_masses 		= new double [2];
	current_resonance_decay_masses[0] 	= 0.0;
	current_resonance_decay_masses[1] 	= 0.0;
	previous_resonance_particle_id 		= -1;
	// this is different for each decay channel,
	// not just each resonance
	previous_decay_channel_idx 			= -1;
	previous_resonance_mass 			= 0.0;
	previous_resonance_Gamma 			= 0.0;
	previous_resonance_total_br 		= 0.0;

	// thermal pions only
	if (chosen_resonances.size() == 0)
	{
		n_decay_channels = 1;
		n_resonance = 0;
		thermal_pions_only = true;
		if (VERBOSE > 0)
			*out << "Thermal pion(+) only!" << endl;
		decay_channels = new decay_info [n_decay_channels];
		decay_channels[0].resonance_decay_masses
			= new double [Maxdecaypart];
	}
	// if there are resonances to consider
	else
	{
		// n_decay_channels is total number of decay channels which
		// contribute to pion yield from chosen decay_channels
		n_decay_channels
			= get_number_of_decay_channels(chosen_resonances, all_particles);
		n_resonance
			= (int)chosen_resonances.size();

		if (VERBOSE > 0)
			*out << "Computed n_decay_channels = "
					<< n_decay_channels << endl
					<< "Computed n_resonance = "
					<< n_resonance << endl;

		decay_channels
			= new decay_info [n_decay_channels];

		// now loop over all resonances and load decay channel info
		int temp_idx = 0;
		for (int icr = 0; icr < n_resonance; icr++)
		{
			particle_info particle_temp
				= all_particles[chosen_resonances[icr]];

			if (VERBOSE > 0)
				*out << "Loading resonance: "
						<< particle_temp.name
						<< ", chosen_resonances[" << icr << "] = "
						<< chosen_resonances[icr] << endl;

			for (int idecay = 0; idecay < particle_temp.decays; idecay++)
			{
				if (VERBOSE > 0)
					*out << "Current temp_idx = " << temp_idx << endl;

				// if all contributing decay channels have been loaded
				if (temp_idx == n_decay_channels)
					break;

				decay_info decay_to_load;

				// set name of resonance
				decay_to_load.resonance_name = particle_temp.name;

				// check if effective branching is too small for inclusion
				// in source variances
				bool effective_br_is_too_small = false;
				if (particle_temp.decays_effective_branchratio[idecay]
						<= 1.e-12)
					effective_br_is_too_small = true;

				// set index of resonance in all_particles
				decay_to_load.resonance_particle_id
						= chosen_resonances[icr];
				// set index of resonance in chosen_resonances
				decay_to_load.resonance_idx
						= icr;
				decay_to_load.resonance_decay_monvals
						= new int [Maxdecaypart];
				decay_to_load.resonance_decay_masses
						= new double [Maxdecaypart];
				decay_to_load.resonance_decay_Gammas
						= new double [Maxdecaypart];

				// Set resonance decay masses
				for (int ii = 0; ii < Maxdecaypart; ii++)
				{
					decay_to_load.resonance_decay_monvals[ii]
						= particle_temp.decays_part[idecay][ii];

					if (particle_temp.decays_part[idecay][ii] == 0)
					{
						decay_to_load.resonance_decay_masses[ii]
							= 0.0;
						decay_to_load.resonance_decay_Gammas[ii]
							= 0.0;
					}
					else
					{
						int tempID
							= lookup_particle_id_from_monval(
								all_particles,
								Nparticle,
								particle_temp.decays_part[idecay][ii]);
						decay_to_load.resonance_decay_masses[ii]
							= all_particles[tempID].mass;
						decay_to_load.resonance_decay_Gammas[ii]
							= all_particles[tempID].width;
					}
				}
				decay_to_load.resonance_mu
							= particle_temp.mu;
				decay_to_load.resonance_gspin
							= particle_temp.gspin;
				decay_to_load.resonance_sign
							= particle_temp.sign;
				decay_to_load.resonance_mass
							= particle_temp.mass;
				decay_to_load.nbody
							= abs(particle_temp.decays_Npart[idecay]);
				decay_to_load.resonance_Gamma
							= particle_temp.width;
				decay_to_load.resonance_total_br
							= particle_temp.decays_effective_branchratio[idecay];
				decay_to_load.resonance_direct_br
							= particle_temp.decays_branchratio[idecay];
				
				//check if particle lifetime is too long for
				//inclusion in source variances
				bool lifetime_is_too_long = false;
				//if (decay_channels[temp_idx].resonance_Gamma
				//		< hbarC / max_lifetime)
				//	lifetime_is_too_long = true;
				//i.e., for lifetimes longer than 100 fm/c, skip decay channel

				if (VERBOSE > 0)
					*out << "Resonance = "
							<< decay_to_load.resonance_name
							<< ", decay channel " << idecay + 1
							<< ": mu=" << decay_to_load.resonance_mu
							<< ", gs=" << decay_to_load.resonance_gspin
							<< ", sign=" << decay_to_load.resonance_sign
							<< ", M=" << decay_to_load.resonance_mass
							<< ", nbody=" << decay_to_load.nbody
							<< ", Gamma=" << decay_to_load.resonance_Gamma
							<< ", total br=" << decay_to_load.resonance_total_br
							<< ", direct br="
							<< decay_to_load.resonance_direct_br << endl;

				if (VERBOSE > 0)
					*out << "Resonance = "
							<< decay_to_load.resonance_name << ": ";

				for (int dec_pt_idx = 0;
						dec_pt_idx < decay_to_load.nbody;
						dec_pt_idx++)
				if (VERBOSE > 0)
					*out << "m[" << dec_pt_idx << "] = "
							<< decay_to_load.resonance_decay_masses[dec_pt_idx]
							<< "   "
							<< decay_to_load.resonance_decay_monvals[dec_pt_idx]
							<< "   ";

				if (VERBOSE > 0) *out << endl << endl;

				// if decay channel parent resonance is not too long-lived
				// and decay channel contains at least one target daughter
				// particle, include channel
				decay_to_load.include_channel
					= !lifetime_is_too_long && !effective_br_is_too_small;

				decay_channels[temp_idx++]
					= decay_to_load;
			}
		}
	}

	/////////////////////////
	//Create flattened arrays

	// - space-time moments in a given (qt,qz)-loop
	const int giant_flat_array_size
		= n_pT_pts * n_pphi_pts * n_pY_pts
			* qxnpts * qynpts * ntrig;
	thermal_target_dN_dypTdpTdphi_moments
		= new double [giant_flat_array_size];
	current_dN_dypTdpTdphi_moments
		= new double [giant_flat_array_size];
	for (int i = 0; i < giant_flat_array_size; ++i)
	{
		thermal_target_dN_dypTdpTdphi_moments[i] 	= 0.0;
		current_dN_dypTdpTdphi_moments[i] 			= 0.0;
	}

	// - all space-time moments
	const int full_size
		= n_pT_pts * n_pphi_pts
			* qtnpts * qxnpts * qynpts * qznpts * ntrig;
	thermal_target_Yeq0_moments
		= new double [full_size];
	full_target_Yeq0_moments
		= new double [full_size];
	for (int i = 0; i < full_size; ++i)
	{
		thermal_target_Yeq0_moments[i] 	= 0.0;
		full_target_Yeq0_moments[i] 	= 0.0;
	}

	// Set the list of q points
	int qidx = 0;
	qlist = new double * [qxnpts*qynpts];
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	{
		qlist[qidx] = new double [4];
		++qidx;
	}

	// - flattened arrays for resonance interpolation
	res_log_info
		= new double * [n_pT_pts * n_pphi_pts * n_pY_pts];
	res_sign_info
		= new double * [n_pT_pts * n_pphi_pts * n_pY_pts];
	res_moments_info
		= new double * [n_pT_pts * n_pphi_pts * n_pY_pts];
	for (int i = 0; i < n_pT_pts * n_pphi_pts * n_pY_pts; ++i)
	{
		res_log_info[i]
			= new double [qspace_cs_slice_length];
		res_sign_info[i]
			= new double [qspace_cs_slice_length];
		res_moments_info[i]
			= new double [qspace_cs_slice_length];
		for (int ii = 0; ii < qspace_cs_slice_length; ++ii)
		{
			res_log_info[i][ii] 	= 0.0;
			res_sign_info[i][ii] 	= 0.0;
			res_moments_info[i][ii] = 0.0;
		}
	}

	// initialize spectra and correlation function arrays
	spectra 		= new double ** [Nparticle];
	thermal_spectra = new double ** [Nparticle];
	log_spectra 	= new double ** [Nparticle];
	sign_spectra 	= new double ** [Nparticle];
	for (int ir = 0; ir < Nparticle; ++ir)
	{
		spectra[ir] 		= new double * [n_pT_pts];
		thermal_spectra[ir] = new double * [n_pT_pts];
		log_spectra[ir] 	= new double * [n_pT_pts];
		sign_spectra[ir] 	= new double * [n_pT_pts];
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		{
			spectra[ir][ipT] 			= new double [n_pphi_pts];
			thermal_spectra[ir][ipT] 	= new double [n_pphi_pts];
			log_spectra[ir][ipT] 		= new double [n_pphi_pts];
			sign_spectra[ir][ipT] 		= new double [n_pphi_pts];
			for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
			{
				spectra[ir][ipT][ipphi] 		= 0.0;
				thermal_spectra[ir][ipT][ipphi] = 0.0;
				log_spectra[ir][ipT][ipphi] 	= 0.0;
				sign_spectra[ir][ipT][ipphi] 	= 0.0;
			}
		}
	}

	// set-up integration points for resonance integrals
	v_pts 		= new double [n_v_pts];
	v_wts 		= new double [n_v_pts];
	zeta_pts 	= new double [n_zeta_pts];
	zeta_wts 	= new double [n_zeta_pts];

	//initialize all gaussian points for resonance integrals
	//syntax: int gauss_quadrature(
	//				int order, int kind,
	//				double alpha, double beta,
	//				double a, double b,
	//				double x[], double w[])
	// - points for zeta integral
	gauss_quadrature(n_zeta_pts, 1, 0.0, 0.0,
						zeta_min, zeta_max, zeta_pts, zeta_wts);
	// - points for v integral
	gauss_quadrature(n_v_pts, 1, 0.0, 0.0,
						v_min, v_max, v_pts, v_wts);

	/////////////////////////////
	//set pT, pphi, and pY points

	// - pT points
	SP_pT 		= new double [n_pT_pts];
	SP_pT_wts 	= new double [n_pT_pts];
	gauss_quadrature(n_pT_pts, 5, 0.0, 0.0, 0.0,
						13.0*n_pT_pts/15.0, SP_pT, SP_pT_wts);
	//stratify_npts( 1e-6, 1.0+1e-6, 4.0, int(0.75*n_pT_pts), n_pT_pts, SP_pT );

	//for (int ipt = 0; ipt < n_pT_pts; ++ipt)
	//	cout << ipt << "   " << SP_pT[ipt] << endl;

	// - pphi points
	SP_pphi 	= new double [n_pphi_pts];
	SP_pphi_wts = new double [n_pphi_pts];
	sin_SP_pphi = new double [n_pphi_pts];
	cos_SP_pphi = new double [n_pphi_pts];

	gauss_quadrature(n_pphi_pts, 1, 0.0, 0.0,
						Kphi_min, Kphi_max, SP_pphi, SP_pphi_wts);

	for(int ipphi=0; ipphi < n_pphi_pts; ipphi++)
	{
		sin_SP_pphi[ipphi] = sin(SP_pphi[ipphi]);
		cos_SP_pphi[ipphi] = cos(SP_pphi[ipphi]);
		//cout << ipphi << "   " << SP_pphi[ipphi] << endl;
	}

	// - pY points (Chebyshev nodes useful for resonance interpolation)
	SP_Del_pY 		= new double [n_pY_pts];
	// pY points are fixed by two options:
	//	- USE_RAPIDITY_SYMMETRY determines whether to assume special symmetry
	//	  of space-time moments under certain reflections
	//	- USE_ADJUSTED_MINIMUM places the smallest pY node at Del_pY == 0
	double tmp_tan 	= tan(M_PI / (4.0*n_pY_pts));
	adjusted_SP_Del_pY_minimum
					= ( USE_RAPIDITY_SYMMETRY and USE_ADJUSTED_MINIMUM ) ?
						- SP_Del_pY_max*tmp_tan*tmp_tan
						: SP_Del_pY_min;

	// needed to define interpolation nodes
	double local_scale 	= 0.5 * (adjusted_SP_Del_pY_minimum - SP_Del_pY_max);
	double local_center = 0.5 * (adjusted_SP_Del_pY_minimum + SP_Del_pY_max);
	
	// set interpolation nodes
	double znodes[n_pY_pts];
	for (int ipY = 0; ipY < n_pY_pts; ++ipY)
	{
		znodes[ipY] 	= - cos( M_PI * ( 2.*(ipY+1.) - 1. )
								/ (2.*n_pY_pts) );
		SP_Del_pY[ipY] 	= local_center
							+ local_scale
							* cos( M_PI * ( 2.*(ipY+1.) - 1. )
								/ (2.*n_pY_pts) );
		//cout << ipY << "   " << SP_Del_pY[ipY] << endl;
	}

	// initialize here
	ch_SP_pY = new double [n_pY_pts];
	sh_SP_pY = new double [n_pY_pts];
	for (int ipY = 0; ipY < n_pY_pts; ++ipY)
	{
		ch_SP_pY[ipY] = 0.0;
		sh_SP_pY[ipY] = 0.0;
	}

	// perform necessary Chebyshev function evaluations here
	double chebTnorm[n_pY_pts];
	chebTcfs = new double [n_pY_pts * n_pY_pts];
	for (int ii = 0; ii < n_pY_pts; ++ii)
	{
		chebTnorm[ii] = 0.0;
		for (int kk = 0; kk < n_pY_pts; ++kk)
		{
			double tf = csf::Tfun(ii,znodes[kk]);
			chebTcfs[ii * n_pY_pts + kk] = tf;
			chebTnorm[ii] += tf*tf;
		}
		for (int kk = 0; kk < n_pY_pts; ++kk)
		{
			chebTcfs[ii * n_pY_pts + kk] /= chebTnorm[ii];
			if (ii==0) chebTcfs[ii * n_pY_pts + kk] *= 2.0;
		}
	}
	// Usage: contract appropriate segment of chebTcfs vector with weighted resonance spectra
	//	at fixed pT, pphi to get Chebyshev coeffcients to use in interpolation
	//  CAUTION: chebTcfs are NOT (yet) the coefficients which enter the final Chebyshev sum!!!


	// set spatial rapidity grid
	eta_s = new double [eta_s_npts];
	eta_s_weight = new double [eta_s_npts];
	gauss_quadrature(eta_s_npts, 1, 0.0, 0.0,
						eta_s_i, eta_s_f, eta_s, eta_s_weight);
	ch_eta_s = new double [eta_s_npts];
	sh_eta_s = new double [eta_s_npts];
	for (int ieta = 0; ieta < eta_s_npts; ieta++)
	{
		ch_eta_s[ieta] = cosh(eta_s[ieta]);
		sh_eta_s[ieta] = sinh(eta_s[ieta]);
	}

	//alternative spatial rapidity grid (POSSIBLY DELETE)
	base_Del_eta_s = new double [base_Del_eta_s_npts];
	base_Del_eta_s_weight = new double [base_Del_eta_s_npts];
	gauss_quadrature(base_Del_eta_s_npts, 6, 0.0, 0.0, 0.0, 1.0,
						base_Del_eta_s, base_Del_eta_s_weight);

	// vector for flow plane angles (0 by default)
	plane_angle = new double [n_order];

	//pair momentum
	// - KT
	K_T 			= new double [nKT];
	double dK_T 	= (KT_max - KT_min)/(nKT - 1 + 1e-100);
	for (int i = 0; i < nKT; ++i) K_T[i] = KT_min + i*dK_T;

	// - Kphi
	K_phi 			= new double [nKphi];
	K_phi_weight 	= new double [nKphi];
	gauss_quadrature(nKphi, 1, 0.0, 0.0,
						Kphi_min, Kphi_max, K_phi, K_phi_weight);

	// - KY
	K_y 			= 0.;
	ch_K_y 			= cosh(K_y);
	sh_K_y 			= sinh(K_y);
	beta_l 			= sh_K_y/ch_K_y;

	// set HBT radii
	R2_side_GF 			= new double * [n_pT_pts];
	R2_out_GF 			= new double * [n_pT_pts];
	R2_long_GF 			= new double * [n_pT_pts];
	R2_outside_GF 		= new double * [n_pT_pts];
	R2_sidelong_GF 		= new double * [n_pT_pts];
	R2_outlong_GF 		= new double * [n_pT_pts];

	R2_side_GF_C 		= new double * [nKT];
	R2_out_GF_C 		= new double * [nKT];
	R2_long_GF_C 		= new double * [nKT];
	R2_outside_GF_C 	= new double * [nKT];
	R2_sidelong_GF_C 	= new double * [nKT];
	R2_outlong_GF_C 	= new double * [nKT];

	R2_side_GF_S 		= new double * [nKT];
	R2_out_GF_S 		= new double * [nKT];
	R2_long_GF_S 		= new double * [nKT];
	R2_outside_GF_S 	= new double * [nKT];
	R2_sidelong_GF_S 	= new double * [nKT];
	R2_outlong_GF_S 	= new double * [nKT];

	R2_side_err 		= new double * [n_pT_pts];
	R2_out_err 			= new double * [n_pT_pts];
	R2_long_err 		= new double * [n_pT_pts];
	R2_outside_err 		= new double * [n_pT_pts];
	R2_sidelong_err 	= new double * [n_pT_pts];
	R2_outlong_err 		= new double * [n_pT_pts];

	R2_side_QM 			= new double * [n_pT_pts];
	R2_out_QM 			= new double * [n_pT_pts];
	R2_long_QM 			= new double * [n_pT_pts];
	R2_outside_QM 		= new double * [n_pT_pts];
	R2_sidelong_QM 		= new double * [n_pT_pts];
	R2_outlong_QM 		= new double * [n_pT_pts];

	R2_side_QM_C 		= new double * [nKT];
	R2_out_QM_C 		= new double * [nKT];
	R2_long_QM_C 		= new double * [nKT];
	R2_outside_QM_C 	= new double * [nKT];
	R2_sidelong_QM_C 	= new double * [nKT];
	R2_outlong_QM_C 	= new double * [nKT];

	R2_side_QM_S 		= new double * [nKT];
	R2_out_QM_S 		= new double * [nKT];
	R2_long_QM_S 		= new double * [nKT];
	R2_outside_QM_S 	= new double * [nKT];
	R2_sidelong_QM_S 	= new double * [nKT];
	R2_outlong_QM_S 	= new double * [nKT];

	lambda_Correl 		= new double * [n_pT_pts];
	lambda_Correl_err 	= new double * [n_pT_pts];

	lambda_QM 			= new double * [n_pT_pts];

	for(int ipt = 0; ipt < n_pT_pts; ++ipt)
	{
		R2_side_GF[ipt] 		= new double [n_pphi_pts];
		R2_out_GF[ipt] 			= new double [n_pphi_pts];
		R2_outside_GF[ipt] 		= new double [n_pphi_pts];
		R2_long_GF[ipt] 		= new double [n_pphi_pts];
		R2_sidelong_GF[ipt] 	= new double [n_pphi_pts];
		R2_outlong_GF[ipt] 		= new double [n_pphi_pts];

		R2_side_QM[ipt] 		= new double [n_pphi_pts];
		R2_out_QM[ipt] 			= new double [n_pphi_pts];
		R2_outside_QM[ipt] 		= new double [n_pphi_pts];
		R2_long_QM[ipt] 		= new double [n_pphi_pts];
		R2_sidelong_QM[ipt] 	= new double [n_pphi_pts];
		R2_outlong_QM[ipt] 		= new double [n_pphi_pts];

		R2_side_err[ipt] 		= new double [n_pphi_pts];
		R2_out_err[ipt] 		= new double [n_pphi_pts];
		R2_long_err[ipt] 		= new double [n_pphi_pts];
		R2_outside_err[ipt] 	= new double [n_pphi_pts];
		R2_sidelong_err[ipt] 	= new double [n_pphi_pts];
		R2_outlong_err[ipt] 	= new double [n_pphi_pts];

		lambda_Correl[ipt] 		= new double [n_pphi_pts];
		lambda_Correl_err[ipt] 	= new double [n_pphi_pts];

		lambda_QM[ipt] 			= new double [n_pphi_pts];

		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			R2_side_GF[ipt][ipphi] 			= 0.0;
			R2_out_GF[ipt][ipphi] 			= 0.0;
			R2_long_GF[ipt][ipphi] 			= 0.0;
			R2_outside_GF[ipt][ipphi] 		= 0.0;
			R2_sidelong_GF[ipt][ipphi] 		= 0.0;
			R2_outlong_GF[ipt][ipphi] 		= 0.0;

			R2_side_err[ipt][ipphi] 		= 0.0;
			R2_out_err[ipt][ipphi] 			= 0.0;
			R2_long_err[ipt][ipphi] 		= 0.0;
			R2_outside_err[ipt][ipphi] 		= 0.0;
			R2_sidelong_err[ipt][ipphi] 	= 0.0;
			R2_outlong_err[ipt][ipphi] 		= 0.0;

			lambda_Correl[ipt][ipphi] 		= 0.0;
			lambda_Correl_err[ipt][ipphi] 	= 0.0;
		}
	}

	for (int iKT = 0; iKT < nKT; ++iKT)
	{
		R2_side_GF_C[iKT] 		= new double [n_order];
		R2_out_GF_C[iKT] 		= new double [n_order];
		R2_outside_GF_C[iKT] 	= new double [n_order];
		R2_long_GF_C[iKT] 		= new double [n_order];
		R2_sidelong_GF_C[iKT] 	= new double [n_order];
		R2_outlong_GF_C[iKT] 	= new double [n_order];

		R2_side_GF_S[iKT] 		= new double [n_order];
		R2_out_GF_S[iKT] 		= new double [n_order];
		R2_outside_GF_S[iKT] 	= new double [n_order];
		R2_long_GF_S[iKT] 		= new double [n_order];
		R2_sidelong_GF_S[iKT] 	= new double [n_order];
		R2_outlong_GF_S[iKT] 	= new double [n_order];

		R2_side_QM_C[iKT] 		= new double [n_order];
		R2_out_QM_C[iKT] 		= new double [n_order];
		R2_outside_QM_C[iKT] 	= new double [n_order];
		R2_long_QM_C[iKT] 		= new double [n_order];
		R2_sidelong_QM_C[iKT] 	= new double [n_order];
		R2_outlong_QM_C[iKT] 	= new double [n_order];

		R2_side_QM_S[iKT] 		= new double [n_order];
		R2_out_QM_S[iKT] 		= new double [n_order];
		R2_outside_QM_S[iKT] 	= new double [n_order];
		R2_long_QM_S[iKT] 		= new double [n_order];
		R2_sidelong_QM_S[iKT] 	= new double [n_order];
		R2_outlong_QM_S[iKT] 	= new double [n_order];

		for (int in = 0; in < n_order; ++in)
		{
			R2_side_GF_C[iKT][in] 		= 0.;
			R2_out_GF_C[iKT][in] 		= 0.;
			R2_long_GF_C[iKT][in] 		= 0.;
			R2_outside_GF_C[iKT][in] 	= 0.;
			R2_sidelong_GF_C[iKT][in] 	= 0.;
			R2_outlong_GF_C[iKT][in] 	= 0.;

			R2_side_GF_S[iKT][in] 		= 0.;
			R2_out_GF_S[iKT][in] 		= 0.;
			R2_long_GF_S[iKT][in] 		= 0.;
			R2_outside_GF_S[iKT][in] 	= 0.;
			R2_sidelong_GF_S[iKT][in] 	= 0.;
			R2_outlong_GF_S[iKT][in] 	= 0.;

			R2_side_QM_C[iKT][in] 		= 0.;
			R2_out_QM_C[iKT][in] 		= 0.;
			R2_long_QM_C[iKT][in] 		= 0.;
			R2_outside_QM_C[iKT][in] 	= 0.;
			R2_sidelong_QM_C[iKT][in] 	= 0.;
			R2_outlong_QM_C[iKT][in] 	= 0.;

			R2_side_QM_S[iKT][in] 		= 0.;
			R2_out_QM_S[iKT][in] 		= 0.;
			R2_long_QM_S[iKT][in] 		= 0.;
			R2_outside_QM_S[iKT][in] 	= 0.;
			R2_sidelong_QM_S[iKT][in] 	= 0.;
			R2_outlong_QM_S[iKT][in] 	= 0.;
		}
	}

	return;
}


void CorrelationFunction::Set_qlist(int iqt, int iqz)
{
	int qidx = 0;
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	{
		qlist[qidx][0] = qt_pts[iqt];
		qlist[qidx][1] = qx_pts[iqx];
		qlist[qidx][2] = qy_pts[iqy];
		qlist[qidx][3] = qz_pts[iqz];
		qidx++;
	}
}

void CorrelationFunction::Set_eiqx_matrices()
{
	oscx = new double * [FO_length];
	oscy = new double * [FO_length];

	for (int isurf = 0; isurf < FO_length; ++isurf)
	{
		FO_surf * surf 	= &FOsurf_ptr[isurf];
		oscx[isurf] 	= new double [qxnpts * 2];
		oscy[isurf] 	= new double [qynpts * 2];

		double tau 		= surf->tau;
		double xpt 		= surf->xpt;
		double ypt 		= surf->ypt;

		for (int iqx = 0; iqx < qxnpts; ++iqx)
		{
			oscx[isurf][iqx * 2 + 0] = cos(hbarCm1*qx_pts[iqx]*xpt);
			oscx[isurf][iqx * 2 + 1] = sin(hbarCm1*qx_pts[iqx]*xpt);
		}

		for (int iqy = 0; iqy < qynpts; ++iqy)
		{
			oscy[isurf][iqy * 2 + 0] = cos(hbarCm1*qy_pts[iqy]*ypt);
			oscy[isurf][iqy * 2 + 1] = sin(hbarCm1*qy_pts[iqy]*ypt);
		}
	}

	return;
}

CorrelationFunction::~CorrelationFunction()
{
   delete [] K_T;
   delete [] K_phi;
   delete [] K_phi_weight;

	for(int ipt=0; ipt<n_pT_pts; ipt++)
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
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
			delete [] spectra[ir][ipT];
		delete [] spectra[ir];
	}
	delete [] spectra;

   return;
}

void CorrelationFunction::Allocate_CFvals()
{
	CFvals 				= new double **** [n_pT_pts];
	thermalCFvals 		= new double **** [n_pT_pts];
	crosstermCFvals 	= new double **** [n_pT_pts];
	resonancesCFvals 	= new double **** [n_pT_pts];
	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	{
		CFvals[ipT] 			= new double *** [n_pphi_pts];
		thermalCFvals[ipT] 		= new double *** [n_pphi_pts];
		crosstermCFvals[ipT] 	= new double *** [n_pphi_pts];
		resonancesCFvals[ipT] 	= new double *** [n_pphi_pts];
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			CFvals[ipT][ipphi] 				= new double ** [qxnpts];
			thermalCFvals[ipT][ipphi] 		= new double ** [qxnpts];
			crosstermCFvals[ipT][ipphi] 	= new double ** [qxnpts];
			resonancesCFvals[ipT][ipphi] 	= new double ** [qxnpts];
			for (int iqx = 0; iqx < qxnpts; ++iqx)
			{
				CFvals[ipT][ipphi][iqx] 			= new double * [qynpts];
				thermalCFvals[ipT][ipphi][iqx] 		= new double * [qynpts];
				crosstermCFvals[ipT][ipphi][iqx] 	= new double * [qynpts];
				resonancesCFvals[ipT][ipphi][iqx] 	= new double * [qynpts];
				for (int iqy = 0; iqy < qynpts; ++iqy)
				{
					CFvals[ipT][ipphi][iqx][iqy]
						= new double [qznpts];
					thermalCFvals[ipT][ipphi][iqx][iqy]
						= new double [qznpts];
					crosstermCFvals[ipT][ipphi][iqx][iqy]
						= new double [qznpts];
					resonancesCFvals[ipT][ipphi][iqx][iqy]
						= new double [qznpts];
					for (int iqz = 0; iqz < qznpts; ++iqz)
					{
						CFvals[ipT][ipphi][iqx][iqy][iqz] 			= 0.0;
						thermalCFvals[ipT][ipphi][iqx][iqy][iqz] 	= 0.0;
						crosstermCFvals[ipT][ipphi][iqx][iqy][iqz] 	= 0.0;
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
	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	{
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
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

void CorrelationFunction::Allocate_decay_channel_info()
{
	if (VERBOSE > 2)
		*out << "Reallocating memory for decay channel information..."
				<< endl;

	const int n_spacetime_dims = 4;

	VEC_n2_v_factor 		= new double [n_v_pts];
	VEC_n2_P_Y 				= new double [n_v_pts];
	VEC_n2_zeta_factor 		= new double [n_v_pts * n_zeta_pts];
	VEC_n2_PPhi_tilde 		= new double [n_v_pts * n_zeta_pts];
	VEC_n2_PPhi_tildeFLIP 	= new double [n_v_pts * n_zeta_pts];
	VEC_n2_PT 				= new double [n_v_pts * n_zeta_pts];
	VEC_n2_Ppm 				= new double * [n_v_pts * n_zeta_pts * 2];

	for(int ii = 0; ii < n_v_pts * n_zeta_pts * 2; ++ii)
		VEC_n2_Ppm[ii] 		= new double [n_spacetime_dims];

	s_pts 					= new double [n_s_pts];
	s_wts 					= new double [n_s_pts];
	VEC_n3_g_s 				= new double [n_s_pts];
	VEC_n3_s_factor 		= new double [n_s_pts];
	VEC_n3_v_factor 		= new double [n_s_pts * n_v_pts];
	VEC_n3_zeta_factor 		= new double [n_s_pts * n_v_pts * n_zeta_pts];
	VEC_n3_P_Y 				= new double [n_s_pts * n_v_pts];
	VEC_n3_PPhi_tilde 		= new double [n_s_pts * n_v_pts * n_zeta_pts];
	VEC_n3_PPhi_tildeFLIP 	= new double [n_s_pts * n_v_pts * n_zeta_pts];
	VEC_n3_PT 				= new double [n_s_pts * n_v_pts * n_zeta_pts];
	VEC_n3_Ppm 				= new double * [n_s_pts * n_v_pts * n_zeta_pts * 2];

	for(int ii = 0; ii < n_s_pts * n_v_pts * n_zeta_pts * 2; ++ii)
		VEC_n3_Ppm[ii] 		= new double [n_spacetime_dims];

	if (VERBOSE > 2)
		*out << "Reallocated memory for decay channel information."
				<< endl;

	return;
}

void CorrelationFunction::Delete_decay_channel_info()
{
	if (VERBOSE > 2)
		*out << "Deleting memory for decay channel information..."
				<< endl;

	for(int ii = 0; ii < n_v_pts * n_zeta_pts * 2; ++ii)
		delete [] VEC_n2_Ppm[ii];
	for(int ii = 0; ii < n_s_pts * n_v_pts * n_zeta_pts * 2; ++ii)
		delete [] VEC_n3_Ppm[ii];

	delete [] VEC_n2_v_factor;
	delete [] VEC_n2_zeta_factor;
	delete [] VEC_n2_P_Y;
	delete [] VEC_n2_PPhi_tilde;
	delete [] VEC_n2_PPhi_tildeFLIP;
	delete [] VEC_n2_PT;
	delete [] VEC_n2_Ppm;

	delete [] s_pts;
	delete [] s_wts;
	delete [] VEC_n3_g_s;
	delete [] VEC_n3_s_factor;
	delete [] VEC_n3_v_factor;
	delete [] VEC_n3_zeta_factor;
	delete [] VEC_n3_P_Y;
	delete [] VEC_n3_PPhi_tilde;
	delete [] VEC_n3_PPhi_tildeFLIP;
	delete [] VEC_n3_PT;
	delete [] VEC_n3_Ppm;

	if (VERBOSE > 2)
		*out << "Deleted memory for decay channel information."
				<< endl;

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
	qt_pts 				= new double [qtnpts];
	qx_pts 				= new double [qxnpts];
	qy_pts 				= new double [qynpts];
	qz_pts 				= new double [qznpts];

	double local_pT_max = SP_pT[n_pT_pts-1];	//max pT value

	//q0 maximized by maximizing pT, maximizing qo, and setting qs==ql==0
	double mtarget 		= all_particles[target_particle_id].mass;
	double qxmax 		= abs(init_qx);
	double qymax 		= abs(init_qy);
	qx_max 				= abs(init_qx);
	qy_max 				= abs(init_qy);
	qz_max 				= abs(init_qz);

	double qxymax 		= sqrt(qxmax*qxmax+qymax*qymax);
	//pretend that Kphi == 0, qx == qo and qs == ql == 0, to maximize qtmax
	double xi2 			= mtarget*mtarget
							+ SP_pT_max*SP_pT_max
							+ 0.25*qxymax*qxymax;
	//THIS IS THE INCORRECT VERSION
	// - ONLY NEEDED FOR COMPARISONS WITH PREVIOUS CODE VERSIONS!!!
	//qtmax = sqrt(xi2 + SP_pT_max*qxymax)
	//			- sqrt(xi2 - SP_pT_max*qxymax)
	//			+ 1.e-10;

	//just choose the biggest value (THIS IS THE CORRECT VERSION!!!!!)
	qtmax = max( qxymax, qz_max );

	Fill_out_pts(qt_pts, qtnpts, qtmax, QT_POINTS_SPACING);
	Fill_out_pts(qx_pts, qxnpts, abs(init_qx), QX_POINTS_SPACING);
	Fill_out_pts(qy_pts, qynpts, abs(init_qy), QY_POINTS_SPACING);
	Fill_out_pts(qz_pts, qznpts, abs(init_qz), QZ_POINTS_SPACING);

	iqt0 = (qtnpts - 1) / 2;
	iqx0 = (qxnpts - 1) / 2;
	iqy0 = (qynpts - 1) / 2;
	iqz0 = (qznpts - 1) / 2;

	return;
}

void CorrelationFunction::Fill_out_pts(
		double * pointsarray, int numpoints,
		double max_val, int spacing_type)
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
				pointsarray[iqd]
					= - max_val + (double)iqd * 2.*max_val
									/ double(numpoints - 1+1e-100);
		}
		// else, use Chebyshev nodes instead...
		else if (spacing_type == 1)
		{
			//double local_scale = max_val / cos(M_PI / (2.*qtnpts));
			double local_scale = -max_val;
			for (int iqd = 0; iqd < numpoints; ++iqd)
				pointsarray[iqd]
					= local_scale * cos( M_PI * ( 2.*(iqd+1.) - 1. )
											/ (2.*numpoints) );
		}
		else if (spacing_type == 2)
		{
			//break range of points into two sections (+ve and -ve)
			const int n 	= (numpoints + 1) / 2;

			double tmptan 	= tan(M_PI / (4.0*n));
			//choose definitions to guarantee that middle point is 0
			double a1 		= -max_val, b1 = max_val*tmptan*tmptan;	
			double a2 		= -max_val*tmptan*tmptan, b2 = max_val;
			double hw1 		= 0.5 * (b1 - a1), hw2 = 0.5 * (b2 - a2);
			
			double z1[n], z2[n-1];
			for (int iqd = 0; iqd < n; ++iqd)
			{
				z1[iqd] = -cos( M_PI*(2.*(iqd+1.) - 1.) / (2.*n) );
				pointsarray[iqd] = (1.0 + z1[iqd]) * hw1 + a1;
			}
			for (int iqd = 0; iqd < n; ++iqd)
			{
				z2[iqd] = -cos( M_PI*(2.*(iqd+1.) - 1.) / (2.*n) );
				pointsarray[n + iqd - 1] = (1.0 + z2[iqd]) * hw2 + a2;
			}
		}
		else if (spacing_type == 3)
		{
			//use Chebyshev on infinite interval
			double tmptan = 25.0*tan(M_PI/(2.0*numpoints));

			for (int iqd = 0; iqd < numpoints; ++iqd)
			{
				double ii = iqd + 1.0;
				pointsarray[iqd]
					= - max_val * tmptan
						* cot(M_PI*(2.0*ii-1.0)/(2.0*numpoints));
			}
		}
		else if (spacing_type == 4)
		{
			//break range of points into two sections (+ve and -ve)
			const int n 	= (numpoints + 1) / 2;

			double tmpsin 	= sin(3.0*M_PI / (4.0*n));
			//double tmpcos1 = cos(3.0*M_PI / (4.0*n));
			double tmpcos2 	= cos(M_PI / (2.0*n));
			double tmpcos3 	= cos(3.0*M_PI / (2.0*n));
			
			for (int iqd = 0; iqd < n; ++iqd)
			{
				double ii 			= iqd + 1.0;
				double numerator 	= 2.0 * max_val * tmpsin * tmpsin
										* ( cos(M_PI*(2.0*ii-1.0)/(2.0*n))
											+ tmpcos2 );
				double denominator 	= ( cos(M_PI*(2.0*ii-1.0)/(2.0*n)) - 1.0 )
										* (tmpcos2 + tmpcos3);
				pointsarray[iqd] 	= numerator / denominator;
				//pointsarray[iqd] = -pointsarray[numpoints-iqd-1];
				pointsarray[numpoints-iqd-1]
									= -pointsarray[iqd];
			}
		}
	}
	return;
}



// returns points in q-space for computing weighted spectra grid
// corresponding to to given q and K choices
// weighted spectra grid thus needs to be interpolated at the point
// which isreturned in qgridpts
void CorrelationFunction::Get_q_points(
		double q1, double q2, double q3,
		double pT, double pphi, double * qgridpts)
{
	double mtarget 	= all_particles[target_particle_id].mass;
	double xi2 		= mtarget*mtarget + pT*pT + 0.25*(q1*q1 + q2*q2 + q3*q3);
	double ckp 		= cos(pphi);
	double skp 		= sin(pphi);

	double qo 		= ckp * q1 + skp * q2;
		
	// set qpts at which to interpolate spectra
	qgridpts[0] 	= sqrt(xi2 + qo*pT) - sqrt(xi2 - qo*pT);	//qt
	qgridpts[1] 	= q1;										//qx
	qgridpts[2] 	= q2;										//qy
	qgridpts[3] 	= q3;										//qz

	return;
}

bool CorrelationFunction::fexists(const char *filename)
{
	ifstream ifile(filename);
	return ifile;
}

//print output to output filestream, one line at a time
void CorrelationFunction::Set_path(string path_in)
{
	path = path_in;
	return;
}

void CorrelationFunction::Set_use_delta_f()
{
	use_delta_f = INCLUDE_DELTA_F;
	if (!INCLUDE_DELTA_F)
		no_df_stem = "_no_df";
	return;
}

void CorrelationFunction::Set_FOsurf_ptr(FO_surf* FOsurf_ptr_in, int FO_length_in)
{
	FOsurf_ptr 	= FOsurf_ptr_in;
	FO_length 	= FO_length_in;

	Set_eiqx_matrices();

	// 3D vector to keep track of which FOcells to include
	FOcells_to_include.resize(FO_length*n_pT_pts, vector<int>(n_pphi_pts));

	//initialize to -1; ensures that all cells get set correctly
	Reset_FOcells_array();

	return;
}

void CorrelationFunction::Get_current_decay_string(
		int dc_idx, string * decay_string)
{
	// N.B. - dc_idx == 0 is thermal pion(+)s in calling loop,
	// dc_idx > 0 gives resonance decays
	//      ==> need to use dc_idx - 1 here
	decay_info this_decay = decay_channels[dc_idx - 1];

	// initialize decay string
	*decay_string = this_decay.resonance_name + " --->> ";

	int temp_monval, tempID;
	for (int decay_part_idx = 0;
			decay_part_idx < this_decay.nbody;
			decay_part_idx++)
	{
		temp_monval
			= this_decay.resonance_decay_monvals[decay_part_idx];
		if (VERBOSE > 3)
			*out << "Get_current_decay_string(): temp_monval = "
					<< temp_monval << endl;
		if (temp_monval == 0)
			continue;
		else
		{
			tempID 			= lookup_particle_id_from_monval(
								all_particles, Nparticle, temp_monval);
			*decay_string  += all_particles[tempID].name;
			if (decay_part_idx < this_decay.nbody - 1)
				*decay_string += " + ";
		}
	}
	return;
}

int CorrelationFunction::Set_daughter_list(int parent_pid)
{
	// reset list
	daughter_resonance_indices.clear();
	
	//for (int i = 0; i < (int)chosen_resonances.size(); ++i)
	//	cout << chosen_resonances[i] << "   "
	//			<< all_particles[chosen_resonances[i]].name << endl;

	// then re-populate the list
	particle_info parent = all_particles[parent_pid];

	// no daughters to worry about
	// if parent resonance is actually stable
	if (parent.stable == 1 && parent.decays_Npart[0] == 1)
		return (0);

	// otherwise, see how many decays to consider
	int number_of_decays = parent.decays;

	// loop through decays for parent resonance
	for (int k = 0; k < number_of_decays; k++)
	{
		// for each decay, nb is the number of daughter particles
		int nb = abs(parent.decays_Npart[k]);

		// loop through each daughter particle
		for (int l = 0; l < nb; l++)
		{
			int pid = lookup_particle_id_from_monval(
						all_particles, Nparticle, parent.decays_part[k][l]);
			//cout << "Searching for " << pid
			//		<< "   (" << all_particles[pid].name
			//		<< ", " << all_particles[pid].effective_branchratio
			//		<< ")"<< endl;

			// N.B. - using a <set> object will automatically remove duplicates
			//		  and keep pid's in a fixed order
			if ( all_particles[pid].effective_branchratio >= 1.e-12
					or pid == target_particle_id )
				daughter_resonance_indices.insert(pid);
		}
	}

	//cout << endl << "Ended up with n_daughter = "
	//		<< daughter_resonance_indices.size() << " for "
	//		<< all_particles[parent_pid].name << " ("
	//		<< parent_pid << ")" << endl;
	//
	//int i = 0;
	//for ( set<int>::iterator it = daughter_resonance_indices.begin();
	//		it != daughter_resonance_indices.end();
	//		++it)
	//	cout << "parent = " << all_particles[parent_pid].name
	//			<< ": " << i++ << "   " << *it << "   "
	//			<< all_particles[*it].name << "   "
	//			<< all_particles[*it].effective_branchratio << endl;

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

	// if pid is not one of the chosen_resonances,
	// is not the target daughter (pion(+)),
	// is not stable and has a non-zero effective branching ratio,
	// then return a warning message and crash
	if ( result < 0
			and pid != particle_id && all_particles[pid].stable == 0
			and all_particles[pid].effective_branchratio >= 1.e-12 )
	{
		*out << " *** lookup_resonance_idx_from_particle_id(): Particle_id = "
				<< pid << " (" << all_particles[pid].name
				<< ") not found in chosen_resonances!" << endl
				<< " *** br = " << all_particles[pid].effective_branchratio
				<< endl << " *** Can only choose from: " << endl;
		for (int ii = 0; ii < (int)chosen_resonances.size(); ii++)
			*out << all_particles[chosen_resonances[ii]].name
				<< ": pid = " << chosen_resonances[ii] << endl;
		exit(8);
	}

	return (result);
}

void CorrelationFunction::Allocate_resonance_running_sum_vectors()
{
    ssum_vec 	= new double [qspace_cs_slice_length];
    vsum_vec 	= new double [qspace_cs_slice_length];
    zetasum_vec = new double [qspace_cs_slice_length];
    Csum_vec 	= new double [qspace_cs_slice_length];
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

void CorrelationFunction
		::Setup_current_daughters_dN_dypTdpTdphi_moments(int n_daughter)
{
	const int giant_flat_array_size
		= n_pT_pts * n_pphi_pts * n_pY_pts * qxnpts * qynpts * ntrig;
	current_daughters_dN_dypTdpTdphi_moments
		= new double * [n_daughter];

	for (int id = 0; id < n_daughter; ++id)
	{
		current_daughters_dN_dypTdpTdphi_moments[id]
			= new double [giant_flat_array_size];
		for (int i = 0; i < giant_flat_array_size; ++i)
			current_daughters_dN_dypTdpTdphi_moments[id][i] = 0.0;
	}

	return;
}

void CorrelationFunction
		::Cleanup_current_daughters_dN_dypTdpTdphi_moments(int n_daughter)
{
	for (int id = 0; id < n_daughter; ++id)
		delete [] current_daughters_dN_dypTdpTdphi_moments[id];
	delete [] current_daughters_dN_dypTdpTdphi_moments;

	return;
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
void CorrelationFunction::Load_decay_channel_info_nb2(
		int dc_idx, double K_T_local, double K_phi_local, double K_y_local)
{
	Mres 				= current_resonance_mass;
	Gamma 				= current_resonance_Gamma;
	//keeps calculation safe when Gamma == 0
	one_by_Gamma_Mres 	= 1./(Gamma*Mres + 1.e-25);
	//N.B. - no need for hbarc, since this will only
	//       multiply something with GeV^2 units in the end
	mass 				= current_daughter_mass;
	//doesn't depend on target daughter particle,
	// just parent resonance and decay channel
	br 					= current_resonance_direct_br;
	m2 					= current_resonance_decay_masses[0];
	m3 					= current_resonance_decay_masses[1];

	//bad choice of notation throughout...update to a new namespace!
	pT 					= K_T_local;
	current_K_phi 		= K_phi_local;
	p_y 				= K_y_local;

	n_body 				= current_reso_nbody;

	// some particles may decay to particles with
	// more total mass than originally
	// --> broaden with resonance widths
	while ((mass + m2) > Mres)
	{
		Mres 	+= 0.25 * current_resonance_Gamma;
		mass 	-= 0.5 * current_daughter_Gamma;
		m2 		-= 0.5 * current_m2_Gamma;
	}

	// transverse mass
	mT 					= sqrt(mass*mass + pT*pT);

	//set up vectors of points to speed-up integrals...
	double s_loc 		= m2*m2;
	double pstar_loc 	= sqrt( ((Mres+mass)*(Mres+mass) - s_loc)
								*((Mres-mass)*(Mres-mass) - s_loc)
								)/(2.0*Mres);
	//for n_body == 2, doesn't actually use s_loc since the result
	// is just a factor * delta(...) and g(s) just returns factor
	double g_s_loc 		= g(s_loc);
	double Estar_loc 	= sqrt(mass*mass + pstar_loc*pstar_loc);
	double psBmT 		= pstar_loc / mT;
	double DeltaY_loc 	= log(psBmT + sqrt(1.+psBmT*psBmT));

	VEC_n2_s_factor 	= br/(4.*M_PI*pstar_loc);	//==g_s_loc

	for(int iv = 0; iv < n_v_pts; ++iv)
	{
		double v_loc 			= v_pts[iv];
		double P_Y_loc 			= p_y + v_loc*DeltaY_loc;
		double mT_ch_P_Y_p_y 	= mT*cosh(v_loc*DeltaY_loc);
		double x2 				= mT_ch_P_Y_p_y*mT_ch_P_Y_p_y - pT*pT;
		double MTbar_loc 		= Estar_loc*Mres*mT_ch_P_Y_p_y/x2;
		double DeltaMT_loc 		= Mres*pT*sqrt(Estar_loc*Estar_loc - x2)/x2;

		VEC_n2_P_Y[iv] 			= P_Y_loc;
		VEC_n2_v_factor[iv] 	= v_wts[iv]*DeltaY_loc/sqrt(x2);

		for(int izeta = 0; izeta < n_zeta_pts; ++izeta)
		{
			double zeta_loc
				= zeta_pts[izeta];
			double MT_loc
				= MTbar_loc + cos(zeta_loc)*DeltaMT_loc;
			double PT_loc
				= sqrt(MT_loc*MT_loc - Mres*Mres);
			double temp_cos_PPhi_tilde_loc
				= ( mT*MT_loc*cosh(P_Y_loc-p_y)
					- Estar_loc*Mres)/(pT*PT_loc);
			//assume that PPhi_tilde is positive in the next step...
			double temp_sin_PPhi_tilde_loc
				= sqrt(1. - temp_cos_PPhi_tilde_loc
							* temp_cos_PPhi_tilde_loc);
			double PPhi_tilde_loc
				= place_in_range( atan2(temp_sin_PPhi_tilde_loc,
									temp_cos_PPhi_tilde_loc),
									Kphi_min, Kphi_max);

			VEC_n2_zeta_factor[NB2_indexer(iv, izeta)]
				= zeta_wts[izeta]*MT_loc;
			VEC_n2_PPhi_tilde[NB2_indexer(iv, izeta)]
				= place_in_range( K_phi_local + PPhi_tilde_loc,
									Kphi_min, Kphi_max);
			VEC_n2_PPhi_tildeFLIP[NB2_indexer(iv, izeta)]
				= place_in_range( K_phi_local - PPhi_tilde_loc,
									Kphi_min, Kphi_max);

			VEC_n2_PT[NB2_indexer(iv, izeta)]
				= PT_loc;

			//set P^+ components
			VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+0][0]
				= MT_loc * cosh(P_Y_loc);
			VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+0][1]
				= PT_loc * cos(K_phi_local + PPhi_tilde_loc);
			VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+0][2]
				= PT_loc * sin(K_phi_local + PPhi_tilde_loc);
			VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+0][3]
				= MT_loc * sinh(P_Y_loc);
			//set P^- components
			VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+1][0]
				= MT_loc * cosh(P_Y_loc);
			VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+1][1]
				= PT_loc * cos(K_phi_local - PPhi_tilde_loc);
			VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+1][2]
				= PT_loc * sin(K_phi_local - PPhi_tilde_loc);
			VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+1][3]
				= MT_loc * sinh(P_Y_loc);
		}
	}

	return;
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
void CorrelationFunction::Load_decay_channel_info_nb3(
		int dc_idx, double K_T_local, double K_phi_local, double K_y_local)
{
	Mres 				= current_resonance_mass;
	Gamma 				= current_resonance_Gamma;
	one_by_Gamma_Mres 	= 1./(Gamma*Mres + 1.e-25);
	//N.B. - no need for hbarc, since this will only multiply
	// something with GeV^2 units in the end
	mass 				= current_daughter_mass;
	//doesn't depend on target daughter particle,
	//just parent resonance and decay channel
	br 					= current_resonance_direct_br;
	m2 					= current_resonance_decay_masses[0];
	m3 					= current_resonance_decay_masses[1];

	pT 					= K_T_local;
	current_K_phi 		= K_phi_local;
	p_y 				= K_y_local;

	n_body 				= current_reso_nbody;

	mT 					= sqrt(mass*mass + pT*pT);
	double s_min_temp 	= (m2 + m3)*(m2 + m3);
	double s_max_temp 	= (Mres - mass)*(Mres - mass);
	gauss_quadrature(n_s_pts, 1, 0.0, 0.0,
						s_min_temp, s_max_temp, s_pts, s_wts);
	// compute Q function to get denominator
	Qfunc 				= get_Q();

	// s-loop
	for (int is = 0; is < n_s_pts; ++is)
	{

		double s_loc 		= s_pts[is];
		double g_s_loc 		= g(s_loc);
		double pstar_loc 	= sqrt( ((Mres+mass)*(Mres+mass) - s_loc)
									*((Mres-mass)*(Mres-mass) - s_loc)
									)/(2.0*Mres);
		double Estar_loc 	= sqrt(mass*mass + pstar_loc*pstar_loc);
		double psBmT 		= pstar_loc / mT;
		double DeltaY_loc 	= log(psBmT + sqrt(1.+psBmT*psBmT));

		VEC_n3_s_factor[is] = s_wts[is] * g_s_loc;

		// v-loop
		for(int iv = 0; iv < n_v_pts; ++iv)
		{
			double v_loc 			= v_pts[iv];
			double P_Y_loc 			= p_y + v_loc*DeltaY_loc;
			double mT_ch_P_Y_p_y 	= mT*cosh(v_loc*DeltaY_loc);
			double x2 				= mT_ch_P_Y_p_y*mT_ch_P_Y_p_y - pT*pT;
			double MTbar_loc 		= Estar_loc*Mres*mT_ch_P_Y_p_y/x2;
			double DeltaMT_loc 		= Mres*pT*sqrt(Estar_loc*Estar_loc - x2)/x2;

			VEC_n3_P_Y[is * n_v_pts + iv]
				= P_Y_loc;
			VEC_n3_v_factor[is * n_v_pts + iv]
				= v_wts[iv]*DeltaY_loc/sqrt(x2);

			// zeta-loop
			for(int izeta = 0; izeta < n_zeta_pts; ++izeta)
			{
				double zeta_loc
					= zeta_pts[izeta];
				double MT_loc
					= MTbar_loc + cos(zeta_loc)*DeltaMT_loc;
				double PT_loc
					= sqrt(MT_loc*MT_loc - Mres*Mres);
				double temp_cos_PPhi_tilde_loc
					= ( mT*MT_loc*cosh(P_Y_loc-p_y)
						- Estar_loc*Mres )/(pT*PT_loc);
				//assume that PPhi_tilde is positive in next step...
				double temp_sin_PPhi_tilde_loc
					= sqrt( 1. - temp_cos_PPhi_tilde_loc
								* temp_cos_PPhi_tilde_loc );
				double PPhi_tilde_loc
					= place_in_range( atan2(temp_sin_PPhi_tilde_loc,
										temp_cos_PPhi_tilde_loc),
										Kphi_min, Kphi_max);

				VEC_n3_zeta_factor[NB3_indexer(is, iv, izeta)]
					= zeta_wts[izeta]*MT_loc;
				VEC_n3_PPhi_tilde[NB3_indexer(is, iv, izeta)]
					= place_in_range( K_phi_local + PPhi_tilde_loc,
										Kphi_min, Kphi_max);
				VEC_n3_PPhi_tildeFLIP[NB3_indexer(is, iv, izeta)]
					= place_in_range( K_phi_local - PPhi_tilde_loc,
										Kphi_min, Kphi_max);

				VEC_n3_PT[NB3_indexer(is, iv, izeta)]
					= PT_loc;

				//set P^+ components
				VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+0][0]
					= MT_loc * cosh(P_Y_loc);
				VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+0][1]
					= PT_loc * cos(K_phi_local + PPhi_tilde_loc);
				VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+0][2]
					= PT_loc * sin(K_phi_local + PPhi_tilde_loc);
				VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+0][3]
					= MT_loc * sinh(P_Y_loc);
				//set P^- components
				VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+1][0]
					= MT_loc * cosh(P_Y_loc);
				VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+1][1]
					= PT_loc * cos(K_phi_local - PPhi_tilde_loc);
				VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+1][2]
					= PT_loc * sin(K_phi_local - PPhi_tilde_loc);
				VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+1][3]
					= MT_loc * sinh(P_Y_loc);
			}
		}
	}

	return;
}

void CorrelationFunction::Allocate_fleshed_out_CF()
{
	fleshed_out_CF 			= new double ** [new_nqxpts];
	fleshed_out_thermal 	= new double ** [new_nqxpts];
	fleshed_out_crossterm 	= new double ** [new_nqxpts];
	fleshed_out_resonances 	= new double ** [new_nqxpts];
	for (int iqx = 0; iqx < new_nqxpts; ++iqx)
	{
		fleshed_out_CF[iqx] 		= new double * [new_nqypts];
		fleshed_out_thermal[iqx] 	= new double * [new_nqypts];
		fleshed_out_crossterm[iqx] 	= new double * [new_nqypts];
		fleshed_out_resonances[iqx] = new double * [new_nqypts];
		for (int iqy = 0; iqy < new_nqypts; ++iqy)
		{
			fleshed_out_CF[iqx][iqy] 			= new double [new_nqzpts];
			fleshed_out_thermal[iqx][iqy] 		= new double [new_nqzpts];
			fleshed_out_crossterm[iqx][iqy] 	= new double [new_nqzpts];
			fleshed_out_resonances[iqx][iqy] 	= new double [new_nqzpts];
			for (int iqz = 0; iqz < new_nqzpts; ++iqz)
			{
				fleshed_out_CF[iqx][iqy][iqz] 			= 0.0;
				fleshed_out_thermal[iqx][iqy][iqz] 		= 0.0;
				fleshed_out_crossterm[iqx][iqy][iqz] 	= 0.0;
				fleshed_out_resonances[iqx][iqy][iqz] 	= 0.0;
			}
		}
	}

	qx_fleshed_out_pts = new double [new_nqxpts];
	qy_fleshed_out_pts = new double [new_nqypts];
	qz_fleshed_out_pts = new double [new_nqzpts];

	return;
}

void CorrelationFunction::Delete_fleshed_out_CF()
{
	for (int iqx = 0; iqx < new_nqxpts; ++iqx)
	{
		for (int iqy = 0; iqy < new_nqypts; ++iqy)
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

/*void CorrelationFunction::Set_Y_eq_0_Bessel_grids(
		int iqt, int iqz, double * BC_chunk)
{
	const std::complex<double> i(0, 1);
	int n_coeffs = n_alpha_points_PIONS;
	double alpha_min = 0.5, alpha_max = 400.0;
	double coeffs_array[n_alpha_points_PIONS];
	double * alpha_pts = new double [n_alpha_points_PIONS];
	double * x_pts = new double [n_alpha_points_PIONS];

	for (int k = 0; k < n_alpha_points_PIONS; ++k)
	{
		x_pts[k] = - cos( M_PI*(2.*(k+1.) - 1.) / (2.*n_alpha_points_PIONS) );
		alpha_pts[k] = 0.5*(x_pts[k] + 1.0)*(alpha_max - alpha_min) + alpha_min;
	}
	double nums[n_alpha_points_PIONS*n_alpha_points_PIONS];
	double dens[n_alpha_points_PIONS];
	int na = n_alpha_points_PIONS;
	for (int j = 0; j < na; ++j)
	{
		dens[j] = 0.0;
		for (int k = 0; k < na; ++k)
		{
			double Tjk = csf::Tfun(j, x_pts[k]);
			dens[j] += Tjk*Tjk;
			nums[j*na+k] = Tjk;
		}
	}

	double expBesselK0re[na];
	double expBesselK0im[na];
	double expBesselK1re[na];
	double expBesselK1im[na];

	double loc_qt = qt_pts[iqt];
	double loc_qz = qz_pts[iqz];
	current_pY_shift = 0.5 * log(abs((loc_qt+loc_qz + 1.e-100)/(loc_qt-loc_qz + 1.e-100)));

	Stopwatch sw_loop;

	int iBC = 0;
	for (int isurf = 0; isurf < FO_length; ++isurf)
	{
		double tau = (&FOsurf_ptr[isurf])->tau;

		double beta = tau * hbarCm1 * loc_qt;
		double gamma = tau * hbarCm1 * loc_qz;
		double gsq = gamma*gamma;

		for (int ia = 0; ia < na; ++ia)
		{
			double loc_alpha = alpha_pts[ia];
			complex<double> ci0, ci1, ck0, ck1, ci0p, ci1p, ck0p, ck1p;
			complex<double> z0 = loc_alpha - i*beta;
			complex<double> z0sq = z0 * z0;
			complex<double> z = sqrt(z0sq + gsq);
			int errorCode = bessf::cbessik01(z, ci0, ci1, ck0, ck1, ci0p, ci1p, ck0p, ck1p);
			double ea = exp(loc_alpha);

			expBesselK0re[ia] = ea * ck0.real();
			expBesselK0im[ia] = ea * ck0.imag();
			expBesselK1re[ia] = ea * ck1.real();
			expBesselK1im[ia] = ea * ck1.imag();
		}

		//////////////////////////////////
		//exp(x) * K_0(x), real part
		//separate out 0th coefficient for additional factor of 2.0
		coeffs_array[0] = 0.0;
		for (int k = 0; k < na; ++k)
			coeffs_array[0] += 2.0*expBesselK0re[k] * nums[0*na+k];
		BC_chunk[iBC++] = coeffs_array[0] / dens[0];
		for (int j = 1; j < na; ++j)
		{
			coeffs_array[j] = 0.0;
			for (int k = 0; k < na; ++k)
				coeffs_array[j] += expBesselK0re[k] * nums[j*na+k];
			BC_chunk[iBC++] = coeffs_array[j] / dens[j];
		}

		//////////////////////////////////
		//exp(x) * K_0(x), imaginary part
		//separate out 0th coefficient for additional factor of 2.0
		coeffs_array[0] = 0.0;
		for (int k = 0; k < na; ++k)
			coeffs_array[0] += 2.0*expBesselK0im[k] * nums[0*na+k];
		BC_chunk[iBC++] = coeffs_array[0] / dens[0];
		for (int j = 1; j < na; ++j)
		{
			coeffs_array[j] = 0.0;
			for (int k = 0; k < na; ++k)
				coeffs_array[j] += expBesselK0im[k] * nums[j*na+k];
			BC_chunk[iBC++] = coeffs_array[j] / dens[j];
		}

		//////////////////////////////////
		//exp(x) * K_1(x), real part
		//separate out 0th coefficient for additional factor of 2.0
		coeffs_array[0] = 0.0;
		for (int k = 0; k < na; ++k)
			coeffs_array[0] += 2.0*expBesselK1re[k] * nums[0*na+k];
		BC_chunk[iBC++] = coeffs_array[0] / dens[0];
		for (int j = 1; j < na; ++j)
		{
			coeffs_array[j] = 0.0;
			for (int k = 0; k < na; ++k)
				coeffs_array[j] += expBesselK1re[k] * nums[j*na+k];
			BC_chunk[iBC++] = coeffs_array[j] / dens[j];
		}

		//////////////////////////////////
		//exp(x) * K_1(x), imaginary part
		//separate out 0th coefficient for additional factor of 2.0
		coeffs_array[0] = 0.0;
		for (int k = 0; k < na; ++k)
			coeffs_array[0] += 2.0*expBesselK1im[k] * nums[0*na+k];
		BC_chunk[iBC++] = coeffs_array[0] / dens[0];
		for (int j = 1; j < na; ++j)
		{
			coeffs_array[j] = 0.0;
			for (int k = 0; k < na; ++k)
				coeffs_array[j] += expBesselK1im[k] * nums[j*na+k];
			BC_chunk[iBC++] = coeffs_array[j] / dens[j];
		}
	}

	delete [] alpha_pts;
	delete [] x_pts;

	cout << "Finished setting Bessel grids successfully." << endl;
	return;
}
*/


void CorrelationFunction::Set_all_Bessel_grids(
		int iqt, int iqz, int particle_mode /*==0*/)
{
	const std::complex<double> i(0, 1);
	int na 				= n_alpha_points;
	double alpha_min 	= 4.0;
	double alpha_max 	= 200.0;

	int n_coeffs 		= na;
	double coeffs_array[na];
	double * dummy 		= new double [na];
	double * alpha_pts 	= new double [na];
	double * x_pts 		= new double [na];

	for (int k = 0; k < na; ++k)
	{
		x_pts[k] 		= - cos( M_PI*(2.*(k+1.) - 1.) / (2.*na) );
		alpha_pts[k] 	= 0.5 * (x_pts[k] + 1.0)
							* (alpha_max - alpha_min) + alpha_min;
	}

	double nums[na*na];
	double dens[na];

	for (int j = 0; j < na; ++j)
	{
		dens[j] 			= 0.0;
		for (int k = 0; k < na; ++k)
		{
			double Tjk 		= csf::Tfun(j, x_pts[k]);
			dens[j] 	   += Tjk*Tjk;
			nums[j*na+k] 	= Tjk;
		}
	}

	double expBesselK0re[na];
	double expBesselK0im[na];
	double expBesselK1re[na];
	double expBesselK1im[na];

	//initialize
	int HDFcode = Administrate_besselcoeffs_HDF_array(0, particle_mode);
	//if initialization unsuccessful, crash
	if (HDFcode < 0)
		exit(1);

	///////////////////////////////////
	// Loop over pY points
	///////////////////////////////////
	Stopwatch sw_loop;
	double loc_qt 		= qt_pts[iqt];
	double loc_qz 		= qz_pts[iqz];
	current_pY_shift 	= 0.5 * log( abs( (loc_qt+loc_qz + 1.e-100)
										/(loc_qt-loc_qz + 1.e-100) ) );

	double * BC_chunk 	= new double [4 * FO_length * na];
	for (int ipY = 0; ipY < n_pY_pts; ++ipY)
	{
		sw_loop.Reset();
		sw_loop.Start();
		//cout << "Starting loop = "
		//		<< iqt << "   " << iqz << "   " << ipY << endl;

		ch_SP_pY[ipY] 	= cosh(SP_Del_pY[ipY] + current_pY_shift);
		sh_SP_pY[ipY] 	= sinh(SP_Del_pY[ipY] + current_pY_shift);
		double ch_pY 	= ch_SP_pY[ipY];
		double sh_pY 	= sh_SP_pY[ipY];

		int iBC 		= 0;
		for (int isurf = 0; isurf < FO_length; ++isurf)
		{
			double tau 		= (&FOsurf_ptr[isurf])->tau;

			double beta 	= tau * hbarCm1 * ( loc_qt*ch_pY - loc_qz*sh_pY );
			double gamma 	= tau * hbarCm1 * ( loc_qz*ch_pY - loc_qt*sh_pY );
			double gsq 		= gamma*gamma;
			//cout << "beta and gamma = " << beta << "   " << gamma
			//		<< "   " << loc_qt << "   " << loc_qz << "   "
			//		<< ch_pY << "   " << sh_pY << endl;

			for (int ia = 0; ia < na; ++ia)
			{
				complex<double> ci0, ci1, ck0, ck1,
								ci0p, ci1p, ck0p, ck1p;

				double loc_alpha 		= alpha_pts[ia];
				complex<double> z0 		= loc_alpha - i*beta;
				complex<double> z0sq 	= z0 * z0;
				complex<double> z 		= sqrt(z0sq + gsq);
				int errorCode 			= bessf::cbessik01(
											z,
											ci0, ci1, ck0, ck1,
											ci0p, ci1p, ck0p, ck1p);
				if ( isnan(ck0.real())
						or isnan(ck0.imag())
						or isnan(ck1.real())
						or isnan(ck1.imag()) )
				{
					cerr << "WARNING: Obtained NaNs in "
							<< "complex Bessel functions: "
							<< loc_alpha << "   " << beta << "   "
							<< gamma << "   " << z << "   "
							<< ck0 << "   " << ck1 << endl;
				}

				double ea = exp(loc_alpha);

				expBesselK0re[ia] = ea * ck0.real();
				expBesselK0im[ia] = ea * ck0.imag();
				expBesselK1re[ia] = ea * ck1.real();
				expBesselK1im[ia] = ea * ck1.imag();
			}

			//////////////////////////////////
			//exp(x) * K_0(x), real part
			//separate out 0th coefficient for
			//additional factor of 2.0
			coeffs_array[0] = 0.0;
			for (int k = 0; k < na; ++k)
				coeffs_array[0] += 2.0*expBesselK0re[k]
									* nums[0*na+k];
			BC_chunk[iBC++] = coeffs_array[0] / dens[0];

			for (int j = 1; j < na; ++j)
			{
				coeffs_array[j] = 0.0;
				for (int k = 0; k < na; ++k)
					coeffs_array[j] += expBesselK0re[k]
										* nums[j*na+k];
				BC_chunk[iBC++] = coeffs_array[j] / dens[j];
			}

			//////////////////////////////////
			//exp(x) * K_0(x), imaginary part
			//separate out 0th coefficient for additional factor of 2.0
			coeffs_array[0] = 0.0;
			for (int k = 0; k < na; ++k)
				coeffs_array[0] += 2.0*expBesselK0im[k]
									* nums[0*na+k];
			BC_chunk[iBC++] = coeffs_array[0] / dens[0];
			for (int j = 1; j < na; ++j)
			{
				coeffs_array[j] = 0.0;
				for (int k = 0; k < na; ++k)
					coeffs_array[j] += expBesselK0im[k]
										* nums[j*na+k];
				BC_chunk[iBC++] = coeffs_array[j] / dens[j];
			}

			//////////////////////////////////
			//exp(x) * K_1(x), real part
			//separate out 0th coefficient for additional factor of 2.0
			coeffs_array[0] = 0.0;
			for (int k = 0; k < na; ++k)
				coeffs_array[0] += 2.0*expBesselK1re[k]
									* nums[0*na+k];
			BC_chunk[iBC++] = coeffs_array[0] / dens[0];
			for (int j = 1; j < na; ++j)
			{
				coeffs_array[j] = 0.0;
				for (int k = 0; k < na; ++k)
					coeffs_array[j] += expBesselK1re[k]
										* nums[j*na+k];
				BC_chunk[iBC++] = coeffs_array[j] / dens[j];
			}

			//////////////////////////////////
			//exp(x) * K_1(x), imaginary part
			//separate out 0th coefficient for additional factor of 2.0
			coeffs_array[0] = 0.0;
			for (int k = 0; k < na; ++k)
				coeffs_array[0] += 2.0*expBesselK1im[k]
									* nums[0*na+k];
			BC_chunk[iBC++] = coeffs_array[0] / dens[0];
			for (int j = 1; j < na; ++j)
			{
				coeffs_array[j] = 0.0;
				for (int k = 0; k < na; ++k)
					coeffs_array[j] += expBesselK1im[k]
										* nums[j*na+k];
				BC_chunk[iBC++] = coeffs_array[j] / dens[j];
			}
		}

		//finally, store the results
		// 0 - set
		HDFcode = Access_besselcoeffs_in_HDF_array(ipY, 0, BC_chunk,
													particle_mode);

		sw_loop.Stop();
		*out << "Finished in "
				<< sw_loop.printTime()
				<< " seconds." << endl;
	}

	// 2 - close
	HDFcode = Administrate_besselcoeffs_HDF_array(2, particle_mode);
	if (HDFcode < 0)
		exit(2);

	delete [] BC_chunk;
	delete [] alpha_pts;
	delete [] x_pts;
	delete [] dummy;

	*out << "Finished setting Bessel grids successfully." << endl;

	return;
}

//End of file
