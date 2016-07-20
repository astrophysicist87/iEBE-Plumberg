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
#include "Stopwatch.h"
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

double SourceVariances::place_in_range(double phi, double min, double max)
{
	while (phi < min || phi > max)
	{
		if (phi < min) phi += twopi;
		else phi -= twopi;
	}

	return (phi);
}

// ************************************************************
// NEW AND IMPROVED version of Analyze_sourcefunction()
// ************************************************************
void SourceVariances::Analyze_sourcefunction()
{
	Stopwatch BIGsw;
	*global_out_stream_ptr << "Plane angle calculations..." << endl;
	BIGsw.Start();
	double plane_psi = 0.0;
	int iorder = USE_PLANE_PSI_ORDER;
	if (USE_PLANE_PSI_ORDER)
	{
		if (VERBOSE > 0) *global_out_stream_ptr << "Determine nth-order plane angles..." << endl;
		Determine_plane_angle();
		if (VERBOSE > 0) *global_out_stream_ptr << "Analyzing source function w.r.t. " << iorder << " th-order participant plane angle..." << endl;
		if (VERBOSE > 0) *global_out_stream_ptr << "psi = " << plane_psi << endl;
		plane_psi = plane_angle[iorder];
	}
	else
	{
		if (VERBOSE > 0) *global_out_stream_ptr << "Analyzing source function w.r.t. psi_0 = " << plane_psi << endl;
	}
	global_plane_psi = plane_psi;
	BIGsw.Stop();
	*global_out_stream_ptr << "\t ...finished in " << BIGsw.printTime() << " seconds." << endl;

	int decay_channel_loop_cutoff = n_decay_channels;			//loop over direct pions and decay_channels

	if (read_in_all_dN_dypTdpTdphi)	//read in spectra if already calculated
	{
		Read_in_all_dN_dypTdpTdphi();
		if (VERBOSE > 0) *global_out_stream_ptr << "************************************************************" << endl
												<< "* Read in all (thermal) space-time moments!" << endl
												<< "************************************************************" << endl << endl;

	}
	else	// calculate necessary spectra from scratch
	{
		*global_out_stream_ptr << "Setting spacetime moments grid..." << endl;
		BIGsw.Start();
		// ************************************************************
		// loop over decay_channels (idc == 0 corresponds to thermal pions)
		// ************************************************************
		for (int idc = 0; idc <= decay_channel_loop_cutoff; ++idc)				//this is inefficient, but will do the job for now
		{
			// ************************************************************
			// check whether to do this decay channel
			// ************************************************************
			if (idc > 0 && thermal_pions_only)
				break;
			else if (!Do_this_decay_channel(idc))
				continue;
	
			// ************************************************************
			// if so, set decay channel info
			// ************************************************************
			Set_current_particle_info(idc);
	
			// ************************************************************
			// decide whether to recycle old moments or calculate new moments
			// ************************************************************
			Get_spacetime_moments(idc);
		}	//computing all resonances' spacetime moments here first
			//THEN do phase-space integrals
	
		if (VERBOSE > 0) *global_out_stream_ptr << endl << "************************************************************"
												<< endl << "* Computed all (thermal) space-time moments!" << endl
												<< "************************************************************" << endl << endl;
		BIGsw.Stop();
		*global_out_stream_ptr << "\t ...finished all (thermal) space-time moments in " << BIGsw.printTime() << " seconds." << endl;
	}
	if (output_all_dN_dypTdpTdphi)
	{
		Output_all_dN_dypTdpTdphi();
		if (VERBOSE > 0) *global_out_stream_ptr << endl << "************************************************************"
												<< endl << "* Output all (thermal) space-time moments!" << endl
												<< "************************************************************" << endl << endl;
	}
	
	//if (SPACETIME_MOMENTS_ONLY)
	//	return;

	*global_out_stream_ptr << "Computing all phase-space integrals..." << endl;
    BIGsw.Reset();
	BIGsw.Start();
	// ************************************************************
	// Compute feeddown with heaviest resonances first
	// ************************************************************
	for (int idc = 1; idc <= decay_channel_loop_cutoff; ++idc)
	{
		// ************************************************************
		// check whether to do this decay channel
		// ************************************************************
		if (decay_channels[idc-1].resonance_particle_id == target_particle_id || thermal_pions_only)
			break;
		else if (!Do_this_decay_channel(idc))
			continue;

		// ************************************************************
		// if so, set decay channel info
		// ************************************************************
		Set_current_particle_info(idc);

		// ************************************************************
		// begin source variances calculations here...
		// ************************************************************
		Allocate_decay_channel_info();				// allocate needed memory
		for (int idc_DI = 0; idc_DI < current_reso_nbody; ++idc_DI)
		{
			int daughter_resonance_particle_id = -1;
			if (!Do_this_daughter_particle(idc, idc_DI, &daughter_resonance_particle_id))
				continue;
			Set_current_daughter_info(idc, idc_DI);
			Do_resonance_integrals(current_resonance_particle_id, daughter_resonance_particle_id, idc);
		}
		Delete_decay_channel_info();				// free up memory
	}											// END of decay channel loop

	// *************************************************************
	// now get HBT radii from source integrals and Fourier transform
	// *************************************************************

	//fill out grid of source_variances
	Set_source_variances_grid();

	//compute HBT radii on grid
	for(int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		double m_perp = sqrt(SPinterp_pT[ipt]*SPinterp_pT[ipt] + particle_mass*particle_mass);
		beta_perp = SPinterp_pT[ipt]/(m_perp*cosh(K_y));

		for(int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			Calculate_R2_side(ipt, ipphi);
			Calculate_R2_out(ipt, ipphi);
			Calculate_R2_long(ipt, ipphi);
			Calculate_R2_outside(ipt, ipphi);
			Calculate_R2_sidelong(ipt, ipphi);
			Calculate_R2_outlong(ipt, ipphi);
		}
	}

	//then interpolate to desired pair-momentum values and Fourier transform
	for(int iKT = 0; iKT < n_localp_T; ++iKT)
	{
		if (abs(K_T[iKT]) < 1.e-10)
			continue;

		for(int iKphi = 0; iKphi < n_localp_phi; ++iKphi)
		{
			Interpolate_HBT_radii(iKT, iKphi);
			R2_Fourier_transform(iKT, plane_psi);
		}
	}

	/*double tmp_R2_side[n_localp_T];
	double tmp_R2_out[n_localp_T];
	for (int iKT = 0; iKT < n_localp_T; ++iKT)
	{
		tmp_R2_side[iKT] = R2_side_C[iKT][0];
		tmp_R2_out[iKT] = R2_out_C[iKT][0];
	}*/

   return;
}

bool SourceVariances::Do_this_decay_channel(int dc_idx)
{
	string local_name = "Thermal pion(+)";
	if (dc_idx == 0)
	{
		if (VERBOSE > 0) *global_out_stream_ptr << endl << local_name << ": doing this one." << endl;
		return true;
	}
	else
	{
		local_name = decay_channels[dc_idx-1].resonance_name;
		Get_current_decay_string(dc_idx, &current_decay_channel_string);
	}
	if (decay_channels[dc_idx-1].include_channel)
	{
		;
	}
	else
	{
		if (VERBOSE > 0) *global_out_stream_ptr << endl << local_name << ": skipping decay " << current_decay_channel_string << "." << endl;
	}

	return (decay_channels[dc_idx-1].include_channel);
}

// ************************************************************
// Checks whether to do daughter particle for given decay channel
// ************************************************************
bool SourceVariances::Do_this_daughter_particle(int dc_idx, int daughter_idx, int * daughter_resonance_pid)
{
	// assume dc_idx > 0
	string local_name = decay_channels[dc_idx-1].resonance_name;

	// look up daughter particle info
	int temp_monval = decay_channels[dc_idx-1].resonance_decay_monvals[daughter_idx];

	if (temp_monval == 0)
		return false;

	int temp_ID = lookup_particle_id_from_monval(all_particles, Nparticle, temp_monval);
	//*daughter_resonance_idx = lookup_resonance_idx_from_particle_id(temp_ID) + 1;
	*daughter_resonance_pid = temp_ID;
	// if daughter was found in chosen_resonances or is pion(+), this is correct
	particle_info temp_daughter = all_particles[temp_ID];

	if (*daughter_resonance_pid < 0 && temp_daughter.monval != particle_monval && temp_daughter.effective_branchratio >= 1.e-12)
		*global_out_stream_ptr << "Couldn't find " << temp_daughter.name << " in chosen_resonances!  Results are probably not reliable..." << endl;

	//bool daughter_does_not_contribute = ( (temp_daughter.stable == 1 || temp_daughter.effective_branchratio < 1.e-12) && temp_daughter.monval != particle_monval );
	bool daughter_does_not_contribute = ( (temp_daughter.decays_Npart[0] == 1 || temp_daughter.effective_branchratio < 1.e-12) && temp_daughter.monval != particle_monval );

	// if daughter particle gives no contribution to final pion spectra
	if (daughter_does_not_contribute)
	{
		if (VERBOSE > 0) *global_out_stream_ptr << "\t * " << local_name << ": in decay " << current_decay_channel_string << ", skipping " << temp_daughter.name
												<< " (daughter_resonance_pid = " << *daughter_resonance_pid << ")." << endl;
		return false;
	}
	else
	{
		if (VERBOSE > 0) *global_out_stream_ptr << "\t * " << local_name << ": in decay " << current_decay_channel_string << ", doing " << temp_daughter.name
												<< " (daughter_resonance_pid = " << *daughter_resonance_pid << ")." << endl;
		return true;
	}
}

void SourceVariances::Set_current_particle_info(int dc_idx)
{
	if (dc_idx == 0)
	{
		muRES = particle_mu;
		signRES = particle_sign;
		gRES = particle_gspin;
		current_resonance_particle_id = target_particle_id;
		
		return;
	}
	else
	{
		// assume dc_idx > 0
		string local_name = decay_channels[dc_idx-1].resonance_name;

		if (VERBOSE > 0) *global_out_stream_ptr << local_name << ": doing decay " << current_decay_channel_string << "." << endl
			<< "\t * " << local_name << ": setting information for this decay channel..." << endl;

		if (dc_idx > 1)
		{
			//cerr << "Setting previous decay channel information for dc_idx = " << dc_idx << endl;
			previous_resonance_particle_id = current_resonance_particle_id;		//for look-up in all_particles
			previous_decay_channel_idx = current_decay_channel_idx;			//different for each decay channel
			previous_resonance_idx = current_resonance_idx;				//different for each decay channel
			previous_resonance_mass = current_resonance_mass;
			previous_resonance_Gamma = current_resonance_Gamma;
			previous_resonance_total_br = current_resonance_total_br;
			previous_resonance_direct_br = current_resonance_direct_br;
			previous_reso_nbody = current_reso_nbody;
		}
		//cerr << "Setting current decay channel information for dc_idx = " << dc_idx << endl;
		current_decay_channel_idx = dc_idx;
		current_resonance_idx = decay_channels[dc_idx-1].resonance_idx;
		current_resonance_particle_id = decay_channels[dc_idx-1].resonance_particle_id;
		current_resonance_mass = decay_channels[dc_idx-1].resonance_mass;
		current_resonance_Gamma = decay_channels[dc_idx-1].resonance_Gamma;
		current_resonance_total_br = decay_channels[dc_idx-1].resonance_total_br;
		current_resonance_direct_br = decay_channels[dc_idx-1].resonance_direct_br;
		current_reso_nbody = decay_channels[dc_idx-1].nbody;
		
		// might want to rename these for notational consistency...
		muRES = decay_channels[dc_idx-1].resonance_mu;
		signRES = decay_channels[dc_idx-1].resonance_sign;
		gRES = decay_channels[dc_idx-1].resonance_gspin;
		
		if (dc_idx > 1)
		{
			int similar_particle_idx = -1;
			int temp_reso_idx = decay_channels[dc_idx-1].resonance_idx;
			
			if ( current_resonance_particle_id == previous_resonance_particle_id )
			{
				//previous resonance is the same as this one...
				recycle_previous_moments = true;
				recycle_similar_moments = false;
			}
			else if ( Search_for_similar_particle( temp_reso_idx, &similar_particle_idx ) )
			{
				//previous resonance is NOT the same as this one BUT this one is sufficiently similar to some preceding one...
				recycle_previous_moments = false;
				recycle_similar_moments = true;
				reso_particle_id_of_moments_to_recycle = chosen_resonances[similar_particle_idx];
			}
			else
			{
				recycle_previous_moments = false;
				recycle_similar_moments = false;
				reso_particle_id_of_moments_to_recycle = -1;	//guarantees it won't be used spuriously
			}
		}
	}
	
	return;
}

void SourceVariances::Set_current_daughter_info(int dc_idx, int daughter_idx)
{
	if (dc_idx > 1)
	{
		previous_resonance_particle_id = current_resonance_particle_id;		//for look-up in all_particles
		previous_decay_channel_idx = current_decay_channel_idx;			//different for each decay channel
		previous_resonance_idx = current_resonance_idx;
		previous_resonance_mass = current_resonance_mass;
		previous_resonance_Gamma = current_resonance_Gamma;
		previous_m2_Gamma = current_m2_Gamma;
		previous_m3_Gamma = current_m3_Gamma;
		previous_resonance_total_br = current_resonance_total_br;
		previous_resonance_direct_br = current_resonance_direct_br;
		previous_reso_nbody = current_reso_nbody;
		previous_daughter_mass = current_daughter_mass;
		previous_daughter_Gamma = current_daughter_Gamma;
	}
	current_decay_channel_idx = dc_idx;
	current_resonance_idx = decay_channels[dc_idx-1].resonance_idx;
	current_resonance_particle_id = decay_channels[dc_idx-1].resonance_particle_id;
	current_resonance_mass = decay_channels[dc_idx-1].resonance_mass;
	current_resonance_Gamma = decay_channels[dc_idx-1].resonance_Gamma;
	current_resonance_total_br = decay_channels[dc_idx-1].resonance_total_br;
	current_resonance_direct_br = decay_channels[dc_idx-1].resonance_direct_br;
	current_reso_nbody = decay_channels[dc_idx-1].nbody;
	current_daughter_mass = decay_channels[dc_idx-1].resonance_decay_masses[daughter_idx];
	current_daughter_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[daughter_idx];

	// might want to rename these for notational consistency...
	muRES = decay_channels[dc_idx-1].resonance_mu;
	signRES = decay_channels[dc_idx-1].resonance_sign;
	gRES = decay_channels[dc_idx-1].resonance_gspin;

	// set non-daughter decay masses for computing contributions to spectra of daughter
	double m2ex = 0.0, m3ex = 0.0, m4ex = 0.0;
	switch(current_reso_nbody)
	{
		case 1:
			break;
		case 2:
			current_resonance_decay_masses[1] = 0.0;
			current_m3_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[0];
			if (daughter_idx == 0)
			{
				current_resonance_decay_masses[0] = decay_channels[dc_idx-1].resonance_decay_masses[1];
				current_m2_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[1];
			}
			else
			{
				current_resonance_decay_masses[0] = decay_channels[dc_idx-1].resonance_decay_masses[0];
				current_m2_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[0];
			}
			break;
		case 3:
			if (daughter_idx == 0)
			{
				current_resonance_decay_masses[0] = decay_channels[dc_idx-1].resonance_decay_masses[1];
				current_resonance_decay_masses[1] = decay_channels[dc_idx-1].resonance_decay_masses[2];
				current_m2_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[1];
				current_m3_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[2];
			}
			else if (daughter_idx == 1)
			{
				current_resonance_decay_masses[0] = decay_channels[dc_idx-1].resonance_decay_masses[0];
				current_resonance_decay_masses[1] = decay_channels[dc_idx-1].resonance_decay_masses[2];
				current_m2_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[0];
				current_m3_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[2];
			}
			else
			{
				current_resonance_decay_masses[0] = decay_channels[dc_idx-1].resonance_decay_masses[0];
				current_resonance_decay_masses[1] = decay_channels[dc_idx-1].resonance_decay_masses[1];
				current_m2_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[0];
				current_m3_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[1];
			}
			break;
		case 4:
			if (daughter_idx == 0)
			{
				m2ex = decay_channels[dc_idx-1].resonance_decay_masses[1];
				m3ex = decay_channels[dc_idx-1].resonance_decay_masses[2];
				m4ex = decay_channels[dc_idx-1].resonance_decay_masses[3];
			}
			else if (daughter_idx == 1)
			{
				m2ex = decay_channels[dc_idx-1].resonance_decay_masses[0];
				m3ex = decay_channels[dc_idx-1].resonance_decay_masses[2];
				m4ex = decay_channels[dc_idx-1].resonance_decay_masses[3];
			}
			else if (daughter_idx == 2)
			{
				m2ex = decay_channels[dc_idx-1].resonance_decay_masses[0];
				m3ex = decay_channels[dc_idx-1].resonance_decay_masses[1];
				m4ex = decay_channels[dc_idx-1].resonance_decay_masses[3];
			}
			else
			{
				m2ex = decay_channels[dc_idx-1].resonance_decay_masses[0];
				m3ex = decay_channels[dc_idx-1].resonance_decay_masses[1];
				m4ex = decay_channels[dc_idx-1].resonance_decay_masses[2];
			}
			current_resonance_decay_masses[0] = m2ex;
			current_resonance_decay_masses[1] = 0.5 * (m3ex + m4ex + current_resonance_mass - current_daughter_mass - m2ex);
			// approximation obtained from earlier resonances code
			//*global_out_stream_ptr << "Current decay " << current_decay_channel_string << ", br = " << current_resonance_direct_br
			//						<< ": {m2ex, m3ex, m4ex, m3eff} = {"
			//						<< m2ex << ", " << m3ex << ", " << m4ex << ", " << current_resonance_decay_masses[1] << "}" << endl;
			break;
		default:
			cerr << "Set_current_daughter_info(): shouldn't have ended up here, bad value of current_reso_nbody!" << endl;
			exit(1);
	}
}

bool SourceVariances::Search_for_similar_particle(int reso_idx, int * result)
{
	// for the timebeing, just search from beginning of decay_channels until similar particle is found;
	// should be more careful, since could lead to small numerical discrepancies if similar particle was
	// already recycled by some other (dissimilar) particle, but ignore this complication for now...
	*result = -1;
	
	for (int local_ir = 0; local_ir < reso_idx; ++local_ir)
	{// only need to search decay_channels that have already been calculated
		if (particles_are_the_same(local_ir, reso_idx))
		{
			*result = local_ir;
			break;
		}
	}
	
	return (*result >= 0);
}

//**********************************************************************************************
bool SourceVariances::particles_are_the_same(int reso_idx1, int reso_idx2)
{
	int icr1 = chosen_resonances[reso_idx1];
	int icr2 = chosen_resonances[reso_idx2];
	//int icr1 = reso_idx1;
	//int icr2 = reso_idx2;
	if (all_particles[icr1].sign != all_particles[icr2].sign)
		return false;
	if (abs(all_particles[icr1].mass-all_particles[icr2].mass) / (all_particles[icr2].mass+1.e-30) > PARTICLE_DIFF_TOLERANCE)
		return false;
	//assume chemical potential mu is constant over entire FO surface
	double chem1 = all_particles[icr1].mu, chem2 = all_particles[icr2].mu;
	if (2.*abs(chem1 - chem2)/(chem1 + chem2 + 1.e-30) > PARTICLE_DIFF_TOLERANCE)
		return false;

	return true;
}

void SourceVariances::Recycle_spacetime_moments()
{
	for(int wfi = 0; wfi < n_weighting_functions; ++wfi)
	for(int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for(int iphi = 0; iphi < n_interp_pphi_pts; ++iphi)
	{
		dN_dypTdpTdphi_moments[current_resonance_particle_id][wfi][ipt][iphi] = dN_dypTdpTdphi_moments[reso_particle_id_of_moments_to_recycle][wfi][ipt][iphi];
		ln_dN_dypTdpTdphi_moments[current_resonance_particle_id][wfi][ipt][iphi] = ln_dN_dypTdpTdphi_moments[reso_particle_id_of_moments_to_recycle][wfi][ipt][iphi];
		sign_of_dN_dypTdpTdphi_moments[current_resonance_particle_id][wfi][ipt][iphi] = sign_of_dN_dypTdpTdphi_moments[reso_particle_id_of_moments_to_recycle][wfi][ipt][iphi];
	}

	return;
}

void SourceVariances::Get_spacetime_moments(int dc_idx)
{
//**************************************************************
//Set resonance name
//**************************************************************
	string local_name = "Thermal pion(+)";
	if (dc_idx > 0)
		local_name = decay_channels[dc_idx-1].resonance_name;
//**************************************************************
//Decide what to do with this resonance / decay channel
//**************************************************************
	if (recycle_previous_moments && dc_idx > 1)	// similar (but different) earlier resonance
	{
		if (VERBOSE > 0) *global_out_stream_ptr << local_name
			<< ": new parent resonance (" << decay_channels[current_decay_channel_idx-1].resonance_name << ", dc_idx = " << current_decay_channel_idx
			<< " of " << n_decay_channels << ") same as preceding parent resonance \n\t\t--> reusing old dN_dypTdpTdphi_moments!" << endl;
	}
	else if (recycle_similar_moments && dc_idx > 1)	// similar (but different) earlier resonance
	{
		if (VERBOSE > 0) *global_out_stream_ptr << local_name
			<< ": new parent resonance (" << decay_channels[current_decay_channel_idx-1].resonance_name << ", dc_idx = " << current_decay_channel_idx
			<< " of " << n_decay_channels << ") sufficiently close to preceding parent resonance (" << all_particles[reso_particle_id_of_moments_to_recycle].name
			<< ", reso_particle_id = " << reso_particle_id_of_moments_to_recycle << ") \n\t\t--> reusing old dN_dypTdpTdphi_moments!" << endl;
		Recycle_spacetime_moments();
	}
	else
	{
		if (dc_idx == 0)	//if it's thermal pions
		{
			if (VERBOSE > 0) *global_out_stream_ptr << "  --> Computing dN_dypTdpTdphi_moments for thermal pion(+)!" << endl;
		}
		else if (dc_idx == 1)	//if it's the first resonance
		{
			if (VERBOSE > 0) *global_out_stream_ptr << "  --> Computing dN_dypTdpTdphi_moments for " << local_name << endl;
		}
		else			//if it's a later resonance
		{
			if (VERBOSE > 0 && !recycle_previous_moments && !recycle_similar_moments) *global_out_stream_ptr << local_name
				<< ": new parent resonance (" << decay_channels[current_decay_channel_idx-1].resonance_name << ", dc_idx = " << current_decay_channel_idx
				<< " of " << n_decay_channels << ") dissimilar from all preceding decay_channels \n\t\t--> calculating new dN_dypTdpTdphi_moments!" << endl;
			else
			{
				cerr << "You shouldn't have ended up here!" << endl;
				exit(1);
			}
		}
		Set_dN_dypTdpTdphi_moments(current_resonance_particle_id);
	}
//**************************************************************
//Spacetime moments now set
//**************************************************************
	return;
}

void SourceVariances::Set_dN_dypTdpTdphi_moments(int local_pid)
{
	double localmass = all_particles[local_pid].mass;
	string local_name = "Thermal pion(+)";
	if (local_pid != target_particle_id)
		local_name = all_particles[local_pid].name;

	for(int i=0; i<eta_s_npts; ++i)
	{
		double local_eta_s = eta_s[i];
		double local_cosh = cosh(SP_p_y - local_eta_s);
		double local_sinh = sinh(SP_p_y - local_eta_s);

		for(int ipt=0; ipt<n_interp_pT_pts; ++ipt)
		{
			double mT = sqrt(localmass*localmass + SPinterp_pT[ipt]*SPinterp_pT[ipt]);
			SPinterp_p0[ipt][i] = mT*local_cosh;
			SPinterp_pz[ipt][i] = mT*local_sinh;
		}
	}
	//if (local_pid != target_particle_id)						//skip thermal contributions for timebeing...
		Cal_dN_dypTdpTdphi_with_weights(local_pid);

	return;
}

void SourceVariances::Cal_dN_dypTdpTdphi_with_weights(int local_pid)
{
	double * local_temp_moments = new double [n_weighting_functions];
	double *** temp_moments_array = new double ** [n_interp_pT_pts];
	for (int ipt = 0; ipt < n_interp_pT_pts; ipt++)
	{
		temp_moments_array[ipt] = new double * [n_interp_pphi_pts];
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ipphi++)
		{
			temp_moments_array[ipt][ipphi] = new double [n_weighting_functions];
			for (int wfi = 0; wfi < n_weighting_functions; wfi++)
				temp_moments_array[ipt][ipphi][wfi] = 0.0;
		}
	}

	// set particle information
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;
	double deltaf_prefactor = 0.;
	if (INCLUDE_DELTA_F)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 2.0;

	for(int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		double pT = SPinterp_pT[ipt];
		double * p0_pTslice = SPinterp_p0[ipt];
		double * pz_pTslice = SPinterp_pz[ipt];

		for(int iphi = 0; iphi < n_interp_pphi_pts; ++iphi)
		{
			double sin_pphi = sin_SPinterp_pphi[iphi];
			double cos_pphi = cos_SPinterp_pphi[iphi];
			double px = pT*cos_pphi;
			double py = pT*sin_pphi;

			for (int wfi = 0; wfi < n_weighting_functions; wfi++)
				local_temp_moments[wfi] = 0.0;

			for(int isurf=0; isurf<FO_length; ++isurf)
			{
				FO_surf*surf = &FOsurf_ptr[isurf];

				double tau = surf->tau;

				double vx = surf->vx;
				double vy = surf->vy;
				double gammaT = surf->gammaT;

				double da0 = surf->da0;
				double da1 = surf->da1;
				double da2 = surf->da2;

				double pi00 = surf->pi00;
				double pi01 = surf->pi01;
				double pi02 = surf->pi02;
				double pi11 = surf->pi11;
				double pi12 = surf->pi12;
				double pi22 = surf->pi22;
				double pi33 = surf->pi33;

				double x = surf->xpt;
				double y = surf->ypt;

				for(int ieta=0; ieta < eta_s_npts; ++ieta)
				{
					double p0 = p0_pTslice[ieta];
					double pz = pz_pTslice[ieta];
					double f0;

					//now get distribution function, emission function, etc.
					f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
	
					//viscous corrections
					double deltaf = 0.;
					//if (INCLUDE_DELTA_F)
						deltaf = deltaf_prefactor * (1. - sign*f0)
								* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

					//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
					double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

					//ignore points where delta f is large or emission function goes negative from pdsigma
					if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
					{
						S_p = 0.0;
						continue;
					}

					double S_p_withweight = S_p*tau*eta_s_weight[ieta];
					local_temp_moments[0] += eta_even_factor*S_p_withweight;					//<1>
					//if (INCLUDE_SOURCE_VARIANCES)
					//{
						double t = tau*ch_eta_s[ieta];
						double z = tau*sh_eta_s[ieta];
						local_temp_moments[1] += eta_even_factor*S_p_withweight*x;
						local_temp_moments[2] += eta_even_factor*S_p_withweight*x*x;
						local_temp_moments[3] += eta_even_factor*S_p_withweight*y;
						local_temp_moments[4] += eta_even_factor*S_p_withweight*y*y;
						local_temp_moments[5] += eta_odd_factor*S_p_withweight*z;
						local_temp_moments[6] += eta_even_factor*S_p_withweight*z*z;
						local_temp_moments[7] += eta_even_factor*S_p_withweight*t;
						local_temp_moments[8] += eta_even_factor*S_p_withweight*t*t;
						local_temp_moments[9] += eta_even_factor*S_p_withweight*x*y;
						local_temp_moments[10] += eta_odd_factor*S_p_withweight*x*z;
						local_temp_moments[11] += eta_even_factor*S_p_withweight*x*t;
						local_temp_moments[12] += eta_odd_factor*S_p_withweight*y*z;
						local_temp_moments[13] += eta_even_factor*S_p_withweight*y*t;
						local_temp_moments[14] += eta_odd_factor*S_p_withweight*z*t;
					//}
				}
			}
			for (int wfi = 0; wfi < n_weighting_functions; wfi++)
				temp_moments_array[ipt][iphi][wfi] = local_temp_moments[wfi];
		}
	}

	//set log of dN_dypTdpTdphi_moments...
	for(int wfi = 0; wfi < n_weighting_functions; ++wfi)
	for(int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for(int iphi = 0; iphi < n_interp_pphi_pts; ++iphi)
	{
		double temp = temp_moments_array[ipt][iphi][wfi];
		dN_dypTdpTdphi_moments[local_pid][wfi][ipt][iphi] = temp;
		ln_dN_dypTdpTdphi_moments[local_pid][wfi][ipt][iphi] = log(abs(temp) + 1.e-100);
		sign_of_dN_dypTdpTdphi_moments[local_pid][wfi][ipt][iphi] = sgn(temp);
	}

	for(int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		for(int iphi = 0; iphi < n_interp_pphi_pts; ++iphi)
			delete [] temp_moments_array[ipt][iphi];
		delete [] temp_moments_array[ipt];
	}
	delete [] temp_moments_array;
	delete [] local_temp_moments;

	return;
}

void SourceVariances::Determine_plane_angle()
{
	double localmass = particle_mass;
	double* mT = new double [n_SP_pT];
	double** px = new double* [n_SP_pT];
	double** py = new double* [n_SP_pT];
	double** p0 = new double* [n_SP_pT];
	double** pz = new double* [n_SP_pT];
	for(int ipt=0; ipt<n_SP_pT; ++ipt)
	{
		px[ipt] = new double [n_SP_pphi];
		py[ipt] = new double [n_SP_pphi];
		p0[ipt] = new double [eta_s_npts];
		pz[ipt] = new double [eta_s_npts];
	}
   
	for(int ipt=0; ipt<n_SP_pT; ++ipt)
		mT[ipt] = sqrt(localmass*localmass + SP_pT[ipt]*SP_pT[ipt]);
	for(int iphi = 0; iphi<n_SP_pphi; ++iphi)
	{
		double cos_phi = cos(SP_pphi[iphi]);
		double sin_phi = sin(SP_pphi[iphi]);
		for(int ipt=0; ipt<n_SP_pT; ++ipt)
		{
			px[ipt][iphi] = SP_pT[ipt]*cos_phi;
			py[ipt][iphi] = SP_pT[ipt]*sin_phi;
		}
	}

	for(int i=0; i<eta_s_npts; ++i)
	{
		double local_eta_s = eta_s[i];
		double local_cosh = cosh(SP_p_y - local_eta_s);
		double local_sinh = sinh(SP_p_y - local_eta_s);
		for(int ipt=0; ipt<n_SP_pT; ++ipt)
		{
			p0[ipt][i] = mT[ipt]*local_cosh;
			pz[ipt][i] = mT[ipt]*local_sinh;
		}
	}

	Cal_dN_dypTdpTdphi(p0, px, py, pz);

	for(int ipt=0; ipt<n_SP_pT; ++ipt)
  	{
		for(int iphi=0; iphi<n_SP_pphi; ++iphi)
		{
			dN_dydphi[iphi] += dN_dypTdpTdphi[ipt][iphi]*SP_pT[ipt]*SP_pT_weight[ipt];
			pTdN_dydphi[iphi] += dN_dypTdpTdphi[ipt][iphi]*SP_pT[ipt]*SP_pT[ipt]*SP_pT_weight[ipt];
			dN_dypTdpT[ipt] += dN_dypTdpTdphi[ipt][iphi]*SP_pphi_weight[iphi];
		}
	}
   	double norm = 0.0e0;
   	for(int iphi=0; iphi<n_SP_pphi; ++iphi)
		norm += dN_dydphi[iphi]*SP_pphi_weight[iphi];
   	for(int iorder=0; iorder < n_order; iorder++)
   	{
		double cosine = 0.0e0;
		double sine = 0.0e0;
		for(int iphi=0; iphi<n_SP_pphi; ++iphi)
		{
			cosine += dN_dydphi[iphi]*cos(iorder*SP_pphi[iphi])*SP_pphi_weight[iphi];
			sine += dN_dydphi[iphi]*sin(iorder*SP_pphi[iphi])*SP_pphi_weight[iphi];
			//get pT-differential v_n here
			for(int ipt=0; ipt<n_SP_pT; ++ipt)
			{
				cosine_iorder[ipt][iorder] += dN_dypTdpTdphi[ipt][iphi]*cos(iorder*SP_pphi[iphi])*SP_pphi_weight[iphi];
				sine_iorder[ipt][iorder] += dN_dypTdpTdphi[ipt][iphi]*sin(iorder*SP_pphi[iphi])*SP_pphi_weight[iphi];
			}
		}
		for(int ipt=0; ipt<n_SP_pT; ++ipt)
		{
			cosine_iorder[ipt][iorder] /= dN_dypTdpT[ipt];
			sine_iorder[ipt][iorder] /= dN_dypTdpT[ipt];
		}
		cosine = cosine/norm;
		sine = sine/norm;
		if( sqrt(sine*sine + cosine*cosine) < 1e-8)
			plane_angle[iorder] = 0.0e0;
		else
			plane_angle[iorder] = atan2(sine, cosine)/double(iorder);
   	}
	
	for(int ipt=0; ipt<n_SP_pT; ++ipt)
		dN_dypTdpT[ipt] /= (2.*M_PI);
	
	mean_pT = 0.;
	for(int iphi=0; iphi<n_SP_pphi; ++iphi)
		mean_pT += pTdN_dydphi[iphi]*SP_pphi_weight[iphi];
	mean_pT /= norm;
	plane_angle[0] = norm;

	delete[] mT;
	for(int ipt=0; ipt<n_SP_pT; ++ipt)
	{
		delete [] px[ipt];
		delete [] py[ipt];
		delete [] p0[ipt];
		delete [] pz[ipt];
	}
	delete [] px;
	delete [] py;
	delete [] p0;
	delete [] pz;

	return;
}

void SourceVariances::Load_decay_channel_info(int dc_idx, double K_T_local, double K_phi_local)
{
	Mres = current_resonance_mass;
	Gamma = current_resonance_Gamma;
	double one_by_Gamma_Mres = hbarC/(Gamma*Mres);
	mass = current_daughter_mass;
	br = current_resonance_direct_br;	//doesn't depend on target daughter particle, just parent resonance and decay channel
	m2 = current_resonance_decay_masses[0];
	m3 = current_resonance_decay_masses[1];
	pT = K_T_local;
	current_K_phi = K_phi_local;
	n_body = current_reso_nbody;
	if (n_body == 2)
	{
		// some particles may decay to particles with more total mass than originally
		// --> broaden with resonance widths
		while ((mass + m2) > Mres)
		{
			Mres += 0.25 * current_resonance_Gamma;
			mass -= 0.5 * current_daughter_Gamma;
			m2 -= 0.5 * current_m2_Gamma;
		}

		mT = sqrt(mass*mass + pT*pT);

		//set up vectors of points to speed-up integrals...
		double s_loc = m2*m2;
		VEC_n2_spt = s_loc;
		double pstar_loc = sqrt( ((Mres+mass)*(Mres+mass) - s_loc)*((Mres-mass)*(Mres-mass) - s_loc) )/(2.0*Mres);
		VEC_n2_pstar = pstar_loc;
		check_for_NaNs("pstar_loc", pstar_loc, *global_out_stream_ptr);
		double g_s_loc = g(s_loc);	//for n_body == 2, doesn't actually use s_loc since result is just a factor * delta(...); just returns factor
		VEC_n2_s_factor = br/(4.*M_PI*VEC_n2_pstar);	//==g_s_loc
		double Estar_loc = sqrt(mass*mass + pstar_loc*pstar_loc);
		VEC_n2_Estar = Estar_loc;
		double psBmT = pstar_loc / mT;
		VEC_n2_psBmT = psBmT;
		double DeltaY_loc = log(psBmT + sqrt(1.+psBmT*psBmT));
		VEC_n2_DeltaY = DeltaY_loc;
		p_y = 0.0;
		VEC_n2_Yp = p_y + DeltaY_loc;
		VEC_n2_Ym = p_y - DeltaY_loc;
		check_for_NaNs("DeltaY_loc", DeltaY_loc, *global_out_stream_ptr);
		for(int iv = 0; iv < n_v_pts; ++iv)
		{
			double v_loc = v_pts[iv];
			double P_Y_loc = p_y + v_loc*DeltaY_loc;
			VEC_n2_P_Y[iv] = P_Y_loc;
			double mT_ch_P_Y_p_y = mT*cosh(v_loc*DeltaY_loc);
			double x2 = mT_ch_P_Y_p_y*mT_ch_P_Y_p_y - pT*pT;
			VEC_n2_v_factor[iv] = v_wts[iv]*DeltaY_loc/sqrt(x2);
			double MTbar_loc = Estar_loc*Mres*mT_ch_P_Y_p_y/x2;
			VEC_n2_MTbar[iv] = MTbar_loc;
			double DeltaMT_loc = Mres*pT*sqrt(Estar_loc*Estar_loc - x2)/x2;
			VEC_n2_DeltaMT[iv] = DeltaMT_loc;
			VEC_n2_MTp[iv] = MTbar_loc + DeltaMT_loc;
			VEC_n2_MTm[iv] = MTbar_loc - DeltaMT_loc;
			check_for_NaNs("MTbar_loc", MTbar_loc, *global_out_stream_ptr);
			check_for_NaNs("DeltaMT_loc", DeltaMT_loc, *global_out_stream_ptr);

			for(int izeta = 0; izeta < n_zeta_pts; ++izeta)
			{
				double zeta_loc = zeta_pts[izeta];
				double MT_loc = MTbar_loc + cos(zeta_loc)*DeltaMT_loc;
				VEC_n2_MT[iv][izeta] = MT_loc;
				VEC_n2_zeta_factor[iv][izeta] = zeta_wts[izeta]*MT_loc;
				///////////////////////////////////////////////////////////////
				// Get \tilde{Phi}...
				double PT_loc = sqrt(MT_loc*MT_loc - Mres*Mres);
				double temp_cos_PPhi_tilde_loc = (mT*MT_loc*cosh(P_Y_loc-p_y) - Estar_loc*Mres)/(pT*PT_loc);
				//assume that PPhi_tilde is +ve in next step...
				double temp_sin_PPhi_tilde_loc = sqrt(1. - temp_cos_PPhi_tilde_loc*temp_cos_PPhi_tilde_loc);
				double PPhi_tilde_loc = place_in_range( atan2(temp_sin_PPhi_tilde_loc, temp_cos_PPhi_tilde_loc), interp_pphi_min, interp_pphi_max);
				VEC_n2_PPhi_tilde[iv][izeta] = place_in_range( K_phi_local + PPhi_tilde_loc, interp_pphi_min, interp_pphi_max);
				VEC_n2_PPhi_tildeFLIP[iv][izeta] = place_in_range( K_phi_local - PPhi_tilde_loc, interp_pphi_min, interp_pphi_max);
				VEC_n2_PT[iv][izeta] = PT_loc;
				check_for_NaNs("PT_loc", PT_loc, *global_out_stream_ptr);
				check_for_NaNs("PPhi_tilde_loc", PPhi_tilde_loc, *global_out_stream_ptr);
				///////////////////////////////////////////////////////////////
				// Finally, set P^{+/-}
				VEC_n2_Pp[iv][izeta][0] = MT_loc * cosh(P_Y_loc);
				VEC_n2_Pp[iv][izeta][1] = PT_loc * cos(K_phi_local + PPhi_tilde_loc);
				VEC_n2_Pp[iv][izeta][2] = PT_loc * sin(K_phi_local + PPhi_tilde_loc);
				VEC_n2_Pp[iv][izeta][3] = MT_loc * sinh(P_Y_loc);
				VEC_n2_Pm[iv][izeta][0] = VEC_n2_Pp[iv][izeta][0];
				VEC_n2_Pm[iv][izeta][1] = PT_loc * cos(K_phi_local - PPhi_tilde_loc);
				VEC_n2_Pm[iv][izeta][2] = PT_loc * sin(K_phi_local - PPhi_tilde_loc);
				VEC_n2_Pm[iv][izeta][3] = VEC_n2_Pp[iv][izeta][3];
				for (int ii=0; ii<4; ++ii)
				{
					VEC_n2_alpha[iv][izeta][ii] = one_by_Gamma_Mres * VEC_n2_Pp[iv][izeta][ii];
					VEC_n2_alpha_m[iv][izeta][ii] = one_by_Gamma_Mres * VEC_n2_Pm[iv][izeta][ii];
					check_for_NaNs("VEC_n2_alpha[iv][izeta][ii]", VEC_n2_alpha[iv][izeta][ii], *global_out_stream_ptr);
				}
			}
		}
	}
	else
	{
		mT = sqrt(mass*mass + pT*pT);
		double s_min_temp = (m2 + m3)*(m2 + m3);
		double s_max_temp = (Mres - mass)*(Mres - mass);
		gauss_quadrature(n_s_pts, 1, 0.0, 0.0, s_min_temp, s_max_temp, s_pts, s_wts);
		Qfunc = get_Q();
		for (int is = 0; is < n_s_pts; ++is)
		{
			double s_loc = s_pts[is];
			double g_s_loc = g(s_loc);
			VEC_g_s[is] = g_s_loc;
			VEC_s_factor[is] = s_wts[is]*g_s_loc;
			double pstar_loc = sqrt(((Mres+mass)*(Mres+mass) - s_loc)*((Mres-mass)*(Mres-mass) - s_loc))/(2.0*Mres);
			VEC_pstar[is] = pstar_loc;
			double Estar_loc = sqrt(mass*mass + pstar_loc*pstar_loc);
			VEC_Estar[is] = Estar_loc;
			double psBmT = pstar_loc / mT;
			double DeltaY_loc = log(psBmT + sqrt(1.+psBmT*psBmT));
			VEC_DeltaY[is] = DeltaY_loc;
			p_y = 0.0;
			VEC_Yp[is] = p_y + DeltaY_loc;
			VEC_Ym[is] = p_y - DeltaY_loc;
			for(int iv = 0; iv < n_v_pts; ++iv)
			{
				double v_loc = v_pts[iv];
				double P_Y_loc = p_y + v_loc*DeltaY_loc;
				VEC_P_Y[is][iv] = P_Y_loc;
				double mT_ch_P_Y_p_y = mT*cosh(v_loc*DeltaY_loc);
				double x2 = mT_ch_P_Y_p_y*mT_ch_P_Y_p_y - pT*pT;
				VEC_v_factor[is][iv] = v_wts[iv]*DeltaY_loc/sqrt(x2);
				double MTbar_loc = Estar_loc*Mres*mT_ch_P_Y_p_y/x2;
				VEC_MTbar[is][iv] = MTbar_loc;
				double DeltaMT_loc = Mres*pT*sqrt(Estar_loc*Estar_loc - x2)/x2;
				VEC_DeltaMT[is][iv] = DeltaMT_loc;
				VEC_MTp[is][iv] = MTbar_loc + DeltaMT_loc;
				VEC_MTm[is][iv] = MTbar_loc - DeltaMT_loc;
				for(int izeta = 0; izeta < n_zeta_pts; ++izeta)
				{
					double zeta_loc = zeta_pts[izeta];
					double MT_loc = MTbar_loc + cos(zeta_loc)*DeltaMT_loc;
					VEC_MT[is][iv][izeta] = MT_loc;
					VEC_zeta_factor[is][iv][izeta] = zeta_wts[izeta]*MT_loc;
					double PT_loc = sqrt(MT_loc*MT_loc - Mres*Mres);
					///////////////////////////////////////////////////////////////
					// Get \tilde{Phi}...
					double temp_cos_PPhi_tilde_loc = (mT*MT_loc*cosh(P_Y_loc-p_y) - Estar_loc*Mres)/(pT*PT_loc);
					//assume that PPhi_tilde is +ve in next step...
					double temp_sin_PPhi_tilde_loc = sqrt(1. - temp_cos_PPhi_tilde_loc*temp_cos_PPhi_tilde_loc);
					double PPhi_tilde_loc = place_in_range( atan2(temp_sin_PPhi_tilde_loc, temp_cos_PPhi_tilde_loc), interp_pphi_min, interp_pphi_max);
					VEC_PPhi_tilde[is][iv][izeta] = place_in_range( K_phi_local + PPhi_tilde_loc, interp_pphi_min, interp_pphi_max);
					VEC_PPhi_tildeFLIP[is][iv][izeta] = place_in_range( K_phi_local - PPhi_tilde_loc, interp_pphi_min, interp_pphi_max);
					///////////////////////////////////////////////////////////////
					// Finally, set P^{+/-}
					VEC_PT[is][iv][izeta] = PT_loc;
					VEC_Pp[is][iv][izeta][0] = MT_loc * cosh(P_Y_loc);
					VEC_Pp[is][iv][izeta][1] = PT_loc * cos(K_phi_local + PPhi_tilde_loc);
					VEC_Pp[is][iv][izeta][2] = PT_loc * sin(K_phi_local + PPhi_tilde_loc);
					VEC_Pp[is][iv][izeta][3] = MT_loc * sinh(P_Y_loc);
					VEC_Pm[is][iv][izeta][0] = VEC_Pp[is][iv][izeta][0];
					VEC_Pm[is][iv][izeta][1] = PT_loc * cos(K_phi_local - PPhi_tilde_loc);
					VEC_Pm[is][iv][izeta][2] = PT_loc * sin(K_phi_local - PPhi_tilde_loc);
					VEC_Pm[is][iv][izeta][3] = VEC_Pp[is][iv][izeta][3];
					for (int ii=0; ii<4; ++ii)
					{
						VEC_alpha[is][iv][izeta][ii] = one_by_Gamma_Mres * VEC_Pp[is][iv][izeta][ii];
						VEC_alpha_m[is][iv][izeta][ii] = one_by_Gamma_Mres * VEC_Pm[is][iv][izeta][ii];
						check_for_NaNs("VEC_alpha[iv][izeta][ii]", VEC_alpha[is][iv][izeta][ii], *global_out_stream_ptr);
					}
				}
			}
		}
	}

	return;
}


void SourceVariances::Cal_dN_dypTdpTdphi(double** SP_p0, double** SP_px, double** SP_py, double** SP_pz)
{
	double sign = particle_sign;
	double degen = particle_gspin;
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double localmass = particle_mass;

	for(int isurf=0; isurf<FO_length ; ++isurf)
	{
		FO_surf* surf = &FOsurf_ptr[isurf];
		double tau = surf->tau;
		double mu = surf->particle_mu[particle_id];

		double vx = surf->vx;
		double vy = surf->vy;
		double Tdec = surf->Tdec;
		double one_by_Tdec = 1./Tdec;
		double Pdec = surf->Pdec;
		double Edec = surf->Edec;
		double da0 = surf->da0;
		double da1 = surf->da1;
		double da2 = surf->da2;
		double pi00 = surf->pi00;
		double pi01 = surf->pi01;
		double pi02 = surf->pi02;
		double pi11 = surf->pi11;
		double pi12 = surf->pi12;
		double pi22 = surf->pi22;
		double pi33 = surf->pi33;
		double vT = sqrt(vx*vx + vy*vy);
		double gammaT = 1./sqrt(1. - vT*vT);
		double temp_r = surf->r;
		double temp_phi = surf->phi;

		double deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));
      
		for(int ipt = 0; ipt < n_SP_pT; ++ipt)
		{
			double pT = SP_pT[ipt];
			for(int iphi = 0; iphi < n_SP_pphi; ++iphi)
			{
				double px = SP_px[ipt][iphi];
				double py = SP_py[ipt][iphi];
				double cos_phi_m_pphi = cos(temp_phi - SP_pphi[iphi]);
				for(int ieta=0; ieta < eta_s_npts; ++ieta)
				{
					double p0 = SP_p0[ipt][ieta];
					double pz = SP_pz[ipt][ieta];

					//now get distribution function, emission function, etc.
					double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions

					//p^mu d^3sigma_mu: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
					double pdsigma = p0*da0 + px*da1 + py*da2;

					//viscous corrections
					double Wfactor = p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33;
					double deltaf = 0.;
					if (INCLUDE_DELTA_F)
						deltaf = (1. - sign*f0)*Wfactor*deltaf_prefactor;

					//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
					double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

					//ignore points where delta f is large or emission function goes negative from pdsigma
					if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
						S_p = 0.0e0;


					//double S_p = prefactor*pdsigma*f0*(1.+deltaf);
					double symmetry_factor = 2.0;	//eta_s-symmetry
					double S_p_withweight = S_p*tau*eta_s_weight[ieta]*symmetry_factor; //symmetry_factor accounts for the assumed reflection symmetry along eta direction
					dN_dypTdpTdphi[ipt][iphi] += S_p_withweight;
				}	//end of ieta loop
			}		//end of iphi loop
		}			//end of ipt loop
	}				//end of isurf loop
	return;
}

void SourceVariances::Set_source_variances_grid()
{
	for(int ipt = 0; ipt < n_interp_pT_pts; ipt++)
	for(int ipphi = 0; ipphi < n_interp_pphi_pts; ipphi++)
	{
		S_func[ipt][ipphi] = dN_dypTdpTdphi_moments[target_particle_id][0][ipt][ipphi];
		//if (INCLUDE_SOURCE_VARIANCES)
		//{
			x_S[ipt][ipphi] = dN_dypTdpTdphi_moments[target_particle_id][1][ipt][ipphi];
			x2_S[ipt][ipphi] = dN_dypTdpTdphi_moments[target_particle_id][2][ipt][ipphi];
			y_S[ipt][ipphi] = dN_dypTdpTdphi_moments[target_particle_id][3][ipt][ipphi];
			y2_S[ipt][ipphi] = dN_dypTdpTdphi_moments[target_particle_id][4][ipt][ipphi];
			z_S[ipt][ipphi] = dN_dypTdpTdphi_moments[target_particle_id][5][ipt][ipphi];
			z2_S[ipt][ipphi] = dN_dypTdpTdphi_moments[target_particle_id][6][ipt][ipphi];
			t_S[ipt][ipphi] = dN_dypTdpTdphi_moments[target_particle_id][7][ipt][ipphi];
			t2_S[ipt][ipphi] = dN_dypTdpTdphi_moments[target_particle_id][8][ipt][ipphi];
			xy_S[ipt][ipphi] = dN_dypTdpTdphi_moments[target_particle_id][9][ipt][ipphi];
			xz_S[ipt][ipphi] = dN_dypTdpTdphi_moments[target_particle_id][10][ipt][ipphi];
			xt_S[ipt][ipphi] = dN_dypTdpTdphi_moments[target_particle_id][11][ipt][ipphi];
			yz_S[ipt][ipphi] = dN_dypTdpTdphi_moments[target_particle_id][12][ipt][ipphi];
			yt_S[ipt][ipphi] = dN_dypTdpTdphi_moments[target_particle_id][13][ipt][ipphi];
			zt_S[ipt][ipphi] = dN_dypTdpTdphi_moments[target_particle_id][14][ipt][ipphi];

			double cp = cos(SPinterp_pphi[ipphi]), sp = sin(SPinterp_pphi[ipphi]);
			xs_S[ipt][ipphi] = -sp*x_S[ipt][ipphi] + cp*y_S[ipt][ipphi];
			xs2_S[ipt][ipphi] = sp*sp*x2_S[ipt][ipphi] - 2.*sp*cp*xy_S[ipt][ipphi] + cp*cp*y2_S[ipt][ipphi];
			xo_S[ipt][ipphi] = cp*x_S[ipt][ipphi] + sp*y_S[ipt][ipphi];
			xo2_S[ipt][ipphi] = cp*cp*x2_S[ipt][ipphi] + 2.*sp*cp*xy_S[ipt][ipphi] + sp*sp*y2_S[ipt][ipphi];
			xl_S[ipt][ipphi] = z_S[ipt][ipphi];
			xl2_S[ipt][ipphi] = z2_S[ipt][ipphi];
			xo_xs_S[ipt][ipphi] = xy_S[ipt][ipphi]*(cp*cp-sp*sp) + (y2_S[ipt][ipphi] - x2_S[ipt][ipphi])*sp*cp;
			xl_xs_S[ipt][ipphi] = -xz_S[ipt][ipphi] * sp + yz_S[ipt][ipphi] * cp;
			xs_t_S[ipt][ipphi] = -xt_S[ipt][ipphi] * sp + yt_S[ipt][ipphi] * cp;
			xo_xl_S[ipt][ipphi] = xz_S[ipt][ipphi] * cp + yz_S[ipt][ipphi] * sp;
			xo_t_S[ipt][ipphi] = xt_S[ipt][ipphi] * cp + yt_S[ipt][ipphi] * sp;
			xl_t_S[ipt][ipphi] = zt_S[ipt][ipphi];
		//}
	}
	return;
}

void SourceVariances::Interpolate_HBT_radii(int iKT, int iKphi)
{
	double phi_K = K_phi[iKphi];
	double KT = K_T[iKT];

	//int interpolation_mode = 0;		//linear
	int interpolation_mode = 1;		//cubic

	R2_side[iKT][iKphi] = interpolate2D(SPinterp_pT, SPinterp_pphi, R2_side_grid, KT, phi_K, n_interp_pT_pts, n_interp_pphi_pts, interpolation_mode, false, false);
	R2_out[iKT][iKphi] = interpolate2D(SPinterp_pT, SPinterp_pphi, R2_out_grid, KT, phi_K, n_interp_pT_pts, n_interp_pphi_pts, interpolation_mode, false, false);
	R2_long[iKT][iKphi] = interpolate2D(SPinterp_pT, SPinterp_pphi, R2_long_grid, KT, phi_K, n_interp_pT_pts, n_interp_pphi_pts, interpolation_mode, false, false);
	R2_outside[iKT][iKphi] = interpolate2D(SPinterp_pT, SPinterp_pphi, R2_outside_grid, KT, phi_K, n_interp_pT_pts, n_interp_pphi_pts, interpolation_mode, false, false);
	R2_sidelong[iKT][iKphi] = interpolate2D(SPinterp_pT, SPinterp_pphi, R2_sidelong_grid, KT, phi_K, n_interp_pT_pts, n_interp_pphi_pts, interpolation_mode, false, false);
	R2_outlong[iKT][iKphi] = interpolate2D(SPinterp_pT, SPinterp_pphi, R2_outlong_grid, KT, phi_K, n_interp_pT_pts, n_interp_pphi_pts, interpolation_mode, false, false);

	return;
}

void SourceVariances::Calculate_R2_side(int ipt, int ipphi)
{
   double norm = S_func[ipt][ipphi];
   double term1 = xs2_S[ipt][ipphi];
   double term2 = xs_S[ipt][ipphi];

   R2_side_grid[ipt][ipphi] = term1/norm - term2*term2/(norm*norm);

   return;
}

void SourceVariances::Calculate_R2_out(int ipt, int ipphi)
{
   double norm = S_func[ipt][ipphi];
   double term1 = xo2_S[ipt][ipphi] - 2.*beta_perp*xo_t_S[ipt][ipphi] + beta_perp*beta_perp*t2_S[ipt][ipphi];
   double term2 = xo_S[ipt][ipphi] - beta_perp*t_S[ipt][ipphi];

   R2_out_grid[ipt][ipphi] = term1/norm - term2*term2/(norm*norm);

   return;
}

void SourceVariances::Calculate_R2_outside(int ipt, int ipphi)
{
   double norm = S_func[ipt][ipphi];
   double term1 = xo_xs_S[ipt][ipphi] - beta_perp*xs_t_S[ipt][ipphi];
   double term2 = xo_S[ipt][ipphi] - beta_perp*t_S[ipt][ipphi];
   double term3 = xs_S[ipt][ipphi];

   R2_outside_grid[ipt][ipphi] = term1/norm - term2*term3/(norm*norm);

   return;
}

void SourceVariances::Calculate_R2_long(int ipt, int ipphi)
{
   double norm = S_func[ipt][ipphi];
   double term1 = xl2_S[ipt][ipphi] - 2.*beta_l*xl_t_S[ipt][ipphi] + beta_l*beta_l*t2_S[ipt][ipphi];
   double term2 = xl_S[ipt][ipphi] - beta_l*t_S[ipt][ipphi];

   R2_long_grid[ipt][ipphi] = term1/norm - term2*term2/(norm*norm);

   return;
}

void SourceVariances::Calculate_R2_outlong(int ipt, int ipphi)
{
   double norm = S_func[ipt][ipphi];
   double term1 = xo_xl_S[ipt][ipphi] - beta_perp*xl_t_S[ipt][ipphi] - beta_l*xo_t_S[ipt][ipphi] + beta_perp*beta_l*t2_S[ipt][ipphi];
   double term2 = xo_S[ipt][ipphi] - beta_perp*t_S[ipt][ipphi];
   double term3 = xl_S[ipt][ipphi] - beta_l*t_S[ipt][ipphi];

   R2_outlong_grid[ipt][ipphi] = term1/norm - term2*term3/(norm*norm);

   return;
}

void SourceVariances::Calculate_R2_sidelong(int ipt, int ipphi)
{
   double norm = S_func[ipt][ipphi];
   double term1 = xl_xs_S[ipt][ipphi] - beta_l*xs_t_S[ipt][ipphi];
   double term2 = xs_S[ipt][ipphi];
   double term3 = xl_S[ipt][ipphi] - beta_l*t_S[ipt][ipphi];

   R2_sidelong_grid[ipt][ipphi] = term1/norm - term2*term3/(norm*norm);

   return;
}

void SourceVariances::R2_Fourier_transform(int iKT, double plane_psi)
{
	for(int Morder=0; Morder<n_order; ++Morder)
	{
		double cos_mK_phi[n_localp_phi], sin_mK_phi[n_localp_phi];
		for(int i=0; i<n_localp_phi; ++i)
		{
			cos_mK_phi[i] = cos(Morder*(K_phi[i] - plane_psi));
			sin_mK_phi[i] = sin(Morder*(K_phi[i] - plane_psi));
		}
		double temp_sum_side_cos = 0.0e0;
		double temp_sum_side_sin = 0.0e0;
		double temp_sum_out_cos = 0.0e0;
		double temp_sum_out_sin = 0.0e0;
		double temp_sum_outside_cos = 0.0e0;
		double temp_sum_outside_sin = 0.0e0;
		double temp_sum_long_cos = 0.0e0;
		double temp_sum_long_sin = 0.0e0;
		double temp_sum_sidelong_cos = 0.0e0;
		double temp_sum_sidelong_sin = 0.0e0;
		double temp_sum_outlong_cos = 0.0e0;
		double temp_sum_outlong_sin = 0.0e0;
		for(int i=0; i<n_localp_phi; ++i)
		{
			temp_sum_side_cos += R2_side[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
			temp_sum_side_sin += R2_side[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
			temp_sum_out_cos += R2_out[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
			temp_sum_out_sin += R2_out[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
			temp_sum_outside_cos += R2_outside[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
			temp_sum_outside_sin += R2_outside[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
			temp_sum_long_cos += R2_long[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
			temp_sum_long_sin += R2_long[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
			temp_sum_sidelong_cos += R2_sidelong[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
			temp_sum_sidelong_sin += R2_sidelong[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
			temp_sum_outlong_cos += R2_outlong[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
			temp_sum_outlong_sin += R2_outlong[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
		}
		R2_side_C[iKT][Morder] = temp_sum_side_cos/(2.*M_PI);
		R2_side_S[iKT][Morder] = temp_sum_side_sin/(2.*M_PI);
		R2_out_C[iKT][Morder] = temp_sum_out_cos/(2.*M_PI);
		R2_out_S[iKT][Morder] = temp_sum_out_sin/(2.*M_PI);
		R2_outside_C[iKT][Morder] = temp_sum_outside_cos/(2.*M_PI);
		R2_outside_S[iKT][Morder] = temp_sum_outside_sin/(2.*M_PI);
		R2_long_C[iKT][Morder] = temp_sum_long_cos/(2.*M_PI);
		R2_long_S[iKT][Morder] = temp_sum_long_sin/(2.*M_PI);
		R2_sidelong_C[iKT][Morder] = temp_sum_sidelong_cos/(2.*M_PI);
		R2_sidelong_S[iKT][Morder] = temp_sum_sidelong_sin/(2.*M_PI);
		R2_outlong_C[iKT][Morder] = temp_sum_outlong_cos/(2.*M_PI);
		R2_outlong_S[iKT][Morder] = temp_sum_outlong_sin/(2.*M_PI);
	}
	return;
}

//***************************************************************************************************

/*void SourceVariances::test_function(int local_pid)
{
	ostringstream filename_stream_icpl;
	filename_stream_icpl << global_path << "/interpolation_comparison_monval_" << all_particles[local_pid].monval << ".dat";
	ofstream output_icpl(filename_stream_icpl.str().c_str());

	current_level_of_output = 0;

	res_sign_info = sign_of_dN_dypTdpTdphi_moments[local_pid];
	res_log_info = ln_dN_dypTdpTdphi_moments[local_pid];
	res_moments_info = dN_dypTdpTdphi_moments[local_pid];

	double local_pT_min = 0.0, local_pT_max = 10.0, local_pphi_min = 0.0, local_pphi_max = 2.*M_PI;
	double npt = 101., npphi = 1.;
	double Del_pT = (local_pT_max - local_pT_min) / (npt - 1.);
	double Del_pphi = (local_pphi_max - local_pphi_min) / npphi;
	double * result2 = new double [1];
	for (int ipt = 0; ipt < npt; ipt++)
	for (int ipphi = 0; ipphi < npphi; ipphi++)
	{
		double local_pT = local_pT_min + ipt * Del_pT;
		double local_pphi = local_pphi_min + ipphi * Del_pphi;
		double result1 = Cal_dN_dypTdpTdphi_function(local_pid, local_pT, local_pphi);
		result2[0] = 0.0;
if (ipt==(int)npt-1 && ipphi==(int)npphi-1 && local_pid == 1)
	current_level_of_output = 1;
else
	current_level_of_output = 0;
		double result3 = Edndp3_original(local_pT, local_pphi, local_pid, 0);
		Edndp3(local_pT, local_pphi, result2);
		output_icpl << local_pid << "   " << local_pT << "   " << local_pphi << "   " << result1 << "   " << result2[0] << "   " << result3 << endl;
		result2[0] = 0.0;
	}

	delete [] result2;

	output_icpl.close();

	return;
}*/

//End of file
