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
#include<limits>

#include "cfwr.h"
#include "cfwr_lib.h"
#include "Arsenal.h"
#include "Stopwatch.h"
#include "gauss_quadrature.h"
#include "bessel.h"

using namespace std;

const std::complex<double> i(0, 1);
gsl_cheb_series *cs_accel_expK0re, *cs_accel_expK0im, *cs_accel_expK1re, *cs_accel_expK1im;
//bool cheb_set = false;
bool print_stuff = false;

// only need to calculated interpolation grid of spacetime moments for each resonance, NOT each decay channel!
bool recycle_previous_moments = false;
bool recycle_similar_moments = false;
int reso_particle_id_of_moments_to_recycle = -1;
int reso_idx_of_moments_to_recycle = -1;
string reso_name_of_moments_to_recycle = "NULL";
string current_decay_channel_string = "NULL";

//define some parameters for the exact emission function
const double Rad = 5.0, Del_tau = 1.0, tau0 = 5.0, etaf = 0.6;

inline void I(double alpha, double beta, double gamma, complex<double> & I0, complex<double> & I1, complex<double> & I2, complex<double> & I3)
{
	complex<double> ci0, ci1, ck0, ck1, ci0p, ci1p, ck0p, ck1p;
	complex<double> z0 = alpha - i*beta;
	complex<double> z0sq = pow(z0, 2.0);
	double gsq = gamma*gamma;
	complex<double> z = sqrt(z0sq + gsq);
	int errorCode = bessf::cbessik01(z, ci0, ci1, ck0, ck1, ci0p, ci1p, ck0p, ck1p);
	
	I0 = 2.0*ck0;
	I1 = 2.0*z0*ck1 / z;
	I2 = 2.0*z0sq*ck0 / (z*z)
			+ 2.0*(z0sq - gsq)*ck1 / pow(z, 3.0);
	I3 = 2.0*z0*( ( pow(z0, 4.0) - 2.0* z0sq*gsq - 3.0 * pow(gamma, 4.0) ) * ck0 / z
						+ (-6.0*gsq + z0sq*(2.0 + z0sq + gsq)) * ck1
				) / pow(z,5.0);

	return;
}

inline void Ifunc2(double alpha, vector<complex<double> > * I0, vector<complex<double> > * I1, vector<complex<double> > * I2, vector<complex<double> > * I3, int max_n_terms_to_compute)
{
	double beta = 0.0, gamma = 0.0;
	double gsq = gamma*gamma;
	for (int k = 1; k <= max_n_terms_to_compute; ++k)
	{
		complex<double> ci0, ci1, ck0, ck1, ci0p, ci1p, ck0p, ck1p;
		complex<double> z0 = k*alpha - i*beta;
		complex<double> z0sq = z0*z0;
		complex<double> zsq = z0sq + gsq;
		complex<double> z = sqrt(zsq);
		complex<double> zcu = zsq*z;
		complex<double> zqi = zsq*zcu;

		int errorCode = bessf::cbessik01(z, ci0, ci1, ck0, ck1, ci0p, ci1p, ck0p, ck1p);
	
		(*I0).push_back(2.0*ck0);
		(*I1).push_back(2.0*z0*ck1 / z);
		(*I2).push_back(
				2.0*z0sq*ck0 / zsq
				+ 2.0*(z0sq - gsq)*ck1 / zcu
				);
		(*I3).push_back(
				2.0*z0*( ( z0sq*z0sq - 2.0* z0sq*gsq - 3.0 * gsq*gsq ) * ck0 / z
				+ (-6.0*gsq + z0sq*(2.0 + z0sq + gsq)) * ck1 ) / zqi
				);
	}

	return;
}

inline void Iint2(double alpha, double beta, double gamma, double & I0r, double & I1r, double & I2r, double & I3r, double & I0i, double & I1i, double & I2i, double & I3i)
{
	complex<double> z0 = alpha - i*beta;
	complex<double> z0sq = z0*z0;
	double gsq = gamma*gamma;
	complex<double> zsq = z0sq + gsq;
	complex<double> z = sqrt(zsq);
	complex<double> zcu = zsq*z;
	complex<double> zqi = zsq*zcu;
	double ea = exp(-alpha);

//	complex<double> Cci0, Cci1, Cck0, Cck1, Cci0p, Cci1p, Cck0p, Cck1p;
//	int errorCode = bessf::cbessik01(z, Cci0, Cci1, Cck0, Cck1, Cci0p, Cci1p, Cck0p, Cck1p);

	complex<double> ck0(	ea * gsl_cheb_eval (cs_accel_expK0re, alpha),
							ea * gsl_cheb_eval (cs_accel_expK0im, alpha) );
	complex<double> ck1(	ea * gsl_cheb_eval (cs_accel_expK1re, alpha),
							ea * gsl_cheb_eval (cs_accel_expK1im, alpha) );

//cout << "Bessel Check: " << setw(10) << setprecision(8) << ea * gsl_cheb_eval (cs_accel_expK0re, alpha) << "   " << ea * gsl_cheb_eval (cs_accel_expK0im, alpha) << "   " << ea * gsl_cheb_eval (cs_accel_expK1re, alpha) << "   " << ea * gsl_cheb_eval (cs_accel_expK1im, alpha) << "   " << Cck0.real() << "   " << Cck0.imag() << "   " << Cck1.real() << "   " << Cck1.imag() << "   ";
if (print_stuff) cout << "Bessel Check: " << setw(18) << setprecision(16) << alpha << "   " << beta << "   " << gamma << setw(10) << setprecision(8) << endl;
//if (1) exit(8);

//I think this is the fix...
//ck0 = Cck0;
//ck1 = Cck1;

	complex<double> I0 = 2.0*ck0;
	complex<double> I1 = 2.0*z0*ck1 / z;
	complex<double> I2 = 2.0*z0sq*ck0 / zsq
			+ 2.0*(z0sq - gsq)*ck1 / zcu;
	complex<double> I3 = 2.0*z0*( ( z0sq*z0sq - 2.0* z0sq*gsq - 3.0 * gsq*gsq ) * ck0 / z
						+ (-6.0*gsq + z0sq*(2.0 + z0sq + gsq)) * ck1
				) / zqi;

	I0r = I0.real();
	I1r = I1.real();
	I2r = I2.real();
	I3r = I3.real();
	I0i = I0.imag();
	I1i = I1.imag();
	I2i = I2.imag();
	I3i = I3.imag();

	return;
}

inline double Hfactor(double r, double tau)
{
	return (
			exp( -r*r/(2.0*Rad*Rad) - (tau-tau0)*(tau-tau0)/(2.0*Del_tau*Del_tau) ) / (M_PI*Del_tau)
			);
}

inline double eta_t(double r)
{
	return ( etaf*r/Rad );
}

double CorrelationFunction::place_in_range(double phi, double min, double max)
{
	while (phi < min || phi > max)
	{
		if (phi < min) phi += twopi;
		else phi -= twopi;
	}

	return (phi);
}

// ************************************************************
// Compute correlation function at all specified q points for all resonances here
// ************************************************************
void CorrelationFunction::Fourier_transform_emission_function(int iqt, int iqz)
{
	Stopwatch BIGsw;
	global_plane_psi = 0.0;	//for now

	int decay_channel_loop_cutoff = n_decay_channels;			//loop over direct pions and decay_channels
	current_iqt = iqt;
	current_iqz = iqz;
	double loc_qz = qz_pts[iqz];
	double loc_qt = qt_pts[iqt];
	current_pY_shift = 0.5 * log(abs((loc_qt+loc_qz + 1.e-100)/(loc_qt-loc_qz + 1.e-100)));
	//current_pY_shift = 0.0;

	///////
	//on first loop ONLY, initialize necessary HDF files
	if (iqt == 0 && iqz == 0)
	{
		*global_out_stream_ptr << "Initializing HDF files...";
		int HDFInitializationSuccess = Administrate_resonance_HDF_array(0);
		HDFInitializationSuccess = Administrate_target_thermal_HDF_array(0);
		if (!thermal_pions_only && ( !USE_EXACT || USE_CF ) )
			Set_all_Bessel_grids(iqt, iqz);
		//Set_all_Bessel_grids(iqt, iqz, 1);
		*global_out_stream_ptr << "done." << endl << endl;
	}
	else
	{
		*global_out_stream_ptr << "Initializing/opening HDF files...";
		int HDFInitializationSuccess = Administrate_resonance_HDF_array(1);		//open
		HDFInitializationSuccess = Administrate_target_thermal_HDF_array(1);	//open
		if (!thermal_pions_only && ( !USE_EXACT || USE_CF ) )
			Set_all_Bessel_grids(iqt, iqz);
		//Set_all_Bessel_grids(iqt, iqz, 1);
		*global_out_stream_ptr << "done." << endl << endl;
	}
	///////

	*global_out_stream_ptr << "Setting spacetime moments grid..." << endl;
	BIGsw.Start();

	// loop over decay_channels (idc == 0 corresponds to thermal pions)
	for (int idc = 0; idc <= decay_channel_loop_cutoff; ++idc)
	{
		// check whether to do this decay channel
		if (idc > 0 && thermal_pions_only)
			break;
		else if (!Do_this_decay_channel(idc))
			continue;
	
		Set_current_particle_info(idc);

		//if (current_resonance_particle_id != 49)
		//	continue;

		Get_spacetime_moments(idc, iqt, iqz);
	}

	BIGsw.Stop();
	*global_out_stream_ptr << "\t ...finished all (thermal) space-time moments for loop (iqt = " << iqt << ", iqz = " << iqz << ") in " << BIGsw.printTime() << " seconds." << endl;
	
	//only need to calculate spectra, etc. once
	// Now dump all thermal spectra before continuing with resonance decay calculations
	if ( iqt == (qtnpts - 1)/2 && iqz == (qznpts - 1)/2 )
	{
		Dump_spectra_array("thermal_spectra.dat", thermal_spectra);
		//set logs and signs!
		for (int ipid = 0; ipid < Nparticle; ++ipid)
			Set_spectra_logs_and_signs(ipid);
	}

	///////
	*global_out_stream_ptr << "Cleaning up HDF files...";
	{
		int HDFInitializationSuccess = Administrate_resonance_HDF_array(1);		//open
		HDFInitializationSuccess = Administrate_target_thermal_HDF_array(1);	//open

		//get thermal target moments here
		int accessHDFresonanceSpectra = Access_resonance_in_HDF_array(target_particle_id, iqt, iqz, 1, thermal_target_dN_dypTdpTdphi_moments, true);

		/*for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		for (int ipY = 0; ipY < n_pY_pts; ++ipY)
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
			cout << "THERMALDUMP: " << ipT << "   " << ipphi << "   " << ipY << "   " << iqx << "   " << iqy << "   " << current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,0)] << "   " << thermal_target_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,0)] << endl;*/

		if (accessHDFresonanceSpectra < 0)
		{
			cerr << "Failed to get this resonance(local_pid = " << target_particle_id << ") in HDF array!  Exiting..." << endl;
			exit(1);
		}			
		//make sure they are written to separate file here
		accessHDFresonanceSpectra = Access_target_thermal_in_HDF_array(iqt, iqz, 0, thermal_target_dN_dypTdpTdphi_moments);
		if (accessHDFresonanceSpectra < 0)
		{
			cerr << "Failed to set this resonance(local_pid = " << target_particle_id << ") in HDF array!  Exiting..." << endl;
			exit(1);
		}			

		//save thermal moments (without resonance decay feeddown) separately
		int closeHDFresonanceSpectra = Administrate_resonance_HDF_array(2);
		closeHDFresonanceSpectra = Administrate_target_thermal_HDF_array(2);
	}
	*global_out_stream_ptr << "done." << endl << endl;
	///////

   return;
}

void CorrelationFunction::Compute_phase_space_integrals(int iqt, int iqz)
{
	if (thermal_pions_only)
	{
		*global_out_stream_ptr << "Thermal pions only: no phase-space integrals need to be computed." << endl;
		return;
	}

	Stopwatch BIGsw;
	int decay_channel_loop_cutoff = n_decay_channels;			//loop over direct pions and decay_channels

	//set needed q-points
	Set_qlist(iqt, iqz);

	int openHDFresonanceSpectra = Administrate_resonance_HDF_array(1);	//open

	*global_out_stream_ptr << "Computing all phase-space integrals..." << endl;
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

		//if (current_resonance_particle_id != 49)
		//	continue;

		Allocate_decay_channel_info();				// allocate needed memory

		// ************************************************************
		// begin resonance decay calculations here...
		// ************************************************************

		Load_resonance_and_daughter_spectra(decay_channels[idc-1].resonance_particle_id, iqt, iqz);

		Stopwatch decay_sw;
		decay_sw.Start();
		for (int idc_DI = 0; idc_DI < current_reso_nbody; ++idc_DI)
		{
			int daughter_resonance_particle_id = -1;
			if (!Do_this_daughter_particle(idc, idc_DI, &daughter_resonance_particle_id))
				continue;

			Set_current_daughter_info(idc, idc_DI);

			Do_resonance_integrals(current_resonance_particle_id, daughter_resonance_particle_id, idc, iqt, iqz);
		}
		Update_daughter_spectra(decay_channels[idc-1].resonance_particle_id, iqt, iqz);

		Delete_decay_channel_info();				// free up memory
		decay_sw.Stop();
		*global_out_stream_ptr << " - Finished decay loop for " << decay_channels[idc-1].resonance_name << " in " << decay_sw.printTime() << " seconds." << endl;
	}			
								// END of decay channel loop
	BIGsw.Stop();
	*global_out_stream_ptr << "\t ...finished computing all phase-space integrals for loop (iqt = "
							<< iqt << ", iqz = " << iqz << ") in " << BIGsw.printTime() << " seconds." << endl;

	if ( iqt == (qtnpts - 1)/2 && iqz == (qznpts - 1)/2 )
		Dump_spectra_array("full_spectra.dat", spectra);

	return;
}
//////////////////////////////////////////////////////////////////////////////////
// End of main routines for setting up computation of correlation function

bool CorrelationFunction::Do_this_decay_channel(int dc_idx)
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
	bool tmp_bool = decay_channels[dc_idx-1].include_channel;
	if (!tmp_bool && VERBOSE > 0) *global_out_stream_ptr << endl << local_name << ": skipping decay " << current_decay_channel_string << "." << endl;

	return (tmp_bool);
}

// ************************************************************
// Checks whether to do daughter particle for any given decay channel
// ************************************************************
bool CorrelationFunction::Do_this_daughter_particle(int dc_idx, int daughter_idx, int * daughter_resonance_pid)
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

void CorrelationFunction::Set_current_particle_info(int dc_idx)
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

		if (VERBOSE > 0) *global_out_stream_ptr << endl << local_name << ": doing decay " << current_decay_channel_string << "." << endl
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
				if (VERBOSE > 0) *global_out_stream_ptr << "\t * " << decay_channels[dc_idx-1].resonance_name << " (same as the last one)." << endl;
			}
			else if ( Search_for_similar_particle( temp_reso_idx, &similar_particle_idx ) )
			{
				//previous resonance is NOT the same as this one BUT this one is sufficiently similar to some preceding one...
				recycle_previous_moments = false;
				recycle_similar_moments = true;
				reso_particle_id_of_moments_to_recycle = chosen_resonances[similar_particle_idx];
				reso_idx_of_moments_to_recycle = similar_particle_idx;
				if (VERBOSE > 0) *global_out_stream_ptr << "\t * " << decay_channels[dc_idx-1].resonance_name << " (different from the last one, but close enough to "
														<< all_particles[reso_particle_id_of_moments_to_recycle].name << ")." << endl;
			}
			else
			{
				recycle_previous_moments = false;
				recycle_similar_moments = false;
				reso_particle_id_of_moments_to_recycle = -1;	//guarantees it won't be used spuriously
				reso_idx_of_moments_to_recycle = -1;
				if (VERBOSE > 0) *global_out_stream_ptr << "\t * " << decay_channels[dc_idx-1].resonance_name << " (different from the last one --> calculating afresh)." << endl;
			}
		}
	}
	
	return;
}

void CorrelationFunction::Set_current_daughter_info(int dc_idx, int daughter_idx)
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
			break;
		default:
			cerr << "Set_current_daughter_info(): shouldn't have ended up here, bad value of current_reso_nbody!" << endl;
			exit(1);
	}
}

bool CorrelationFunction::Search_for_similar_particle(int reso_idx, int * result)
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
bool CorrelationFunction::particles_are_the_same(int reso_idx1, int reso_idx2)
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

void CorrelationFunction::Recycle_spacetime_moments()
{
	//*global_out_stream_ptr << "PIDs: " << current_resonance_particle_id << "   " << reso_particle_id_of_moments_to_recycle << endl;
	int HDFcopyChunkSuccess = Copy_chunk(current_resonance_particle_id, reso_particle_id_of_moments_to_recycle);
	if (HDFcopyChunkSuccess < 0) exit(1);

	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	{
		spectra[current_resonance_particle_id][ipT][ipphi] = spectra[reso_particle_id_of_moments_to_recycle][ipT][ipphi];
		thermal_spectra[current_resonance_particle_id][ipT][ipphi] = thermal_spectra[reso_particle_id_of_moments_to_recycle][ipT][ipphi];
	}

	return;
}

//**************************************************************
//**************************************************************
void CorrelationFunction::Load_resonance_and_daughter_spectra(int local_pid, int iqt, int iqz)
{
	// get parent resonance spectra, set logs and signs arrays that are needed for interpolation
	int getHDFresonanceSpectra = Access_resonance_in_HDF_array(local_pid, iqt, iqz, 1, current_dN_dypTdpTdphi_moments);

	// get spectra for all daughters, set all of the logs and signs arrays that are needed for interpolation
	int n_daughters = Set_daughter_list(local_pid);

	if (n_daughters > 0)
	{
		int d_idx = 0;

		Setup_current_daughters_dN_dypTdpTdphi_moments(n_daughters);

		for (set<int>::iterator it = daughter_resonance_indices.begin(); it != daughter_resonance_indices.end(); ++it)
		{
			int daughter_pid = *it;		//daughter pid is pointed to by iterator
			getHDFresonanceSpectra = Access_resonance_in_HDF_array(daughter_pid, iqt, iqz, 1, current_daughters_dN_dypTdpTdphi_moments[d_idx]);
			++d_idx;
		}
		/*for (int id = 0; id < n_daughters; ++id)
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
                for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
                for (int ipY = 0; ipY < n_pY_pts; ++ipY)
                for (int iqx = 0; iqx < qxnpts; ++iqx)
                for (int iqy = 0; iqy < qynpts; ++iqy)
			cout << "CHECKDUMP: " << id << "   " << ipT << "   " << ipphi << "   " << ipY << "   " << iqx << "   " << iqy << "   " << current_daughters_dN_dypTdpTdphi_moments[id][fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,0)] << "   " << current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,0)] << "   " << thermal_target_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,0)] << endl;*/
	}
	else
	{
		cerr << "Particle is stable, shouldn't have ended up here!  Something went wrong..." << endl;
		exit(1);
	}

	return;
}

void CorrelationFunction::Update_daughter_spectra(int local_pid, int iqt, int iqz)
{
	int d_idx = 0;
	for (set<int>::iterator it = daughter_resonance_indices.begin(); it != daughter_resonance_indices.end(); ++it)
	{
		int daughter_pid = *it;		//daughter pid is pointed to by iterator
		int setHDFresonanceSpectra = Access_resonance_in_HDF_array(daughter_pid, iqt, iqz, 0, current_daughters_dN_dypTdpTdphi_moments[d_idx]);

		++d_idx;
	}

	// cleanup previous iteration and setup the new one
	Cleanup_current_daughters_dN_dypTdpTdphi_moments(daughter_resonance_indices.size());

	return;
}

void CorrelationFunction::Set_spectra_logs_and_signs(int local_pid)
{
	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	{
		log_spectra[local_pid][ipT][ipphi] = log(abs(spectra[local_pid][ipT][ipphi])+1.e-100);
		sign_spectra[local_pid][ipT][ipphi] = sgn(spectra[local_pid][ipT][ipphi]);
	}

	return;
}

//**************************************************************
//**************************************************************

void CorrelationFunction::Get_spacetime_moments(int dc_idx, int iqt, int iqz)
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
	if (recycle_previous_moments && dc_idx > 1)	// same as earlier resonance
	{
		if (VERBOSE > 0) *global_out_stream_ptr << local_name
			<< ": new parent resonance (" << decay_channels[current_decay_channel_idx-1].resonance_name << ", dc_idx = " << current_decay_channel_idx
			<< " of " << n_decay_channels << ") same as preceding parent resonance \n\t\t--> reusing old dN_dypTdpTdphi_moments!" << endl;
	}
	else if (recycle_similar_moments && dc_idx > 1)	// sufficiently similar (but different) earlier resonance
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
			if (!recycle_previous_moments && !recycle_similar_moments) 
			{
				if (VERBOSE > 0)
					*global_out_stream_ptr << local_name
						<< ": new parent resonance (" << decay_channels[current_decay_channel_idx-1].resonance_name << ", dc_idx = " << current_decay_channel_idx
						<< " of " << n_decay_channels << ") dissimilar from all preceding decay_channels \n\t\t--> calculating new dN_dypTdpTdphi_moments!" << endl;
			}
			else
			{
				cerr << "You shouldn't have ended up here!" << endl;
				exit(1);
			}
		}

		//allows to omit thermal spectra calculations from specified resonances, e.g., all resonances which contribute up to 60% of decay pions
		if (find(osr.begin(), osr.end(), current_resonance_particle_id) != osr.end())
			*global_out_stream_ptr << "  --> ACTUALLY SKIPPING WEIGHTED THERMAL SPECTRA FOR " << local_name << endl;
		else
		{
			*global_out_stream_ptr << "  --> ACTUALLY DOING WEIGHTED THERMAL SPECTRA FOR " << local_name << endl;
			Set_dN_dypTdpTdphi_moments(current_resonance_particle_id, iqt, iqz);
		}
	}
//**************************************************************
//Spacetime moments now set
//**************************************************************
	return;
}

void CorrelationFunction::Reset_FOcells_array()
{
	for (int iFOipT = 0; iFOipT < FO_length * n_pT_pts; ++iFOipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		FOcells_to_include[iFOipT][ipphi] = std::numeric_limits<int>::max();
	return;
}

void CorrelationFunction::Set_dN_dypTdpTdphi_moments(int local_pid, int iqt, int iqz)
{
	double localmass = all_particles[local_pid].mass;
	int local_particle_mode = 1;
	int local_na = n_alpha_points_PIONS;
	string local_name = "Thermal pion(+)";
	if (local_pid != target_particle_id)
	{
		local_name = all_particles[local_pid].name;
		local_particle_mode = 0;
		local_na = n_alpha_points;
	}
	
	// get spectra at each fluid cell, sort by importance
	Stopwatch sw, sw_qtqzpY;
	sw.Start();
	//if ( iqt == (qtnpts - 1)/2 && iqz == (qznpts - 1)/2 )
	if ( iqt == 0 && iqz == 0 )	//do it right away, so we know which FOcells to skip hereafter
	{
		*global_out_stream_ptr << "Computing un-weighted thermal spectra..." << endl;
		if ( USE_EXACT && !USE_CF )
		{
			Cal_dN_dypTdpTdphi_no_weights_toy(local_pid);
		}
		else
		{
			if (local_pid == target_particle_id)
				Cal_dN_dypTdpTdphi_no_weights_Yeq0_alternate();
			else
				Cal_dN_dypTdpTdphi_no_weights(local_pid);
		}
	}
	else						//otherwise, be sure to load the important FOcells from files!
	{
		*global_out_stream_ptr << "Loading important FOcells from file...";
		Load_FOcells(local_pid);
		*global_out_stream_ptr << "done." << endl;
	}

	if (local_particle_mode == 1)	//do pions at the end
		return;

	// get weighted spectra with only most important fluid cells, up to given threshhold
	*global_out_stream_ptr << "Computing weighted thermal spectra..." << endl;

	if ( USE_EXACT && !USE_CF )
	{
		// Loop over pY points
		for (int ipY = 0; ipY < n_pY_pts; ++ipY)
		{
			sw_qtqzpY.Reset();
			sw_qtqzpY.Start();
			ch_SP_pY[ipY] = cosh(SP_Del_pY[ipY] + current_pY_shift);
			sh_SP_pY[ipY] = sinh(SP_Del_pY[ipY] + current_pY_shift);
//cout << "CHECKloop: " << ipY << "   " << iqt << "   " << iqz << "   " << current_pY_shift << "   "
//		<< qt_pts[iqt] << "   " << qz_pts[iqz] << "   " << SP_Del_pY[ipY] << "   "
//		<< ch_SP_pY[ipY] << "   " << sh_SP_pY[ipY] << endl;

			//load appropriate Bessel coefficients
			Cal_dN_dypTdpTdphi_with_weights_toy(local_pid, iqt, iqz, ipY,
												current_dN_dypTdpTdphi_moments);

			sw_qtqzpY.Stop();
			if (VERBOSE > 1) *global_out_stream_ptr
							<< "Finished loop with ( iqt, iqz, ipY ) = ( "
							<< iqt << ", " << iqz << ", " << ipY << " ) in "
							<< sw_qtqzpY.printTime() << " seconds." << endl;
		}
	}
	else
	{
		//prepare for reading...
		int HDFcode = Administrate_besselcoeffs_HDF_array(1, local_particle_mode);	//open
		double * BC_chunk = new double [4 * FO_length * local_na];
		cs_accel_expK0re = gsl_cheb_alloc (local_na - 1);
		cs_accel_expK0im = gsl_cheb_alloc (local_na - 1);
		cs_accel_expK1re = gsl_cheb_alloc (local_na - 1);
		cs_accel_expK1im = gsl_cheb_alloc (local_na - 1);

		///////////////////////////////////
		// Loop over pY points
		///////////////////////////////////
		for (int ipY = 0; ipY < n_pY_pts; ++ipY)
		{
			sw_qtqzpY.Reset();
			sw_qtqzpY.Start();
			ch_SP_pY[ipY] = cosh(SP_Del_pY[ipY] + current_pY_shift);
			sh_SP_pY[ipY] = sinh(SP_Del_pY[ipY] + current_pY_shift);

			//load appropriate Bessel coefficients
			HDFcode = Access_besselcoeffs_in_HDF_array(ipY, 1, BC_chunk, local_particle_mode);
			Cal_dN_dypTdpTdphi_with_weights(local_pid, ipY, iqt, iqz, BC_chunk, local_particle_mode);

			sw_qtqzpY.Stop();
			if (VERBOSE > 1) *global_out_stream_ptr << "Finished loop with ( iqt, iqz, ipY ) = ( "
							<< iqt << ", " << iqz << ", " << ipY << " ) in "
							<< sw_qtqzpY.printTime() << " seconds." << endl;
		}

		HDFcode = Administrate_besselcoeffs_HDF_array(2, local_particle_mode);
		delete [] BC_chunk;
	}


	//HDF5 block
	if (local_pid == target_particle_id && MIDRAPIDITY_PIONS_ONLY)	//nothing to store
		;
	else
	{
		// store in HDF5 file
		int setHDFresonanceSpectra = Access_resonance_in_HDF_array(local_pid, iqt, iqz, 0, current_dN_dypTdpTdphi_moments);
		if (setHDFresonanceSpectra < 0)
		{
			cerr << "Failed to set this resonance(local_pid = "
					<< local_pid << ") in HDF resonance array!  Exiting..." << endl;
			exit(1);
		}

		if (local_pid == target_particle_id)
		{
			int HDFInitializationSuccess = Administrate_target_thermal_HDF_array(1);	//open
			int setHDFtargetSpectra = Access_target_thermal_in_HDF_array(iqt, iqz, 0, current_dN_dypTdpTdphi_moments, true);
			if (setHDFtargetSpectra < 0)
			{
				cerr << "Failed to set this resonance(local_pid = "
						<< local_pid << ") in HDF thermal target array!  Exiting..." << endl;
				exit(1);
			}
			HDFInitializationSuccess = Administrate_target_thermal_HDF_array(2);	//close
		}
	}

	if ( iqt == 0 && iqz == 0 )	//if the FOcells were computed this loop, be sure to dump them to files!
		Dump_FOcells(local_pid);
	else						//otherwise, just clear everything
		Reset_FOcells_array();

	sw.Stop();
	*global_out_stream_ptr << "Took " << sw.printTime() << " seconds to set dN/dypTdpTdphi moments." << endl;

	return;
}

void CorrelationFunction::Cal_dN_dypTdpTdphi_no_weights(int local_pid)
{
	// set particle information
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double loc_Tdec = (&FOsurf_ptr[0])->Tdec;
	double loc_Pdec = (&FOsurf_ptr[0])->Pdec;
	double loc_Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./loc_Tdec;

	double deltaf_prefactor = 0.;
	if (use_delta_f)
		deltaf_prefactor = 1./(2.0*loc_Tdec*loc_Tdec*(loc_Edec+loc_Pdec));

	double eta_s_symmetry_factor = 2.0;

	Stopwatch sw_ThermalResonanceSpectra;

	//Time full calculation of thermal spectra for this resonance
	sw_ThermalResonanceSpectra.Start();

	//////////////////////////////////////////////////
	// Organize calculation of boost-invariant spectra
	// by looping over pT and pphi
	//////////////////////////////////////////////////
	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	{
		double loc_pT = SP_pT[ipT];
		double loc_pphi = SP_pphi[ipphi];
		int pTpphi_index = ipT * n_pphi_pts + ipphi;
		double sin_pphi = sin_SP_pphi[ipphi];
		double cos_pphi = cos_SP_pphi[ipphi];

		double px = loc_pT*cos_pphi;
		double py = loc_pT*sin_pphi;
		double loc_mT = sqrt(loc_pT*loc_pT+localmass*localmass);

		double spectra_at_pTpphi = 0.0;

		////////////////////////////////////////////////
		// Loop over freeze-out surface fluid cells
		////////////////////////////////////////////////
		for (int isurf = 0; isurf < FO_length; ++isurf)
		{
			//////////////////////////////////
			//////////////////////////////////
			FO_surf * surf = &FOsurf_ptr[isurf];

			double tau = surf->tau;
			double rpt = surf->r;
			double phipt = place_in_range(surf->phi, SP_pphi_min, SP_pphi_max);

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

			double A = tau*prefactor*loc_mT*da0;
			double B = tau*prefactor*(px*da1 + py*da2);
			double C = deltaf_prefactor;

			double a = loc_mT*loc_mT*(pi00 + pi33);
			double b = -2.0*loc_mT*(px*pi01 + py*pi02);
			double c = px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 - loc_mT*loc_mT*pi33;

			double alpha = ( USE_EXACT ) ? one_by_Tdec*loc_mT*cosh(eta_t(rpt)) : one_by_Tdec*gammaT*loc_mT;
			double loc_beta = 0.0, loc_gamma = 0.0;
			double transverse_f0 = exp( one_by_Tdec*(gammaT*(px*vx + py*vy) + mu) );

			complex<double> I0_a_b_g, I1_a_b_g, I2_a_b_g, I3_a_b_g;
			I(alpha, loc_beta, loc_gamma, I0_a_b_g, I1_a_b_g, I2_a_b_g, I3_a_b_g);
			double term1 = A*Hfactor(rpt, tau)*exp( one_by_Tdec*loc_pT*sinh(eta_t(rpt))*cos(phipt - loc_pphi) )*I1_a_b_g.real();
			double term2 = 0.0, term3 = 0.0;
			
			if ( !USE_EXACT )
			{
				complex<double> I0_2a_b_g, I1_2a_b_g, I2_2a_b_g, I3_2a_b_g;
				I(2.0*alpha, loc_beta, loc_gamma, I0_2a_b_g, I1_2a_b_g, I2_2a_b_g, I3_2a_b_g);

				term1 = transverse_f0 * (A*I1_a_b_g.real() + B*I0_a_b_g.real());
				term2 = C * transverse_f0
								* ( A*a*I3_a_b_g.real() + (B*a+b*A)*I2_a_b_g.real() + (B*b+c*A)*I1_a_b_g.real() + B*c*I0_a_b_g.real() );
				term3 = -sign * C * transverse_f0 * transverse_f0
								* ( A*a*I3_2a_b_g.real() + (B*a+b*A)*I2_2a_b_g.real() + (B*b+c*A)*I1_2a_b_g.real() + B*c*I0_2a_b_g.real() );
			}

			double FOcell_density = term1 + term2 + term3;
			//////////////////////////////////
			//////////////////////////////////

			//////////////////////////////////
			// Now decide what to do with this FO cell
			//ignore points where delta f is large or emission function goes negative from pdsigma
			if ( abs(term2 + term3) > abs(term1) )
			{
				FOcell_density = term1;	//if viscous corrections large, just ignore them
				FOcells_to_include[isurf * n_pT_pts + ipT][ipphi] = 1;		//1 means I cut out viscous corrections but kept equilibrium distribution
			}
			else
				FOcells_to_include[isurf * n_pT_pts + ipT][ipphi] = 2;		//2 means I kept everything
			if ( flagneg == 1 && FOcell_density < tol )	//whether or not viscous corrections were large, if S goes negative, set it to zero
			{												// and ignore on subsequent loops
				FOcell_density = 0.0;	//if negative, set to 0
				FOcells_to_include[isurf * n_pT_pts + ipT][ipphi] = 0;		//0 means I cut out everything (i.e., this cell ignored)
				continue;
			}

			// add FOdensity into full spectra at this pT, pphi
			spectra_at_pTpphi += FOcell_density;
		}		//end of isurf loop

		//update spectra
		spectra[local_pid][ipT][ipphi] = spectra_at_pTpphi;
		thermal_spectra[local_pid][ipT][ipphi] = spectra_at_pTpphi;
		log_spectra[local_pid][ipT][ipphi] = log(abs(spectra_at_pTpphi) + 1.e-100);
		sign_spectra[local_pid][ipT][ipphi] = sgn(spectra_at_pTpphi);
	}		// end of pT, pphi loop

	sw_ThermalResonanceSpectra.Stop();
	*global_out_stream_ptr << "\t\t\t*** Took " << sw_ThermalResonanceSpectra.printTime() << " seconds for whole function." << endl;

	return;
}

//////////////////////////////////////////
void CorrelationFunction::Cal_dN_dypTdpTdphi_with_weights(int local_pid, int ipY, int iqt, int iqz, double * BC_chunk, int local_part_mode)
{
	Stopwatch sw, sw_FOsurf;
	sw.Start();

	// set particle information
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;

	int local_na = n_alpha_points;
	//double alpha_min = 4.0, alpha_max = 75.0;
	double alpha_min = 4.0, alpha_max = 200.0;
	if (local_part_mode == 1)
	{
		local_na = n_alpha_points_PIONS;
		alpha_min = 0.5;
		alpha_max = 400.0;
	}

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double eta_s_symmetry_factor = 2.0;
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;
	double deltaf_prefactor = 0.;
	if (use_delta_f)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	double * expK0_Bessel_re = new double [local_na];
	double * expK0_Bessel_im = new double [local_na];
	double * expK1_Bessel_re = new double [local_na];
	double * expK1_Bessel_im = new double [local_na];
	cs_accel_expK0re->a = alpha_min;
	cs_accel_expK0re->b = alpha_max;
	cs_accel_expK0im->a = alpha_min;
	cs_accel_expK0im->b = alpha_max;
	cs_accel_expK1re->a = alpha_min;
	cs_accel_expK1re->b = alpha_max;
	cs_accel_expK1im->a = alpha_min;
	cs_accel_expK1im->b = alpha_max;

	for (int ia = 0; ia < local_na; ++ia)
	{
		expK0_Bessel_re[ia] = 0.0;
		expK0_Bessel_im[ia] = 0.0;
		expK1_Bessel_re[ia] = 0.0;
		expK1_Bessel_im[ia] = 0.0;
	}

	double I0_a_b_g_re, I1_a_b_g_re, I2_a_b_g_re, I3_a_b_g_re;
	double I0_2a_b_g_re, I1_2a_b_g_re, I2_2a_b_g_re, I3_2a_b_g_re;
	double I0_a_b_g_im, I1_a_b_g_im, I2_a_b_g_im, I3_a_b_g_im;
	double I0_2a_b_g_im, I1_2a_b_g_im, I2_2a_b_g_im, I3_2a_b_g_im;

	double C = deltaf_prefactor;

	double ** alt_long_array_CR = new double * [qxnpts * qynpts];
	double ** alt_long_array_CI = new double * [qxnpts * qynpts];
	double ** alt_long_array_SR = new double * [qxnpts * qynpts];
	double ** alt_long_array_SI = new double * [qxnpts * qynpts];
	for (int isa = 0; isa < qxnpts * qynpts; ++isa)
	{
		alt_long_array_CR[isa] = new double [n_pT_pts * n_pphi_pts];
		alt_long_array_CI[isa] = new double [n_pT_pts * n_pphi_pts];
		alt_long_array_SR[isa] = new double [n_pT_pts * n_pphi_pts];
		alt_long_array_SI[isa] = new double [n_pT_pts * n_pphi_pts];
		for (int isa2 = 0; isa2 < n_pT_pts * n_pphi_pts; ++isa2)
		{
			alt_long_array_CR[isa][isa2] = 0.0;
			alt_long_array_CI[isa][isa2] = 0.0;
			alt_long_array_SR[isa][isa2] = 0.0;
			alt_long_array_SI[isa][isa2] = 0.0;
		}
	}

	/////////////////////////////////////////////////////////////
	// Loop over all freeze-out surface fluid cells (for now)
	/////////////////////////////////////////////////////////////
	int iBC = 0;
	for (int isurf = 0; isurf < FO_length; ++isurf)
	{
		FO_surf * surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double rpt = surf->r;
		double phipt = place_in_range(surf->phi, SP_pphi_min, SP_pphi_max);

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

		double qt = qt_pts[iqt];
		double qz = qz_pts[iqz];
		double ch_pY = ch_SP_pY[ipY];
		double sh_pY = sh_SP_pY[ipY];
		double beta = tau * hbarCm1 * ( qt*ch_pY - qz*sh_pY );
		double gamma = tau * hbarCm1 * ( qz*ch_pY - qt*sh_pY );

		// Load Bessel Chebyshev coefficients
		for (int ia = 0; ia < local_na; ++ia)
			expK0_Bessel_re[ia] = BC_chunk[iBC++];
		for (int ia = 0; ia < local_na; ++ia)
			expK0_Bessel_im[ia] = BC_chunk[iBC++];
		for (int ia = 0; ia < local_na; ++ia)
			expK1_Bessel_re[ia] = BC_chunk[iBC++];
		for (int ia = 0; ia < local_na; ++ia)
			expK1_Bessel_im[ia] = BC_chunk[iBC++];

		cs_accel_expK0re->c = expK0_Bessel_re;
		cs_accel_expK0im->c = expK0_Bessel_im;
		cs_accel_expK1re->c = expK1_Bessel_re;
		cs_accel_expK1im->c = expK1_Bessel_im;

		double * tmpX = oscx[isurf];
		double * tmpY = oscy[isurf];
		double short_array_C[n_pT_pts * n_pphi_pts], short_array_S[n_pT_pts * n_pphi_pts];
		for (int isa = 0; isa < n_pT_pts * n_pphi_pts; ++isa)
		{
			short_array_C[isa] = 0.0;
			short_array_S[isa] = 0.0;
		}

		/////////////////////////////////////////////////////
		// Loop over pT and pphi points (as fast as possible)
		/////////////////////////////////////////////////////
		int iidx = 0;
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		{
			vector<int> FOcells_to_do = FOcells_to_include[isurf * n_pT_pts + ipT];

			double pT = SP_pT[ipT];
			double mT = sqrt(pT*pT+localmass*localmass);
			double alpha = ( USE_EXACT ) ? one_by_Tdec*mT*cosh(eta_t(rpt)) : one_by_Tdec*gammaT*mT;
			
//print_stuff = bool( /*ipT == n_pT_pts - 1 &&*/ ipY==ipY0 && iqt==iqt0 && iqz==iqz0 );
			Iint2(alpha, beta, gamma, I0_a_b_g_re, I1_a_b_g_re, I2_a_b_g_re, I3_a_b_g_re, I0_a_b_g_im, I1_a_b_g_im, I2_a_b_g_im, I3_a_b_g_im);
			if ( !USE_EXACT )
				Iint2(2.0*alpha, beta, gamma, I0_2a_b_g_re, I1_2a_b_g_re, I2_2a_b_g_re, I3_2a_b_g_re, I0_2a_b_g_im, I1_2a_b_g_im, I2_2a_b_g_im, I3_2a_b_g_im);

			double A = tau*prefactor*mT*da0;
			double a = mT*mT*(pi00 + pi33);

			for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
			{
				int do_this_FOcell = FOcells_to_do[ipphi];
				// initialize transverse momentum information
				double px = pT*cos_SP_pphi[ipphi];
				double py = pT*sin_SP_pphi[ipphi];
				double pphi = SP_pphi[ipphi];

				double B = tau*prefactor*(px*da1 + py*da2);
				double b = -2.0*mT*(px*pi01 + py*pi02);
				double c = px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 - mT*mT*pi33;

				double overall_S_factor = A*Hfactor(rpt, tau)*exp( one_by_Tdec*pT*sinh(eta_t(rpt))*cos(phipt - pphi) );
				double term1_re = overall_S_factor*I1_a_b_g_re;
				double term1_im = overall_S_factor*I1_a_b_g_im;
				double term2_re = 0.0, term3_re = 0.0, term2_im = 0.0, term3_im = 0.0;

				if ( !USE_EXACT )
				{
					double transverse_f0 = exp( one_by_Tdec*(gammaT*(px*vx + py*vy) + mu) );
					term1_re = transverse_f0 * (A*I1_a_b_g_re + B*I0_a_b_g_re);
					term1_im = transverse_f0 * (A*I1_a_b_g_im + B*I0_a_b_g_im);

					double c1 = A*a, c2 = B*a+b*A, c3 = B*b+c*A, c4 = B*c;
					double C1 = C * transverse_f0;
					double C2 = -sign * transverse_f0 * C1;
					term2_re = C1 * ( c1*I3_a_b_g_re + c2*I2_a_b_g_re + c3*I1_a_b_g_re + c4*I0_a_b_g_re );
					term3_re = C2 * ( c1*I3_2a_b_g_re + c2*I2_2a_b_g_re + c3*I1_2a_b_g_re + c4*I0_2a_b_g_re );
					term2_im = C1 * ( c1*I3_a_b_g_im + c2*I2_a_b_g_im + c3*I1_a_b_g_im + c4*I0_a_b_g_im );
					term3_im = C2 * ( c1*I3_2a_b_g_im + c2*I2_2a_b_g_im + c3*I1_2a_b_g_im + c4*I0_2a_b_g_im );
				}

				short_array_C[iidx] = (do_this_FOcell>0)*term1_re + (do_this_FOcell==2)*(term2_re + term3_re);
				short_array_S[iidx++] = (do_this_FOcell>0)*term1_im + (do_this_FOcell==2)*(term2_im + term3_im);
			}
		}

		////////////////////////////////////////
		// Loop over qx, qy, pT, and pphi points
		////////////////////////////////////////
		long idx = 0;
		const long iidx_end = (long)n_pT_pts * (long)n_pphi_pts;
		for (int iqx = 0; iqx < (qxnpts+1)/2; ++iqx)
		{
			double cosAx = tmpX[iqx * 2 + 0], sinAx = tmpX[iqx * 2 + 1];
			for (int iqy = 0; iqy < qynpts; ++iqy)
			{
				double cosAy = tmpY[iqy * 2 + 0], sinAy = tmpY[iqy * 2 + 1];
				double cos_trans_Fourier = cosAx*cosAy - sinAx*sinAy;
				double sin_trans_Fourier = sinAx*cosAy + cosAx*sinAy;
				double * ala_CR = alt_long_array_CR[idx];
				double * ala_CI = alt_long_array_CI[idx];
				double * ala_SR = alt_long_array_SR[idx];
				double * ala_SI = alt_long_array_SI[idx++];
				long iidx_local = 0;
				while ( iidx_local < iidx_end )
				{
					double cos_qx_S_x_K = short_array_C[iidx_local];
					ala_CR[iidx_local] += cos_trans_Fourier * cos_qx_S_x_K;
					ala_CI[iidx_local++] -= sin_trans_Fourier * cos_qx_S_x_K;	//phi_T comes with extra minus sign
				}
				iidx_local = 0;
				while ( iidx_local < iidx_end )
				{
					double sin_qx_S_x_K = short_array_S[iidx_local];
					//ala_SR[iidx_local] += cos_trans_Fourier * sin_qx_S_x_K;
					//ala_SI[iidx_local++] += sin_trans_Fourier * sin_qx_S_x_K;
					ala_SR[iidx_local] += sin_trans_Fourier * sin_qx_S_x_K;		//phi_T comes with extra minus sign which cancels with minus sign on sin*sin
					ala_SI[iidx_local++] += cos_trans_Fourier * sin_qx_S_x_K;
				}
			}
		}
	}

	//use reflection symmetry in transverse plane to get speed-up
	for (int iqx = 0; iqx < (qxnpts-1)/2; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	{
		alt_long_array_CR[(qxnpts-iqx-1) * qynpts + (qynpts-iqy-1)][ipT * n_pphi_pts + ipphi] = alt_long_array_CR[iqx * qynpts + iqy][ipT * n_pphi_pts + ipphi];
		alt_long_array_CI[(qxnpts-iqx-1) * qynpts + (qynpts-iqy-1)][ipT * n_pphi_pts + ipphi] = -alt_long_array_CI[iqx * qynpts + iqy][ipT * n_pphi_pts + ipphi];	//N.B. - sin_trans* odd
		alt_long_array_SR[(qxnpts-iqx-1) * qynpts + (qynpts-iqy-1)][ipT * n_pphi_pts + ipphi] = -alt_long_array_SR[iqx * qynpts + iqy][ipT * n_pphi_pts + ipphi];	//N.B. - sin_trans* odd
		alt_long_array_SI[(qxnpts-iqx-1) * qynpts + (qynpts-iqy-1)][ipT * n_pphi_pts + ipphi] = alt_long_array_SI[iqx * qynpts + iqy][ipT * n_pphi_pts + ipphi];
	}

	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	{
		current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,0)] = alt_long_array_CR[iqx * qynpts + iqy][ipT * n_pphi_pts + ipphi];		//CL,CT
		current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,1)] = alt_long_array_CI[iqx * qynpts + iqy][ipT * n_pphi_pts + ipphi];		//CL,ST
		//current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,2)] = alt_long_array_SR[iqx * qynpts + iqy][ipT * n_pphi_pts + ipphi];
		//current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,3)] = alt_long_array_SI[iqx * qynpts + iqy][ipT * n_pphi_pts + ipphi];
		current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,2)] = alt_long_array_SI[iqx * qynpts + iqy][ipT * n_pphi_pts + ipphi];		//SL,CT
		current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,3)] = alt_long_array_SR[iqx * qynpts + iqy][ipT * n_pphi_pts + ipphi];		//SL,ST
	}
	//////////
	//////////

	delete [] expK0_Bessel_re;
	delete [] expK0_Bessel_im;
	delete [] expK1_Bessel_re;
	delete [] expK1_Bessel_im;

	for (int isa = 0; isa < qxnpts * qynpts; ++isa)
	{
		delete [] alt_long_array_CR[isa];
		delete [] alt_long_array_CI[isa];
		delete [] alt_long_array_SR[isa];
		delete [] alt_long_array_SI[isa];
	}
	delete [] alt_long_array_CR;
	delete [] alt_long_array_CI;
	delete [] alt_long_array_SR;
	delete [] alt_long_array_SI;

	sw.Stop();
	*global_out_stream_ptr << "Total function call took " << sw.printTime() << " seconds." << endl;
	
	return;
}

void CorrelationFunction::Cal_dN_dypTdpTdphi_no_weights_Yeq0_alternate()
{
	//pions
	int local_pid = target_particle_id;

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
	if (use_delta_f)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	//spatial rapidity grid
	//const int eta_s_npts = 31;
	double * eta_s = new double [eta_s_npts];
	double * eta_s_weight = new double [eta_s_npts];
	//double eta_s_i = 0.0, eta_s_f = 4.0;
	gauss_quadrature(eta_s_npts, 1, 0.0, 0.0, eta_s_i, eta_s_f, eta_s, eta_s_weight);
	double * ch_eta_s = new double [eta_s_npts];
	double * sh_eta_s = new double [eta_s_npts];
	for (int ieta = 0; ieta < eta_s_npts; ieta++)
	{
		ch_eta_s[ieta] = cosh(eta_s[ieta]);
		sh_eta_s[ieta] = sinh(eta_s[ieta]);
	}
	// set the rapidity-integration symmetry factor
	double eta_even_factor = 2.0;

	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	{
		double pT = SP_pT[ipT];
		double pphi = SP_pphi[ipphi];
		double px = pT*cos_SP_pphi[ipphi];
		double py = pT*sin_SP_pphi[ipphi];
		double mT = sqrt(pT*pT+localmass*localmass);

		double spectra_at_pTpphi = 0.0;

		for (int isurf = 0; isurf < FO_length; ++isurf)
		{
			FO_surf*surf = &FOsurf_ptr[isurf];

			double tau = surf->tau;
			double xpt = surf->xpt;
			double ypt = surf->ypt;
			double rpt = surf->r;
			double phipt = place_in_range(surf->phi, SP_pphi_min, SP_pphi_max);
			double ch_eta_t = cosh(eta_t(rpt));
			double sh_eta_t = sinh(eta_t(rpt));

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

			for (int ieta = 0; ieta < eta_s_npts; ++ieta)
			{
				double p0 = mT*ch_eta_s[ieta];
				double pz = mT*sh_eta_s[ieta];

				double f0 = ( USE_EXACT ) ? Hfactor(rpt, tau)*exp( one_by_Tdec*pT*sh_eta_t*cos(phipt - pphi) )
												* exp( -one_by_Tdec*p0*ch_eta_t )
							: 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
				//viscous corrections
				double deltaf = 0.;
				if (use_delta_f && !USE_EXACT)
					deltaf = deltaf_prefactor * (1. - sign*f0)
								* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

				//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
				double S_p_with_weight = ( USE_EXACT ) ? eta_s_weight[ieta]*tau*prefactor*da0*f0
											: eta_s_weight[ieta]*tau*prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

				//ignore points where delta f is large or emission function goes negative from pdsigma
				if ( (1. + deltaf < 0.0) || (flagneg == 1 && S_p_with_weight < tol) )
				{
					S_p_with_weight = 0.0;
					continue;
				}

				spectra_at_pTpphi += eta_even_factor * S_p_with_weight;
			}
		}

		//update spectra
		spectra[local_pid][ipT][ipphi] = spectra_at_pTpphi;
		thermal_spectra[local_pid][ipT][ipphi] = spectra_at_pTpphi;
		log_spectra[local_pid][ipT][ipphi] = log(abs(spectra_at_pTpphi) + 1.e-100);
		sign_spectra[local_pid][ipT][ipphi] = sgn(spectra_at_pTpphi);
	}

	//clean up
	delete [] eta_s;
	delete [] eta_s_weight;
	delete [] ch_eta_s;
	delete [] sh_eta_s;

	return;
}


void CorrelationFunction::Cal_dN_dypTdpTdphi_with_weights_Yeq0_alternate(int iqt, int iqz)
{
	//pions
	int local_pid = target_particle_id;

	//q-point info
	double qt = qt_pts[iqt], qz = qz_pts[iqz];

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
	if (use_delta_f)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	//spatial rapidity grid
	//const int eta_s_npts = 101;
	//const int eta_s_npts = 31;	//probably close enough...
	double * eta_s = new double [eta_s_npts];
	double * eta_s_weight = new double [eta_s_npts];
	//double eta_s_i = 0.0, eta_s_f = 4.0;
	gauss_quadrature(eta_s_npts, 1, 0.0, 0.0, eta_s_i, eta_s_f, eta_s, eta_s_weight);
	double * ch_eta_s = new double [eta_s_npts];
	double * sh_eta_s = new double [eta_s_npts];
	for (int ieta = 0; ieta < eta_s_npts; ieta++)
	{
		ch_eta_s[ieta] = cosh(eta_s[ieta]);
		sh_eta_s[ieta] = sinh(eta_s[ieta]);
	}

	for (int isurf = 0; isurf < FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double xpt = surf->xpt;
		double ypt = surf->ypt;
		double rpt = surf->r;
		double phipt = place_in_range(surf->phi, SP_pphi_min, SP_pphi_max);
		double ch_eta_t = cosh(eta_t(rpt));
		double sh_eta_t = sinh(eta_t(rpt));

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

		double * tmpX = oscx[isurf];
		double * tmpY = oscy[isurf];

		for (int ieta = 0; ieta < eta_s_npts; ++ieta)
		{
			double tpt = tau*ch_eta_s[ieta];
			double zpt = tau*sh_eta_s[ieta];
			double phi_L = (tpt*qt-zpt*qz)/hbarC;
			double phi_L_mz = (tpt*qt+zpt*qz)/hbarC;
			double cos_phi_L = cos(phi_L) + cos(phi_L_mz);	//shortcut for eta_s-integral
			double sin_phi_L = sin(phi_L) + sin(phi_L_mz);	//shortcut for eta_s-integral

			for (int ipT = 0; ipT < n_pT_pts; ++ipT)
			for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
			{
				double pT = SP_pT[ipT];
				double pphi = SP_pphi[ipphi];
				double px = pT*cos_SP_pphi[ipphi];
				double py = pT*sin_SP_pphi[ipphi];
				double mT = sqrt(pT*pT+localmass*localmass);
				double p0 = mT*ch_eta_s[ieta];
				double pz = mT*sh_eta_s[ieta];

				double f0 = ( USE_EXACT ) ? Hfactor(rpt, tau)*exp( one_by_Tdec*pT*sh_eta_t*cos(phipt - pphi) )
												* exp( -one_by_Tdec*p0*ch_eta_t )
							: 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
				//viscous corrections
				double deltaf = 0.;
				if (use_delta_f && !USE_EXACT)
					deltaf = deltaf_prefactor * (1. - sign*f0)
								* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

				//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
				double S_p_with_weight = ( USE_EXACT ) ? eta_s_weight[ieta]*tau*prefactor*da0*f0
											: eta_s_weight[ieta]*tau*prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

				//ignore points where delta f is large or emission function goes negative from pdsigma
				if ( (1. + deltaf < 0.0) || (flagneg == 1 && S_p_with_weight < tol) )
				{
					S_p_with_weight = 0.0;
					continue;
				}

				for (int iqx = 0; iqx < qxnpts; ++iqx)
				for (int iqy = 0; iqy < qynpts; ++iqy)
				{
					double cosAx = tmpX[iqx * 2 + 0], sinAx = tmpX[iqx * 2 + 1];
					double cosAy = tmpY[iqy * 2 + 0], sinAy = tmpY[iqy * 2 + 1];
					double cos_trans_Fourier = cosAx*cosAy - sinAx*sinAy;	//==cos(qx x + qy y)
					double sin_trans_Fourier = sinAx*cosAy + cosAx*sinAy;	//==sin(qx x + qy y)
					thermal_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, 0)] += S_p_with_weight * cos_phi_L * cos_trans_Fourier;
					thermal_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, 1)] -= S_p_with_weight * cos_phi_L * sin_trans_Fourier;
					thermal_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, 2)] += S_p_with_weight * sin_phi_L * cos_trans_Fourier;
					thermal_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, 3)] += S_p_with_weight * sin_phi_L * sin_trans_Fourier;
				}
			}
		}
	}

	//clean up
	delete [] eta_s;
	delete [] eta_s_weight;
	delete [] ch_eta_s;
	delete [] sh_eta_s;

	return;
}


//only works at Y==0!!!
void CorrelationFunction::Reflect_in_qz_and_qt()
{
cout << "Currently reflecting in qz and qt!" << endl;
	//for (int iqt = 0; iqt < (qtnpts-1)/2; ++iqt)
	for (int iqt = 0; iqt < (qtnpts+1)/2; ++iqt)
	{
		//reflect in qz first
		//for (int iqz = 0; iqz < (qznpts-1)/2; ++iqz)
		for (int iqz = 0; iqz < (qznpts+1)/2; ++iqz)
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int itrig = 0; itrig < ntrig; ++itrig)
		{
			thermal_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, qznpts - iqz - 1, itrig)]
				= thermal_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, itrig)];			//alternating prefactor still consistent with ntrig == 4
			full_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, qznpts - iqz - 1, itrig)]
				= full_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, itrig)];				//alternating prefactor still consistent with ntrig == 4
			//cout << "Reflection(qz): " << thermal_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, qznpts - iqz - 1, itrig)] << "   "
			//		<< thermal_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, itrig)] << "   "
			//		<< full_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, qznpts - iqz - 1, itrig)] << "   "
			//		<< full_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, itrig)] << endl;
		}
	}

	double trig_parities[4] = {1.0, -1.0, -1.0, 1.0};

	//for (int iqt = 0; iqt < (qtnpts-1)/2; ++iqt)
	for (int iqt = 0; iqt < (qtnpts+1)/2; ++iqt)
	{
		//then reflect everything in qt
		for (int iqz = 0; iqz < qznpts; ++iqz)
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int itrig = 0; itrig < ntrig; ++itrig)
		{
			thermal_target_Yeq0_moments[indexer(ipT, ipphi, qtnpts - iqt - 1, qxnpts - iqx - 1, qynpts - iqy - 1, qznpts - iqz - 1, itrig)]
				= trig_parities[itrig] * thermal_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, itrig)];					//alternating prefactor still consistent with ntrig == 4
			full_target_Yeq0_moments[indexer(ipT, ipphi, qtnpts - iqt - 1, qxnpts - iqx - 1, qynpts - iqy - 1, qznpts - iqz - 1, itrig)]
				= trig_parities[itrig] * full_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, itrig)];						//alternating prefactor still consistent with ntrig == 4
			//cout << "Reflection(qt): " << thermal_target_Yeq0_moments[indexer(ipT, ipphi, qtnpts - iqt - 1, qxnpts - iqx - 1, qynpts - iqy - 1, qznpts - iqz - 1, itrig)] << "   "
			//		<< trig_parities[itrig] * thermal_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, itrig)] << "   "
			//		<< full_target_Yeq0_moments[indexer(ipT, ipphi, qtnpts - iqt - 1, qxnpts - iqx - 1, qynpts - iqy - 1, qznpts - iqz - 1, itrig)] << "   "
			//		<< trig_parities[itrig] * full_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, itrig)] << endl;
		}
	}

cout << "Finished reflection in qz and qt!" << endl;

	return;
}



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
double CorrelationFunction::S_x_p( int local_pid,
									int isurf, double eta_s,
									double pT, double pphi, double pY )
{
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
	if (use_delta_f)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	FO_surf * surf = &FOsurf_ptr[isurf];

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

	double px = pT*cos(pphi);
	double py = pT*sin(pphi);
	double mT = sqrt(pT*pT+localmass*localmass);
	double p0 = mT*cosh(pY-eta_s);
	double pz = mT*sinh(pY-eta_s);	//double-check sign, but okay for now

	double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
	//viscous corrections
	double deltaf = 0.;
	if (use_delta_f)
		deltaf = deltaf_prefactor * (1. - sign*f0)
					* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

	//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
	double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

	//ignore points where delta f is large or emission function goes negative from pdsigma
	if ( (1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol) )
		S_p = 0.0;

	return (S_p);
}


//End of file
