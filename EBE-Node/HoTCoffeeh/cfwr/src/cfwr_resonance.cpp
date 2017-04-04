#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>
#include<time.h>

#include "cfwr.h"
#include "cfwr_lib.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"

using namespace std;

const double PTCHANGE = 1.0;

int local_verbose = 0;

template < typename T >
void check_for_NaNs(string variable_name, const T variable_value, ofstream& localout)
{
	if (isnan(variable_value))
		localout << "ERROR: " << variable_name << " = " << variable_value << endl;
	return;
}

double CorrelationFunction::get_Q()
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

double CorrelationFunction::g(double s)
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

void CorrelationFunction::Do_resonance_integrals(int parent_resonance_particle_id, int daughter_particle_id, int decay_channel)
{
	time_t rawtime;
  	struct tm * timeinfo;

	int daughter_lookup_idx = distance(daughter_resonance_indices.begin(), daughter_resonance_indices.find(daughter_particle_id));

	Allocate_resonance_running_sum_vectors();

	//set these for quick look-up in EdNd3p() routine
	spec_vals_info = spectra[parent_resonance_particle_id];
	spec_log_info = log_spectra[parent_resonance_particle_id];
	spec_sign_info = sign_spectra[parent_resonance_particle_id];

	Flatten_dN_dypTdpTdphi_moments(parent_resonance_particle_id);

	int tmp_parent_monval = all_particles[parent_resonance_particle_id].monval;
	int tmp_daughter_monval = all_particles[daughter_particle_id].monval;
	n_body = current_reso_nbody;
	current_parent_resonance = parent_resonance_particle_id;

	local_verbose = 0;

	int iqt0 = (qtnpts-1)/2;
	int iqx0 = (qxnpts-1)/2;
	int iqy0 = (qynpts-1)/2;
	int iqz0 = (qznpts-1)/2;
	int tmp_qpt_cs_idx = ( ( ( iqt0 * qxnpts + iqx0 ) * qynpts + iqy0 ) * qznpts + iqz0 ) * 2;

	if (n_body == 2)
	{
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		for (int ipY = 0; ipY < n_pY_pts; ++ipY)
		{
			double local_pT = SP_pT[ipT];
			double local_pphi = SP_pphi[ipphi];
			double local_pY = SP_pY[ipphi];
			current_ipT = ipT;
			current_ipphi = ipphi;
			current_ipY = ipY;
			//current_qlist_slice = qlist[ipT];
			Zero_resonance_running_sum_vector(ssum_vec);
			Zero_resonance_running_sum_vector(vsum_vec);
			Zero_resonance_running_sum_vector(zetasum_vec);
			Zero_resonance_running_sum_vector(Csum_vec);
			Load_decay_channel_info(decay_channel, local_pT, local_pphi, local_pY);	// set decay channel information

			//then g(s) is delta-function, skip s-integration entirely
			//double s_loc = m2*m2;
			double vsum = 0.0;
			for (int iv = 0; iv < n_v_pts; ++iv)
			{
				//time (&rawtime);
				//timeinfo = localtime (&rawtime);
				Zero_resonance_running_sum_vector(zetasum_vec);
				double zetasum = 0.0;
				for (int izeta = 0; izeta < n_zeta_pts; ++izeta)
				{
					Zero_resonance_running_sum_vector(Csum_vec);
					double Csum = 0.0;
					double PKT = VEC_n2_PT[iv][izeta];
					double PKY = VEC_n2_P_Y[iv];
					double PKphi = VEC_n2_PPhi_tilde[iv][izeta];
					for (int tempidx = 0; tempidx <= 1; ++tempidx)
					{
						if (tempidx != 0)
							PKphi = VEC_n2_PPhi_tildeFLIP[iv][izeta];		//also takes Pp --> Pm
						currentPpm = VEC_n2_Ppm[iv][izeta][tempidx];
						Edndp3(PKT, PKphi, &Csum);							//set spectra
						if (!IGNORE_LONG_LIVED_RESONANCES || Gamma >= hbarC / max_lifetime)
							eiqxEdndp3(PKT, PKphi, Csum_vec, local_verbose);					//set weights
					}												// end of tempidx sum
					for (int qpt_cs_idx = 0; qpt_cs_idx < qspace_cs_slice_length; ++qpt_cs_idx)
						zetasum_vec[qpt_cs_idx] += VEC_n2_zeta_factor[iv][izeta]*Csum_vec[qpt_cs_idx];
					zetasum += VEC_n2_zeta_factor[iv][izeta]*Csum;
				}													// end of zeta sum
				for (int qpt_cs_idx = 0; qpt_cs_idx < qspace_cs_slice_length; ++qpt_cs_idx)
					vsum_vec[qpt_cs_idx] += VEC_n2_v_factor[iv]*zetasum_vec[qpt_cs_idx];
				vsum += VEC_n2_v_factor[iv]*zetasum;
			}														// end of v sum
			for (int qpt_cs_idx = 0; qpt_cs_idx < qspace_cs_slice_length; ++qpt_cs_idx)
				ssum_vec[qpt_cs_idx] += Mres*VEC_n2_s_factor*vsum_vec[qpt_cs_idx];
			double ssum = Mres*VEC_n2_s_factor*vsum;

			//update all gridpoints for all daughter moments
			int qpt_cs_idx = 0;
			for (int iqt = 0; iqt < qtnpts; ++iqt)
			for (int iqx = 0; iqx < qxnpts; ++iqx)
			for (int iqy = 0; iqy < qynpts; ++iqy)
			for (int iqz = 0; iqz < qznpts; ++iqz)
			for (int itrig = 0; itrig < 2; ++itrig)
			{
				current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipT,ipphi,ipY,iqt,iqx,iqy,iqz,itrig)] += ssum_vec[qpt_cs_idx];
				++qpt_cs_idx;
			}

			//update daughter spectra separately
			spectra[daughter_particle_id][ipT][ipphi] += ssum;
			log_spectra[daughter_particle_id][ipT][ipphi] = log(abs(spectra[daughter_particle_id][ipT][ipphi])+1.e-100);
			sign_spectra[daughter_particle_id][ipT][ipphi] = sgn(spectra[daughter_particle_id][ipT][ipphi]);

//cerr << "CHECK: " << daughter_lookup_idx << "   " << daughter_particle_id << "   " << parent_resonance_particle_id << "   " << decay_channel
//		<< "   " << ipT << "   " << ipphi << "   " << ssum_vec[0] << "   " << ssum_vec[1] << "   " << ssum
//		<< "   " << spectra[daughter_particle_id][ipT][ipphi] << "   " << current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipT,ipphi,0,0,0,0,0)] << "   " << endl;

			//only do this if all qpoints array sizes are odd!
			//now, if ignoring long-lived resonances, take them out of the correlator numerator, but keep them in the denominator (AKA, spectra)
			//////////////////if (IGNORE_LONG_LIVED_RESONANCES && qtnpts%2==1 && qxnpts%2==1 && qynpts%2==1 && qznpts%2==1 && Gamma < hbarC / max_lifetime)
			//////////////////	current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipT,ipphi,iqt0,iqx0,iqy0,iqz0,0)] -= ssum_vec[tmp_qpt_cs_idx];
	
			if (isnan(current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipT,ipphi,ipY,0,0,0,0,0)]
					+ current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipT,ipphi,ipY,0,0,0,0,1)]))
			{
				*global_out_stream_ptr << "ERROR: NaNs encountered!" << endl
										<< "current_daughters_dN_dypTdpTdphi_moments[" << daughter_lookup_idx << "][" << ipT << "][" << ipphi << "][" << ipY << "][0][0][0][0][0] = "
										<< setw(8) << setprecision(15)
										<< current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipT,ipphi,ipY,0,0,0,0,0)] << endl
										<< "current_daughters_dN_dypTdpTdphi_moments[" << daughter_lookup_idx << "][" << ipT << "][" << ipphi << "][" << ipY << "][0][0][0][0][1] = "
										<< setw(8) << setprecision(15)
										<< current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipT,ipphi,ipY,0,0,0,0,1)] << endl
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
				if (local_verbose == 0) exit(1);
			}
		}											// end of pT, pphi loops
	}												// end of nbody == 2
	else
	{
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		for (int ipY = 0; ipY < n_pY_pts; ++ipY)
		{
			double local_pT = SP_pT[ipT];
			double local_pphi = SP_pphi[ipphi];
			double local_pY = SP_pphi[ipY];
			current_ipT = ipT;
			current_ipphi = ipphi;
			current_ipphi = ipY;
			//current_qlist_slice = qlist[ipT];
			Zero_resonance_running_sum_vector(ssum_vec);
			Zero_resonance_running_sum_vector(vsum_vec);
			Zero_resonance_running_sum_vector(zetasum_vec);
			Zero_resonance_running_sum_vector(Csum_vec);
			Load_decay_channel_info(decay_channel, local_pT, local_pphi, local_pY);	// set decay channel information

			double ssum = 0.0;
			for (int is = 0; is < n_s_pts; ++is)
			{
				double vsum = 0.0;
 		  		Zero_resonance_running_sum_vector(vsum_vec);
				for (int iv = 0; iv < n_v_pts; ++iv)
				{
					double zetasum = 0.0;
					Zero_resonance_running_sum_vector(zetasum_vec);
					for (int izeta = 0; izeta < n_zeta_pts; ++izeta)
					{
						double Csum = 0.0;
						Zero_resonance_running_sum_vector(Csum_vec);
						double PKT = VEC_PT[is][iv][izeta];
						double PKY = VEC_P_Y[is][iv];
						double PKphi = VEC_PPhi_tilde[is][iv][izeta];
						for (int tempidx = 0; tempidx <= 1; ++tempidx)
						{
							if (tempidx != 0)
								PKphi = VEC_PPhi_tildeFLIP[is][iv][izeta];		//also takes Pp --> Pm
							currentPpm = VEC_Ppm[is][iv][izeta][tempidx];
							Edndp3(PKT, PKphi, &Csum);							//set spectra
							if (!IGNORE_LONG_LIVED_RESONANCES || Gamma >= hbarC / max_lifetime)
								eiqxEdndp3(PKT, PKphi, Csum_vec, local_verbose);					//set weights
						}										// end of tempidx sum
						for (int qpt_cs_idx = 0; qpt_cs_idx < qspace_cs_slice_length; ++qpt_cs_idx)
							zetasum_vec[qpt_cs_idx] += VEC_zeta_factor[is][iv][izeta]*Csum_vec[qpt_cs_idx];
						zetasum += VEC_zeta_factor[is][iv][izeta]*Csum;
					}											// end of zeta sum
					for (int qpt_cs_idx = 0; qpt_cs_idx < qspace_cs_slice_length; ++qpt_cs_idx)
					    vsum_vec[qpt_cs_idx] += VEC_v_factor[is][iv]*zetasum_vec[qpt_cs_idx];
					vsum += VEC_v_factor[is][iv]*zetasum;
				}												// end of v sum
				for (int qpt_cs_idx = 0; qpt_cs_idx < qspace_cs_slice_length; ++qpt_cs_idx)
					ssum_vec[qpt_cs_idx] += Mres*VEC_s_factor[is]*vsum_vec[qpt_cs_idx];
				ssum += Mres*VEC_s_factor[is]*vsum;
			}													// end of s sum
			//update all gridpoints for daughter moments
			int qpt_cs_idx = 0;
			for (int iqt = 0; iqt < qtnpts; ++iqt)
			for (int iqx = 0; iqx < qxnpts; ++iqx)
			for (int iqy = 0; iqy < qynpts; ++iqy)
			for (int iqz = 0; iqz < qznpts; ++iqz)
			for (int itrig = 0; itrig < 2; ++itrig)
			{
				current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipT,ipphi,ipY,iqt,iqx,iqy,iqz,itrig)] += ssum_vec[qpt_cs_idx];
				++qpt_cs_idx;
			}

			//update daughter spectra separately
			spectra[daughter_particle_id][ipT][ipphi] += ssum;
			log_spectra[daughter_particle_id][ipT][ipphi] = log(abs(spectra[daughter_particle_id][ipT][ipphi])+1.e-100);
			sign_spectra[daughter_particle_id][ipT][ipphi] = sgn(spectra[daughter_particle_id][ipT][ipphi]);

//cerr << "CHECK: " << daughter_lookup_idx << "   " << daughter_particle_id << "   " << parent_resonance_particle_id << "   " << decay_channel
//		<< "   " << ipT << "   " << ipphi << "   " << ssum_vec[0] << "   " << ssum_vec[1] << "   " << ssum
//		<< "   " << spectra[daughter_particle_id][ipT][ipphi] << "   " << current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipT,ipphi,0,0,0,0,0)] << endl;

			//only do this if all qpoints array sizes are odd!
			//now, if ignoring long-lived resonances, take them out of the correlator numerator, but keep them in the denominator (AKA, spectra)
			/////////////////if (IGNORE_LONG_LIVED_RESONANCES && qtnpts%2==1 && qxnpts%2==1 && qynpts%2==1 && qznpts%2==1 && Gamma < hbarC / max_lifetime)
			/////////////////	current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipT,ipphi,iqt0,iqx0,iqy0,iqz0,0)] -= ssum_vec[tmp_qpt_cs_idx];

			if (isnan(current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipT,ipphi,ipY,0,0,0,0,0)]
					+ current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipT,ipphi,ipY,0,0,0,0,1)]))
			{
				*global_out_stream_ptr << "ERROR: NaNs encountered!" << endl
										<< "current_daughters_dN_dypTdpTdphi_moments[" << daughter_lookup_idx << "][" << ipT << "][" << ipphi << "][" << ipY << "][0][0][0][0][0] = "
										<< setw(8) << setprecision(15)
										<< current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipT,ipphi,ipY,0,0,0,0,0)] << endl
										<< "current_daughters_dN_dypTdpTdphi_moments[" << daughter_lookup_idx << "][" << ipT << "][" << ipphi << "][" << ipY << "][0][0][0][0][1] = "
										<< setw(8) << setprecision(15)
										<< current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipT,ipphi,ipY,0,0,0,0,1)] << endl
										<< "  --> pt = " << local_pT << endl
										<< "  --> pphi = " << local_pphi << endl
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
				if (local_verbose == 0) exit(1);
			}
		}								// end of pT, pphi loops
	}										// end of nbody == 3

	// clean up
	Delete_resonance_running_sum_vectors();

	for (int it = 0; it < (int)spectra_resonance_grid_approximator.size(); ++it)
		delete spectra_resonance_grid_approximator[it];
	for (int it = 0; it < (int)real_resonance_grid_approximator.size(); ++it)
		delete real_resonance_grid_approximator[it];
	for (int it = 0; it < (int)imag_resonance_grid_approximator.size(); ++it)
		delete imag_resonance_grid_approximator[it];

	spectra_resonance_grid_approximator.clear();
	real_resonance_grid_approximator.clear();
	imag_resonance_grid_approximator.clear();

	return;
}

void CorrelationFunction::Flatten_dN_dypTdpTdphi_moments(int parent_resonance_particle_id)
{
	const int dim_loc = 2;
	int npts_loc[dim_loc] = { n_pT_pts, n_pphi_pts };
	int os[dim_loc] = { n_pT_pts-1, n_pphi_pts-1 };
	double lls[dim_loc] = { KT_min, Kphi_min };
	double uls[dim_loc] = { KT_max, Kphi_max };
	int modes_loc[dim_loc] = { 0, 0 };

	int momidx = 0;

	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	for (int ipY = 0; ipY < n_pY_pts; ++ipY)
	{
		// set index for looping
		int qpt_cs_idx = 0;
		
		// flatten arrays, since these are quicker to process
		for (int iqt = 0; iqt < qtnpts; ++iqt)
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int iqz = 0; iqz < qznpts; ++iqz)
		{
			for (int itrig = 0; itrig < 2; ++itrig)
			{
				res_sign_info[ipT][ipphi][qpt_cs_idx] = current_sign_of_dN_dypTdpTdphi_moments[indexer(ipT,ipphi,ipY,iqt,iqx,iqy,iqz,itrig)];
				res_log_info[ipT][ipphi][qpt_cs_idx] = current_ln_dN_dypTdpTdphi_moments[indexer(ipT,ipphi,ipY,iqt,iqx,iqy,iqz,itrig)];
				res_moments_info[ipT][ipphi][qpt_cs_idx] = current_dN_dypTdpTdphi_moments[indexer(ipT,ipphi,ipY,iqt,iqx,iqy,iqz,itrig)];
				++qpt_cs_idx;
			}

			tmp_moments_real[iqt][iqx][iqy][iqz][momidx] = current_dN_dypTdpTdphi_moments[indexer(ipT,ipphi,ipY,iqt,iqx,iqy,iqz,0)];
			tmp_moments_imag[iqt][iqx][iqy][iqz][momidx] = current_dN_dypTdpTdphi_moments[indexer(ipT,ipphi,ipY,iqt,iqx,iqy,iqz,1)];
		}

		flat_spectra[momidx] = spectra[parent_resonance_particle_id][ipT][ipphi];

		momidx++;
	}

	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	{
		double * tmpreal = tmp_moments_real[iqt][iqx][iqy][iqz];
		double * tmpimag = tmp_moments_imag[iqt][iqx][iqy][iqz];
		real_resonance_grid_approximator.push_back( new Chebyshev(tmpreal, npts_loc, os, lls, uls, dim_loc, modes_loc) );
		imag_resonance_grid_approximator.push_back( new Chebyshev(tmpimag, npts_loc, os, lls, uls, dim_loc, modes_loc) );
	}

	spectra_resonance_grid_approximator.push_back( new Chebyshev(flat_spectra, npts_loc, os, lls, uls, dim_loc, modes_loc) );

	return;
}

void CorrelationFunction::Edndp3(double ptr, double phir, double * result, int loc_verb /*==0*/)
{
	double phi0, phi1;
	double f1, f2;

	int npphi_max = n_pphi_pts - 1;
	int npT_max = n_pT_pts - 1;

	// locate pT interval
	int npt = 1;
	while ((ptr > SP_pT[npt]) &&
			(npt < npT_max)) ++npt;
	double pT0 = SP_pT[npt-1];
	double pT1 = SP_pT[npt];

	// locate pphi interval
	int nphi = 1, nphim1 = 0;
	if(phir < SP_pphi[0])			//if angle is less than minimum angle grid point
	{
		phi0 = SP_pphi[npphi_max] - 2. * M_PI;
		phi1 = SP_pphi[0];
		nphi = 0;
		nphim1 = npphi_max;
	}
	else if(phir > SP_pphi[npphi_max])	//if angle is greater than maximum angle grid point
	{
		phi0 = SP_pphi[npphi_max];
		phi1 = SP_pphi[0] + 2. * M_PI;
		nphi = 0;
		nphim1 = npphi_max;
	}
	else						//if angle is within grid range
	{
		while ((phir > SP_pphi[nphi]) &&
				(nphi < npphi_max)) ++nphi;
		nphim1 = nphi - 1;
		phi0 = SP_pphi[nphim1];
		phi1 = SP_pphi[nphi];
	}

	if (pT0==pT1 || phi0==phi1)
	{
		cerr << "ERROR in Edndp3(): pT and/or pphi values equal!" << endl;
		exit(1);
	}

	double one_by_pTdiff = 1./(pT1 - pT0), one_by_pphidiff = 1./(phi1 - phi0);

	// choose pt-pphi slice of resonance info arrays
	double log_f11 = spec_log_info[npt-1][nphim1];
	double log_f12 = spec_log_info[npt-1][nphi];
	double log_f21 = spec_log_info[npt][nphim1];
	double log_f22 = spec_log_info[npt][nphi];

	double f11 = spec_vals_info[npt-1][nphim1];
	double f12 = spec_vals_info[npt-1][nphi];
	double f21 = spec_vals_info[npt][nphim1];
	double f22 = spec_vals_info[npt][nphi];

	/////////////////////////////////////////////////////////////////
	// interpolate over pT values first
	/////////////////////////////////////////////////////////////////
	if(ptr > PTCHANGE)				// if pT interpolation point is larger than PTCHANGE (currently 1.0 GeV)
	{
		double sign_of_f11 = spec_sign_info[npt-1][nphim1];
		double sign_of_f12 = spec_sign_info[npt-1][nphi];
		double sign_of_f21 = spec_sign_info[npt][nphim1];
		double sign_of_f22 = spec_sign_info[npt][nphi];
		
		//*******************************************************************************************************************
		// set f1 first
		//*******************************************************************************************************************
		// if using extrapolation and spectra at pT1 has larger magnitude than at pT0 (or the signs are different), just return zero
		if ( ptr > pT1 && ( log_f21 > log_f11 || sign_of_f11 * sign_of_f21 < 0 ) )
			f1 = 0.0;
		else if (sign_of_f11 * sign_of_f21 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
			f1 = sign_of_f11 * exp( lin_int(ptr-pT0, one_by_pTdiff, log_f11, log_f21) );
		else					// otherwise, just interpolate original vals
			f1 = lin_int(ptr-pT0, one_by_pTdiff, f11, f21);
			
		//*******************************************************************************************************************
		// set f2 next
		//*******************************************************************************************************************
		if ( ptr > pT1 && ( log_f22 > log_f12 || sign_of_f12 * sign_of_f22 < 0 ) )
			f2 = 0.0;
		else if (sign_of_f12 * sign_of_f22 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
			f2 = sign_of_f12 * exp( lin_int(ptr-pT0, one_by_pTdiff, log_f12, log_f22) );
		else					// otherwise, just interpolate original vals
			f2 = lin_int(ptr-pT0, one_by_pTdiff, f12, f22);
		//*******************************************************************************************************************
	}
	else						// if pT is smaller than PTCHANGE, just use linear interpolation, no matter what
	{
		f1 = lin_int(ptr-pT0, one_by_pTdiff, f11, f21);
		f2 = lin_int(ptr-pT0, one_by_pTdiff, f12, f22);
	}
				
	// now, interpolate f1 and f2 over the pphi direction
	*result += lin_int(phir-phi0, one_by_pphidiff, f1, f2);

	return;
}

void CorrelationFunction::eiqxEdndp3(double ptr, double phir, double * results, int loc_verb /*==0*/)
{
	double phi0, phi1;
	double f1, f2;

	int npphi_max = n_pphi_pts - 1;
	int npT_max = n_pT_pts - 1;

	// locate pT interval
	int npt = 1;
	while ((ptr > SP_pT[npt]) &&
			(npt < npT_max)) ++npt;
	double pT0 = SP_pT[npt-1];
	double pT1 = SP_pT[npt];

	// locate pphi interval
	int nphi = 1, nphim1 = 0;
	if(phir < SP_pphi[0])			//if angle is less than minimum angle grid point
	{
		phi0 = SP_pphi[npphi_max] - 2. * M_PI;
		phi1 = SP_pphi[0];
		nphi = 0;
		nphim1 = npphi_max;
	}
	else if(phir > SP_pphi[npphi_max])	//if angle is greater than maximum angle grid point
	{
		phi0 = SP_pphi[npphi_max];
		phi1 = SP_pphi[0] + 2. * M_PI;
		nphi = 0;
		nphim1 = npphi_max;
	}
	else						//if angle is within grid range
	{
		while ((phir > SP_pphi[nphi]) &&
				(nphi < npphi_max)) ++nphi;
		nphim1 = nphi - 1;
		phi0 = SP_pphi[nphim1];
		phi1 = SP_pphi[nphi];
	}

	if (pT0==pT1 || phi0==phi1)
	{
		cerr << "ERROR in eiqxEdndp3(): pT and/or pphi values equal!" << endl;
		exit(1);
	}

	double one_by_pTdiff = 1./(pT1 - pT0), one_by_pphidiff = 1./(phi1 - phi0);
	double del_ptr_pt0 = ptr - pT0, del_phir_phi0 = phir - phi0;

	// choose pt-pphi slice of resonance info arrays
	double * sign_of_f11_arr = res_sign_info[npt-1][nphim1];
	double * sign_of_f12_arr = res_sign_info[npt-1][nphi];
	double * sign_of_f21_arr = res_sign_info[npt][nphim1];
	double * sign_of_f22_arr = res_sign_info[npt][nphi];

	double * log_f11_arr = res_log_info[npt-1][nphim1];
	double * log_f12_arr = res_log_info[npt-1][nphi];
	double * log_f21_arr = res_log_info[npt][nphim1];
	double * log_f22_arr = res_log_info[npt][nphi];

	double * f11_arr = res_moments_info[npt-1][nphim1];
	double * f12_arr = res_moments_info[npt-1][nphi];
	double * f21_arr = res_moments_info[npt][nphim1];
	double * f22_arr = res_moments_info[npt][nphi];

	// set index for looping
	int qpt_cs_idx = 0;
	int qlist_idx = 0;

	if (ptr > PTCHANGE)				// if pT interpolation point is larger than PTCHANGE (currently 1.0 GeV)
	{
		for (int iqt = 0; iqt < qtnpts; ++iqt)
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int iqz = 0; iqz < qznpts; ++iqz)
		{
			double arg = one_by_Gamma_Mres * dot_four_vectors(qlist[qlist_idx], currentPpm);
			double akr = 1./(1.+arg*arg);
			double aki = arg/(1.+arg*arg);

			/////////////////////////////////////////////////////////////////
			// DO COSINE PART FIRST
			/////////////////////////////////////////////////////////////////
			// interpolate over pT values first
			/////////////////////////////////////////////////////////////////
			double sign_of_f11 = sign_of_f11_arr[qpt_cs_idx];
			double sign_of_f12 = sign_of_f12_arr[qpt_cs_idx];
			double sign_of_f21 = sign_of_f21_arr[qpt_cs_idx];
			double sign_of_f22 = sign_of_f22_arr[qpt_cs_idx];
			
			//*******************************************************************************************************************
			// set f1 first
			//*******************************************************************************************************************
			// if using extrapolation and spectra at pT1 has larger magnitude than at pT0 (or the signs are different), just return zero
			if ( ptr > pT1 && ( log_f21_arr[qpt_cs_idx] > log_f11_arr[qpt_cs_idx] || sign_of_f11 * sign_of_f21 < 0 ) )
				f1 = 0.0;
			else if (sign_of_f11 * sign_of_f21 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
				f1 = sign_of_f11 * exp( lin_int(del_ptr_pt0, one_by_pTdiff, log_f11_arr[qpt_cs_idx], log_f21_arr[qpt_cs_idx]) );
			else					// otherwise, just interpolate original vals
				f1 = lin_int(del_ptr_pt0, one_by_pTdiff, f11_arr[qpt_cs_idx], f21_arr[qpt_cs_idx]);
				
			//*******************************************************************************************************************
			// set f2 next
			//*******************************************************************************************************************
			if ( ptr > pT1 && ( log_f22_arr[qpt_cs_idx] > log_f12_arr[qpt_cs_idx] || sign_of_f12 * sign_of_f22 < 0 ) )
				f2 = 0.0;
			else if (sign_of_f12 * sign_of_f22 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
				f2 = sign_of_f12 * exp( lin_int(del_ptr_pt0, one_by_pTdiff, log_f12_arr[qpt_cs_idx], log_f22_arr[qpt_cs_idx]) );
			else					// otherwise, just interpolate original vals
				f2 = lin_int(del_ptr_pt0, one_by_pTdiff, f12_arr[qpt_cs_idx], f22_arr[qpt_cs_idx]);
			//*******************************************************************************************************************
			// now, interpolate f1 and f2 over the pphi direction
			double Zkr = lin_int(del_phir_phi0, one_by_pphidiff, f1, f2);

			double f1s = f1, f2s = f2;
	
			/////////////////////////////////////////////////////////////////
			// DO SINE PART NEXT
			/////////////////////////////////////////////////////////////////
			// interpolate over pT values first
			/////////////////////////////////////////////////////////////////
			sign_of_f11 = sign_of_f11_arr[qpt_cs_idx+1];
			sign_of_f12 = sign_of_f12_arr[qpt_cs_idx+1];
			sign_of_f21 = sign_of_f21_arr[qpt_cs_idx+1];
			sign_of_f22 = sign_of_f22_arr[qpt_cs_idx+1];
			
			//*******************************************************************************************************************
			// set f1 first
			//*******************************************************************************************************************
			// if using extrapolation and spectra at pT1 has larger magnitude than at pT0 (or the signs are different), just return zero
			if ( ptr > pT1 && ( log_f21_arr[qpt_cs_idx+1] > log_f11_arr[qpt_cs_idx+1] || sign_of_f11 * sign_of_f21 < 0 ) )
				f1 = 0.0;
			else if (sign_of_f11 * sign_of_f21 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
				f1 = sign_of_f11 * exp( lin_int(del_ptr_pt0, one_by_pTdiff, log_f11_arr[qpt_cs_idx+1], log_f21_arr[qpt_cs_idx+1]) );
			else					// otherwise, just interpolate original vals
				f1 = lin_int(del_ptr_pt0, one_by_pTdiff, f11_arr[qpt_cs_idx+1], f21_arr[qpt_cs_idx+1]);
				
			//*******************************************************************************************************************
			// set f2 next
			//*******************************************************************************************************************
			if ( ptr > pT1 && ( log_f22_arr[qpt_cs_idx+1] > log_f12_arr[qpt_cs_idx+1] || sign_of_f12 * sign_of_f22 < 0 ) )
				f2 = 0.0;
			else if (sign_of_f12 * sign_of_f22 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
				f2 = sign_of_f12 * exp( lin_int(del_ptr_pt0, one_by_pTdiff, log_f12_arr[qpt_cs_idx+1], log_f22_arr[qpt_cs_idx+1]) );
			else					// otherwise, just interpolate original vals
				f2 = lin_int(del_ptr_pt0, one_by_pTdiff, f12_arr[qpt_cs_idx+1], f22_arr[qpt_cs_idx+1]);
			//*******************************************************************************************************************

			// now, interpolate f1 and f2 over the pphi direction
			double Zki = lin_int(del_phir_phi0, one_by_pphidiff, f1, f2);

			/////////////////////////////////////////////////////
			// Finally, update results vectors appropriately
			/////////////////////////////////////////////////////
			//--> update the real part of weighted daughter spectra
			results[qpt_cs_idx] += akr*Zkr-aki*Zki;
			//--> update the imaginary part of weighted daughter spectra
			results[qpt_cs_idx+1] += akr*Zki+aki*Zkr;
	
			qpt_cs_idx += 2;
			qlist_idx++;
		}	//end of all q-loops
	}
	else						// if pT is smaller than PTCHANGE, just use linear interpolation, no matter what
	{
		for (int iqt = 0; iqt < qtnpts; ++iqt)
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int iqz = 0; iqz < qznpts; ++iqz)
		{
			double arg = one_by_Gamma_Mres * dot_four_vectors(qlist[qlist_idx], currentPpm);
			double akr = 1./(1.+arg*arg);
			double aki = arg/(1.+arg*arg);

			/////////////////////////////////////////////////////////////////
			// DO COSINE PART FIRST
			/////////////////////////////////////////////////////////////////
			// interpolate over pT values first
			/////////////////////////////////////////////////////////////////
			f1 = lin_int(del_ptr_pt0, one_by_pTdiff, f11_arr[qpt_cs_idx], f21_arr[qpt_cs_idx]);
			f2 = lin_int(del_ptr_pt0, one_by_pTdiff, f12_arr[qpt_cs_idx], f22_arr[qpt_cs_idx]);

			// now, interpolate f1 and f2 over the pphi direction
			double Zkr = lin_int(del_phir_phi0, one_by_pphidiff, f1, f2);

			double f1s = f1, f2s = f2;
	
			/////////////////////////////////////////////////////////////////
			// DO SINE PART NEXT
			/////////////////////////////////////////////////////////////////
			// interpolate over pT values first
			/////////////////////////////////////////////////////////////////
			f1 = lin_int(del_ptr_pt0, one_by_pTdiff, f11_arr[qpt_cs_idx+1], f21_arr[qpt_cs_idx+1]);
			f2 = lin_int(del_ptr_pt0, one_by_pTdiff, f12_arr[qpt_cs_idx+1], f22_arr[qpt_cs_idx+1]);

			// now, interpolate f1 and f2 over the pphi direction
			double Zki = lin_int(del_phir_phi0, one_by_pphidiff, f1, f2);

			/////////////////////////////////////////////////////
			// Finally, update results vectors appropriately
			/////////////////////////////////////////////////////
			//--> update the real part of weighted daughter spectra
			results[qpt_cs_idx] += akr*Zkr-aki*Zki;
			//--> update the imaginary part of weighted daughter spectra
			results[qpt_cs_idx+1] += akr*Zki+aki*Zkr;
	
			qpt_cs_idx += 2;
			qlist_idx++;
		}	//end of all q-loops
	}

	return;
}

void CorrelationFunction::Load_decay_channel_info(int dc_idx, double K_T_local, double K_phi_local, double K_y_local)
{
	Mres = current_resonance_mass;
	Gamma = current_resonance_Gamma;
	//one_by_Gamma_Mres = hbarC/(Gamma*Mres + 1.e-25);	//keeps calculation safe when Gamma == 0
	one_by_Gamma_Mres = 1./(Gamma*Mres + 1.e-25);	//keeps calculation safe when Gamma == 0
	//N.B. - no need for hbarc, since this will only multiply something with GeV^2 units in the end
	mass = current_daughter_mass;
	br = current_resonance_direct_br;	//doesn't depend on target daughter particle, just parent resonance and decay channel
	m2 = current_resonance_decay_masses[0];
	m3 = current_resonance_decay_masses[1];

	pT = K_T_local;
	current_K_phi = K_phi_local;
	p_y = K_y_local;

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
		//p_y = 0.0;
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
				double PT_loc = sqrt(MT_loc*MT_loc - Mres*Mres);
				double temp_cos_PPhi_tilde_loc = (mT*MT_loc*cosh(P_Y_loc-p_y) - Estar_loc*Mres)/(pT*PT_loc);
				//assume that PPhi_tilde is +ve in next step...
				double temp_sin_PPhi_tilde_loc = sqrt(1. - temp_cos_PPhi_tilde_loc*temp_cos_PPhi_tilde_loc);
				double PPhi_tilde_loc = place_in_range( atan2(temp_sin_PPhi_tilde_loc, temp_cos_PPhi_tilde_loc), Kphi_min, Kphi_max);
				VEC_n2_PPhi_tilde[iv][izeta] = place_in_range( K_phi_local + PPhi_tilde_loc, Kphi_min, Kphi_max);
				VEC_n2_PPhi_tildeFLIP[iv][izeta] = place_in_range( K_phi_local - PPhi_tilde_loc, Kphi_min, Kphi_max);
				VEC_n2_PT[iv][izeta] = PT_loc;
				//set P^+ components
				VEC_n2_Ppm[iv][izeta][0][0] = MT_loc * cosh(P_Y_loc);
				VEC_n2_Ppm[iv][izeta][0][1] = PT_loc * cos(K_phi_local + PPhi_tilde_loc);
				VEC_n2_Ppm[iv][izeta][0][2] = PT_loc * sin(K_phi_local + PPhi_tilde_loc);
				VEC_n2_Ppm[iv][izeta][0][3] = MT_loc * sinh(P_Y_loc);
				//set P^- components
				VEC_n2_Ppm[iv][izeta][1][0] = MT_loc * cosh(P_Y_loc);
				VEC_n2_Ppm[iv][izeta][1][1] = PT_loc * cos(K_phi_local - PPhi_tilde_loc);
				VEC_n2_Ppm[iv][izeta][1][2] = PT_loc * sin(K_phi_local - PPhi_tilde_loc);
				VEC_n2_Ppm[iv][izeta][1][3] = MT_loc * sinh(P_Y_loc);
				check_for_NaNs("PT_loc", PT_loc, *global_out_stream_ptr);
				check_for_NaNs("PPhi_tilde_loc", PPhi_tilde_loc, *global_out_stream_ptr);
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
			//p_y = 0.0;
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
					double PT_loc = sqrt(MT_loc*MT_loc - Mres*Mres);
					VEC_zeta_factor[is][iv][izeta] = zeta_wts[izeta]*MT_loc;
					double temp_cos_PPhi_tilde_loc = (mT*MT_loc*cosh(P_Y_loc-p_y) - Estar_loc*Mres)/(pT*PT_loc);
					//assume that PPhi_tilde is +ve in next step...
					double temp_sin_PPhi_tilde_loc = sqrt(1. - temp_cos_PPhi_tilde_loc*temp_cos_PPhi_tilde_loc);
					double PPhi_tilde_loc = place_in_range( atan2(temp_sin_PPhi_tilde_loc, temp_cos_PPhi_tilde_loc), Kphi_min, Kphi_max);
					VEC_PPhi_tilde[is][iv][izeta] = place_in_range( K_phi_local + PPhi_tilde_loc, Kphi_min, Kphi_max);
					VEC_PPhi_tildeFLIP[is][iv][izeta] = place_in_range( K_phi_local - PPhi_tilde_loc, Kphi_min, Kphi_max);
					VEC_PT[is][iv][izeta] = PT_loc;
					//set P^+ components
					VEC_Ppm[is][iv][izeta][0][0] = MT_loc * cosh(P_Y_loc);
					VEC_Ppm[is][iv][izeta][0][1] = PT_loc * cos(K_phi_local + PPhi_tilde_loc);
					VEC_Ppm[is][iv][izeta][0][2] = PT_loc * sin(K_phi_local + PPhi_tilde_loc);
					VEC_Ppm[is][iv][izeta][0][3] = MT_loc * sinh(P_Y_loc);
					//set P^- components
					VEC_Ppm[is][iv][izeta][1][0] = MT_loc * cosh(P_Y_loc);
					VEC_Ppm[is][iv][izeta][1][1] = PT_loc * cos(K_phi_local - PPhi_tilde_loc);
					VEC_Ppm[is][iv][izeta][1][2] = PT_loc * sin(K_phi_local - PPhi_tilde_loc);
					VEC_Ppm[is][iv][izeta][1][3] = MT_loc * sinh(P_Y_loc);
				}
			}
		}
	}

	return;
}

//***************************************************************************************************

void CorrelationFunction::Cal_dN_dypTdpTdphi_with_weights_function(int local_pid, double pT, double pphi,
					double qt, double qx, double qy, double qz, double * cosqx_dN_dypTdpTdphi, double * sinqx_dN_dypTdpTdphi)
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
	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 2.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double mT = sqrt(pT*pT+localmass*localmass);

	*cosqx_dN_dypTdpTdphi = 0.0;
	*sinqx_dN_dypTdpTdphi = 0.0;

	for (int isurf = 0; isurf < FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double xpt = surf->xpt;
		double ypt = surf->ypt;

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
			double p0 = mT*cosh(SP_p_y - eta_s[ieta]);
			double pz = mT*sinh(SP_p_y - eta_s[ieta]);

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
			{
				S_p = 0.0;
				continue;
			}

			double tpt = tau*ch_eta_s[ieta];
			double zpt = tau*sh_eta_s[ieta];
			double cosfact = 0.0;
			double sinfact = 0.0;
			for (int ii = 0; ii < 2; ++ii)
			{
				zpt *= -1.;
				double arg = tpt*qt-(xpt*qx+ypt*qy+zpt*qz);
				cosfact += cos(arg/hbarC);
				sinfact += sin(arg/hbarC);
				*cosqx_dN_dypTdpTdphi += cos(arg/hbarC)*S_p*tau*eta_s_weight[ieta];
				*sinqx_dN_dypTdpTdphi += sin(arg/hbarC)*S_p*tau*eta_s_weight[ieta];
				/*if (isnan(cosqx_dN_dypTdpTdphi) || isinf(cosqx_dN_dypTdpTdphi))
				{
					cerr << "NaN or Inf at isurf = " << isurf << " and ieta = " << ieta << endl;
					exit(1);
				}*/
			}
/*cout << "Exact: " << local_pid << "   " << isurf*eta_s_npts+ieta << "   " << isurf << "   " << ieta << "   "
		<< S_p*tau*eta_s_weight[ieta] << "   " << cosfact << "   " << sinfact << endl;*/
		}
	}

	return;
}


double CorrelationFunction::Cal_dN_dypTdpTdphi_function(int local_pid, double pT, double pphi)
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

	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 0.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double dN_dypTdpTdphi = 0.0;

	for(int isurf=0; isurf<FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double r = surf->r;
		double sin_temp_phi = surf->sin_phi;
		double cos_temp_phi = surf->cos_phi;

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

		for(int ieta=0; ieta < eta_s_npts; ++ieta)
		{
			double p0 = sqrt(pT*pT+localmass*localmass)*cosh(SP_p_y - eta_s[ieta]);
			double pz = sqrt(pT*pT+localmass*localmass)*sinh(SP_p_y - eta_s[ieta]);

			double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
	
			//viscous corrections
			double deltaf = 0.;
			if (use_delta_f)
				deltaf = deltaf_prefactor * (1. - sign*f0)
							* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

			//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
			double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

			//ignore points where delta f is large or emission function goes negative from pdsigma
			if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
				S_p = 0.0;

			dN_dypTdpTdphi += eta_even_factor*S_p*tau*eta_s_weight[ieta];
		}
	}

	return dN_dypTdpTdphi;
}

void CorrelationFunction::Determine_plane_angle()
{
	double * dN_dydphi = new double [n_pphi_pts];
	for(int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		dN_dydphi[ipphi] = 0.0;
	double ** dN_dypTdpTdphi = spectra[target_particle_id];

	/*double localmass = particle_mass;

	double * mT = new double [n_pT_pts];
	double ** px = new double * [n_pT_pts];
	double ** py = new double * [n_pT_pts];
	double ** p0 = new double * [n_pT_pts];
	double ** pz = new double * [n_pT_pts];
	

	for(int ipT=0; ipT<n_pT_pts; ++ipT)
	{
		px[ipT] = new double [n_pphi_pts];
		py[ipT] = new double [n_pphi_pts];
		p0[ipT] = new double [eta_s_npts];
		pz[ipT] = new double [eta_s_npts];
	}
   
	for(int ipT=0; ipT<n_pT_pts; ++ipT)
		mT[ipT] = sqrt(localmass*localmass + SP_pT[ipT]*SP_pT[ipT]);*/

	for(int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	{
		double cos_phi = cos(SP_pphi[ipphi]);
		double sin_phi = sin(SP_pphi[ipphi]);
		/*for(int ipT=0; ipT<n_pT_pts; ++ipT)
		{
			px[ipT][ipphi] = SP_pT[ipT]*cos_phi;
			py[ipT][ipphi] = SP_pT[ipT]*sin_phi;
		}*/
	}

	/*for(int i = 0; i < eta_s_npts; ++i)
	{
		double local_eta_s = eta_s[i];
		double local_cosh = cosh(SP_p_y - local_eta_s);
		double local_sinh = sinh(SP_p_y - local_eta_s);
		for(int ipT=0; ipT<n_pT_pts; ++ipT)
		{
			p0[ipT][i] = mT[ipT]*local_cosh;
			pz[ipT][i] = mT[ipT]*local_sinh;
		}
	}*/

	for(int ipT = 0; ipT < n_pT_pts; ++ipT)
	for(int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		dN_dydphi[ipphi] += dN_dypTdpTdphi[ipT][ipphi]*SP_pT[ipT]*SP_pT_wts[ipT];

   	double norm = 0.0e0;
   	for(int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		norm += dN_dydphi[ipphi]*SP_pphi_wts[ipphi];

   	for(int iorder = 0; iorder < n_order; iorder++)
   	{
		double cosine = 0.0e0;
		double sine = 0.0e0;
		for(int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			cosine += dN_dydphi[ipphi]*cos(iorder*SP_pphi[ipphi])*SP_pphi_wts[ipphi];
			sine += dN_dydphi[ipphi]*sin(iorder*SP_pphi[ipphi])*SP_pphi_wts[ipphi];
			//get pT-differential v_n here
			/*for(int ipT=0; ipT<n_pT_pts; ++ipT)
			{
				cosine_iorder[ipT][iorder] += dN_dypTdpTdphi[ipT][ipphi]*cos(iorder*SP_pphi[ipphi])*SP_pphi_wts[ipphi];
				sine_iorder[ipT][iorder] += dN_dypTdpTdphi[ipT][ipphi]*sin(iorder*SP_pphi[ipphi])*SP_pphi_wts[ipphi];
			}*/
		}
		/*for(int ipT=0; ipT<n_pT_pts; ++ipT)
		{
			cosine_iorder[ipT][iorder] /= dN_dypTdpT[ipT];
			sine_iorder[ipT][iorder] /= dN_dypTdpT[ipT];
		}*/
		cosine = cosine/norm;
		sine = sine/norm;
		if( sqrt(sine*sine + cosine*cosine) < 1e-8)
			plane_angle[iorder] = 0.0e0;
		else
			plane_angle[iorder] = atan2(sine, cosine)/double(iorder);
   	}
	
	//for(int ipT=0; ipT<n_pT_pts; ++ipT)
	//	dN_dypTdpT[ipT] /= (2.*M_PI);
	
	/*mean_pT = 0.;
	for(int ipphi=0; ipphi<n_pphi_pts; ++ipphi)
		mean_pT += pTdN_dydphi[ipphi]*SP_pphi_wts[ipphi];
	mean_pT /= norm;*/
	plane_angle[0] = norm;

	/*delete[] mT;
	for(int ipT=0; ipT<n_pT_pts; ++ipT)
	{
		delete [] px[ipT];
		delete [] py[ipT];
		delete [] p0[ipT];
		delete [] pz[ipT];
	}
	delete [] px;
	delete [] py;
	delete [] p0;
	delete [] pz;*/

	delete [] dN_dydphi;

	return;
}


void CorrelationFunction::R2_Fourier_transform(int iKT, double plane_psi, int mode)
{
	const int interpMode = 1;
	//int mode: 0 - GF, 1 - QM
	for(int Morder = 0; Morder < n_order; ++Morder)
	{
		double cos_mK_phi[nKphi], sin_mK_phi[nKphi];

		for(int iKphi = 0; iKphi < nKphi; ++iKphi)
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

		for(int iKphi = 0; iKphi < nKphi; ++iKphi)
		{
			double local_R2s = interpolate2D(SP_pT, SP_pphi, R2_side_GF, K_T[iKT], K_phi[iKphi], n_pT_pts, n_pphi_pts, interpMode, false, true);
			double local_R2o = interpolate2D(SP_pT, SP_pphi, R2_out_GF, K_T[iKT], K_phi[iKphi], n_pT_pts, n_pphi_pts, interpMode, false, true);
			double local_R2os = interpolate2D(SP_pT, SP_pphi, R2_outside_GF, K_T[iKT], K_phi[iKphi], n_pT_pts, n_pphi_pts, interpMode, false, true);
			double local_R2l = interpolate2D(SP_pT, SP_pphi, R2_long_GF, K_T[iKT], K_phi[iKphi], n_pT_pts, n_pphi_pts, interpMode, false, true);
			double local_R2sl = interpolate2D(SP_pT, SP_pphi, R2_sidelong_GF, K_T[iKT], K_phi[iKphi], n_pT_pts, n_pphi_pts, interpMode, false, true);
			double local_Rol = interpolate2D(SP_pT, SP_pphi, R2_outlong_GF, K_T[iKT], K_phi[iKphi], n_pT_pts, n_pphi_pts, interpMode, false, true);
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
		else if (mode == 1)
		{
			R2_side_QM_C[iKT][Morder] = temp_sum_side_cos/(2.*M_PI);
			R2_side_QM_S[iKT][Morder] = temp_sum_side_sin/(2.*M_PI);
			R2_out_QM_C[iKT][Morder] = temp_sum_out_cos/(2.*M_PI);
			R2_out_QM_S[iKT][Morder] = temp_sum_out_sin/(2.*M_PI);
			R2_outside_QM_C[iKT][Morder] = temp_sum_outside_cos/(2.*M_PI);
			R2_outside_QM_S[iKT][Morder] = temp_sum_outside_sin/(2.*M_PI);
			R2_long_QM_C[iKT][Morder] = temp_sum_long_cos/(2.*M_PI);
			R2_long_QM_S[iKT][Morder] = temp_sum_long_sin/(2.*M_PI);
			R2_sidelong_QM_C[iKT][Morder] = temp_sum_sidelong_cos/(2.*M_PI);
			R2_sidelong_QM_S[iKT][Morder] = temp_sum_sidelong_sin/(2.*M_PI);
			R2_outlong_QM_C[iKT][Morder] = temp_sum_outlong_cos/(2.*M_PI);
			R2_outlong_QM_S[iKT][Morder] = temp_sum_outlong_sin/(2.*M_PI);
		}
	}

	return;
}

//End of file
