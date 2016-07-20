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
		for (int ipt = 0; ipt < n_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			double local_pT = SP_pT[ipt];
			double local_pphi = SP_pphi[ipphi];
			current_ipt = ipt;
			current_ipphi = ipphi;
			//current_qlist_slice = qlist[ipt];
			Zero_resonance_running_sum_vector(ssum_vec);
			Zero_resonance_running_sum_vector(vsum_vec);
			Zero_resonance_running_sum_vector(zetasum_vec);
			Zero_resonance_running_sum_vector(Csum_vec);
			Load_decay_channel_info(decay_channel, local_pT, local_pphi);	// set decay channel information

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
				current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipt,ipphi,iqt,iqx,iqy,iqz,itrig)] += ssum_vec[qpt_cs_idx];
				++qpt_cs_idx;
			}

			//update daughter spectra separately
			spectra[daughter_particle_id][ipt][ipphi] += ssum;
			log_spectra[daughter_particle_id][ipt][ipphi] = log(abs(spectra[daughter_particle_id][ipt][ipphi])+1.e-100);
			sign_spectra[daughter_particle_id][ipt][ipphi] = sgn(spectra[daughter_particle_id][ipt][ipphi]);

			//only do this if all qpoints array sizes are odd!
			//now, if ignoring long-lived resonances, take them out of the correlator numerator, but keep them in the denominator (AKA, spectra)
			if (IGNORE_LONG_LIVED_RESONANCES && qtnpts%2==1 && qxnpts%2==1 && qynpts%2==1 && qznpts%2==1 && Gamma < hbarC / max_lifetime)
				current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipt,ipphi,iqt0,iqx0,iqy0,iqz0,0)] -= ssum_vec[tmp_qpt_cs_idx];
	
			if (isnan(current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipt,ipphi,0,0,0,0,0)]
					+ current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipt,ipphi,0,0,0,0,1)]))
			{
				*global_out_stream_ptr << "ERROR: NaNs encountered!" << endl
										<< "current_daughters_dN_dypTdpTdphi_moments[" << daughter_lookup_idx << "][" << ipt << "][" << ipphi << "][0][0][0][0][0] = "
										<< setw(8) << setprecision(15)
										<< current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipt,ipphi,0,0,0,0,0)] << endl
										<< "current_daughters_dN_dypTdpTdphi_moments[" << daughter_lookup_idx << "][" << ipt << "][" << ipphi << "][0][0][0][0][1] = "
										<< setw(8) << setprecision(15)
										<< current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipt,ipphi,0,0,0,0,1)] << endl
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
		for (int ipt = 0; ipt < n_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			double local_pT = SP_pT[ipt];
			double local_pphi = SP_pphi[ipphi];
			current_ipt = ipt;
			current_ipphi = ipphi;
			//current_qlist_slice = qlist[ipt];
			Zero_resonance_running_sum_vector(ssum_vec);
			Zero_resonance_running_sum_vector(vsum_vec);
			Zero_resonance_running_sum_vector(zetasum_vec);
			Zero_resonance_running_sum_vector(Csum_vec);
			Load_decay_channel_info(decay_channel, local_pT, local_pphi);	// set decay channel information

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
				current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipt,ipphi,iqt,iqx,iqy,iqz,itrig)] += ssum_vec[qpt_cs_idx];
				++qpt_cs_idx;
			}

			//update daughter spectra separately
			spectra[daughter_particle_id][ipt][ipphi] += ssum;
			log_spectra[daughter_particle_id][ipt][ipphi] = log(abs(spectra[daughter_particle_id][ipt][ipphi])+1.e-100);
			sign_spectra[daughter_particle_id][ipt][ipphi] = sgn(spectra[daughter_particle_id][ipt][ipphi]);

			//only do this if all qpoints array sizes are odd!
			//now, if ignoring long-lived resonances, take them out of the correlator numerator, but keep them in the denominator (AKA, spectra)
			if (IGNORE_LONG_LIVED_RESONANCES && qtnpts%2==1 && qxnpts%2==1 && qynpts%2==1 && qznpts%2==1 && Gamma < hbarC / max_lifetime)
				current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipt,ipphi,iqt0,iqx0,iqy0,iqz0,0)] -= ssum_vec[tmp_qpt_cs_idx];

			if (isnan(current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipt,ipphi,0,0,0,0,0)]
					+ current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipt,ipphi,0,0,0,0,1)]))
			{
				*global_out_stream_ptr << "ERROR: NaNs encountered!" << endl
										<< "current_daughters_dN_dypTdpTdphi_moments[" << daughter_lookup_idx << "][" << ipt << "][" << ipphi << "][0][0][0][0][0] = "
										<< setw(8) << setprecision(15)
										<< current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipt,ipphi,0,0,0,0,0)] << endl
										<< "current_daughters_dN_dypTdpTdphi_moments[" << daughter_lookup_idx << "][" << ipt << "][" << ipphi << "][0][0][0][0][1] = "
										<< setw(8) << setprecision(15)
										<< current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][indexer(ipt,ipphi,0,0,0,0,1)] << endl
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
	double lls[dim_loc] = { interp_pT_min, interp_pphi_min };
	double uls[dim_loc] = { interp_pT_max, interp_pphi_max };
	int modes_loc[dim_loc] = { 0, 0 };
	//double lls[dim_loc] = { 0.0, interp_pphi_min };
	//double uls[dim_loc] = { (1.-sin(M_PI/n_pT_pts))/(1.+sin(M_PI/n_pT_pts)), interp_pphi_max };
	//int modes_loc[dim_loc] = { 1, 0 };

	int momidx = 0;

	for (int ipt = 0; ipt < n_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
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
				//res_sign_info[ipt][ipphi][qpt_cs_idx] = current_sign_of_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][itrig];
				//res_log_info[ipt][ipphi][qpt_cs_idx] = current_ln_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][itrig];
				//res_moments_info[ipt][ipphi][qpt_cs_idx] = current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][itrig];
				res_sign_info[ipt][ipphi][qpt_cs_idx] = current_sign_of_dN_dypTdpTdphi_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,itrig)];
				res_log_info[ipt][ipphi][qpt_cs_idx] = current_ln_dN_dypTdpTdphi_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,itrig)];
				res_moments_info[ipt][ipphi][qpt_cs_idx] = current_dN_dypTdpTdphi_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,itrig)];
				++qpt_cs_idx;
			}

			//tmp_moments_real[iqt][iqx][iqy][iqz][momidx] = current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0];
			//tmp_moments_imag[iqt][iqx][iqy][iqz][momidx] = current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1];
			tmp_moments_real[iqt][iqx][iqy][iqz][momidx] = current_dN_dypTdpTdphi_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,0)];
			tmp_moments_imag[iqt][iqx][iqy][iqz][momidx] = current_dN_dypTdpTdphi_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,1)];
		}

		flat_spectra[momidx] = spectra[parent_resonance_particle_id][ipt][ipphi];

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

//End of file
