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

const int n_refinement_pts = 101;
double Delta_DpY;
const double PTCHANGE = 1.0;
gsl_cheb_series *cs_accel_expEdNd3p;

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

void CorrelationFunction::Tabulate_resonance_Chebyshev_coefficients(int parent_resonance_particle_id)
{
	cs_accel_expEdNd3p = gsl_cheb_alloc (n_pY_pts - 1);
	cs_accel_expEdNd3p->a = SP_Del_pY_min;
	cs_accel_expEdNd3p->b = SP_Del_pY_max;

	long cfs_array_length = n_pT_pts * n_pphi_pts * qxnpts * qynpts * ntrig;
	chebyshev_a_cfs = new double * [cfs_array_length];
	for (int icf = 0; icf < cfs_array_length; ++icf)
		chebyshev_a_cfs[icf] = new double [n_pY_pts];

	int idx = 0;
	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int itrig = 0; itrig < ntrig; ++itrig)
	{
		for (int ipY = 0; ipY < n_pY_pts; ++ipY)
		{
			chebyshev_a_cfs[idx][ipY] = 0.0;
			for (int kpY = 0; kpY < n_pY_pts; ++kpY)
				chebyshev_a_cfs[idx][ipY] += exp(abs(SP_Del_pY[kpY])) * chebTcfs[ipY * n_pY_pts + kpY] 
												* current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,kpY,iqx,iqy,itrig)];
		}
		++idx;
	}

	return;
}

void CorrelationFunction::Refine_resonance_grids(int parent_resonance_particle_id)
{
	Delta_DpY = (SP_Del_pY_max - SP_Del_pY_min) / (double)(n_refinement_pts - 1);

	long cfs_array_length = n_pT_pts * n_pphi_pts * qxnpts * qynpts * ntrig * n_refinement_pts;
	refined_resonance_grids = new double [cfs_array_length];

	int idx = 0;
	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int itrig = 0; itrig < ntrig; ++itrig)
	{
		cs_accel_expEdNd3p->c = chebyshev_a_cfs[idx];
		for (int iii = 0; iii < n_refinement_pts; ++iii)
		{
			double tmp_pY = SP_Del_pY_min + (double)iii * Delta_DpY;
			refined_resonance_grids[idx*n_refinement_pts+iii] = exp(-tmp_pY) * gsl_cheb_eval (cs_accel_expEdNd3p, tmp_pY);
		}
		++idx;
	}

	return;
}

void CorrelationFunction::Do_resonance_integrals(int parent_resonance_particle_id, int daughter_particle_id, int decay_channel, int iqt, int iqz)
{
	time_t rawtime;
  	struct tm * timeinfo;

	int daughter_lookup_idx = distance(daughter_resonance_indices.begin(), daughter_resonance_indices.find(daughter_particle_id));

	Allocate_resonance_running_sum_vectors();

	//set these for quick look-up in EdNd3p() routine
	spec_vals_info = spectra[parent_resonance_particle_id];
	spec_log_info = log_spectra[parent_resonance_particle_id];
	spec_sign_info = sign_spectra[parent_resonance_particle_id];

	Tabulate_resonance_Chebyshev_coefficients(parent_resonance_particle_id);
	Refine_resonance_grids(parent_resonance_particle_id);

	int tmp_parent_monval = all_particles[parent_resonance_particle_id].monval;
	int tmp_daughter_monval = all_particles[daughter_particle_id].monval;
	n_body = current_reso_nbody;
	current_parent_resonance = parent_resonance_particle_id;

	local_verbose = 0;

	double loc_qz = qz_pts[current_iqz];
	double loc_qt = qt_pts[current_iqt];
	current_pY_shift = - double(abs(loc_qz)>1.e-10) * asinh(loc_qz / sqrt(abs(loc_qt*loc_qt-loc_qz*loc_qz) + 1.e-100));

	if (n_body == 2)
	{
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		for (int ipY = 0; ipY < n_pY_pts; ++ipY)
		{
			double local_pT = SP_pT[ipT];
			double local_pphi = SP_pphi[ipphi];
			double local_pY = SP_Del_pY[ipY] + current_pY_shift;
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
					double PKT = VEC_n2_PT[NB2_indexer(iv,izeta)];
					double PKY = VEC_n2_P_Y[iv] - current_pY_shift;
					double PKphi = VEC_n2_PPhi_tilde[NB2_indexer(iv,izeta)];
					for (int tempidx = 0; tempidx <= 1; ++tempidx)
					{
						if (tempidx != 0)
							PKphi = VEC_n2_PPhi_tildeFLIP[NB2_indexer(iv,izeta)];		//also takes Pp --> Pm
						currentPpm = VEC_n2_Ppm[NB2_indexer(iv,izeta)*2 + tempidx];
						//spectra
						if ( iqt == 0 && iqz == 0 )
							Edndp3(PKT, PKphi, &Csum);
						//space-time moments
						if (!IGNORE_LONG_LIVED_RESONANCES || Gamma >= hbarC / max_lifetime)
							eiqxEdndp3(PKT, PKphi, abs(PKY), Csum_vec, local_verbose);
					}												// end of tempidx sum
					for (int qpt_cs_idx = 0; qpt_cs_idx < qspace_cs_slice_length; ++qpt_cs_idx)
						zetasum_vec[qpt_cs_idx] += VEC_n2_zeta_factor[NB2_indexer(iv,izeta)]*Csum_vec[qpt_cs_idx];
					zetasum += VEC_n2_zeta_factor[NB2_indexer(iv,izeta)]*Csum;
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
			for (int iqx = 0; iqx < qxnpts; ++iqx)
			for (int iqy = 0; iqy < qynpts; ++iqy)
			for (int itrig = 0; itrig < 2; ++itrig)
			{
				current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,itrig)] += ssum_vec[qpt_cs_idx];
				++qpt_cs_idx;
			}

			//update daughter spectra separately
			spectra[daughter_particle_id][ipT][ipphi] += ssum;
			log_spectra[daughter_particle_id][ipT][ipphi] = log(abs(spectra[daughter_particle_id][ipT][ipphi])+1.e-100);
			sign_spectra[daughter_particle_id][ipT][ipphi] = sgn(spectra[daughter_particle_id][ipT][ipphi]);
		}											// end of pT, pphi, pY loops
	}												// end of nbody == 2
	else
	{
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		for (int ipY = 0; ipY < n_pY_pts; ++ipY)
		{
			double local_pT = SP_pT[ipT];
			double local_pphi = SP_pphi[ipphi];
			double local_pY = SP_Del_pY[ipY] + current_pY_shift;
			current_ipT = ipT;
			current_ipphi = ipphi;
			current_ipphi = ipY;
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
						double PKT = VEC_n3_PT[NB3_indexer(is,iv,izeta)];
						double PKY = VEC_n3_P_Y[is*n_v_pts+iv] - current_pY_shift;
						double PKphi = VEC_n3_PPhi_tilde[NB3_indexer(is,iv,izeta)];
						for (int tempidx = 0; tempidx <= 1; ++tempidx)
						{
							if (tempidx != 0)
								PKphi = VEC_n3_PPhi_tildeFLIP[NB3_indexer(is,iv,izeta)];		//also takes Pp --> Pm
							currentPpm = VEC_n3_Ppm[NB3_indexer(is,iv,izeta)*2+tempidx];
							//spectra
							if ( iqt == 0 && iqz == 0 )
								Edndp3(PKT, PKphi, &Csum);
							//space-time moments
							if (!IGNORE_LONG_LIVED_RESONANCES || Gamma >= hbarC / max_lifetime)
								eiqxEdndp3(PKT, PKphi, abs(PKY), Csum_vec, local_verbose);
						}										// end of tempidx sum
						for (int qpt_cs_idx = 0; qpt_cs_idx < qspace_cs_slice_length; ++qpt_cs_idx)
							zetasum_vec[qpt_cs_idx] += VEC_n3_zeta_factor[NB3_indexer(is,iv,izeta)]*Csum_vec[qpt_cs_idx];
						zetasum += VEC_n3_zeta_factor[NB3_indexer(is,iv,izeta)]*Csum;
					}											// end of zeta sum
					for (int qpt_cs_idx = 0; qpt_cs_idx < qspace_cs_slice_length; ++qpt_cs_idx)
					    vsum_vec[qpt_cs_idx] += VEC_n3_v_factor[is*n_v_pts+iv]*zetasum_vec[qpt_cs_idx];
					vsum += VEC_n3_v_factor[is*n_v_pts+iv]*zetasum;
				}												// end of v sum
				for (int qpt_cs_idx = 0; qpt_cs_idx < qspace_cs_slice_length; ++qpt_cs_idx)
					ssum_vec[qpt_cs_idx] += Mres*VEC_n3_s_factor[is]*vsum_vec[qpt_cs_idx];
				ssum += Mres*VEC_n3_s_factor[is]*vsum;
			}													// end of s sum
			//update all gridpoints for daughter moments
			int qpt_cs_idx = 0;
			for (int iqx = 0; iqx < qxnpts; ++iqx)
			for (int iqy = 0; iqy < qynpts; ++iqy)
			for (int itrig = 0; itrig < 2; ++itrig)
			{
				current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,itrig)] += ssum_vec[qpt_cs_idx];
				++qpt_cs_idx;
			}

			//update daughter spectra separately
			spectra[daughter_particle_id][ipT][ipphi] += ssum;
			log_spectra[daughter_particle_id][ipT][ipphi] = log(abs(spectra[daughter_particle_id][ipT][ipphi])+1.e-100);
			sign_spectra[daughter_particle_id][ipT][ipphi] = sgn(spectra[daughter_particle_id][ipT][ipphi]);
		}								// end of pT, pphi, pY loops
	}										// end of nbody == 3

	// clean up

	long cfs_array_length = n_pT_pts * n_pphi_pts * qxnpts * qynpts * ntrig;
	for (int icf = 0; icf < cfs_array_length; ++icf)
		delete [] chebyshev_a_cfs[icf];
	delete [] chebyshev_a_cfs;
	delete [] refined_resonance_grids;

	Delete_resonance_running_sum_vectors();

	return;
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
void CorrelationFunction::Edndp3(double ptr, double pphir, double * result, int loc_verb /*==0*/)
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
	if(pphir < SP_pphi[0])			//if angle is less than minimum angle grid point
	{
		phi0 = SP_pphi[npphi_max] - 2. * M_PI;
		phi1 = SP_pphi[0];
		nphi = 0;
		nphim1 = npphi_max;
	}
	else if(pphir > SP_pphi[npphi_max])	//if angle is greater than maximum angle grid point
	{
		phi0 = SP_pphi[npphi_max];
		phi1 = SP_pphi[0] + 2. * M_PI;
		nphi = 0;
		nphim1 = npphi_max;
	}
	else						//if angle is within grid range
	{
		while ((pphir > SP_pphi[nphi]) &&
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
	*result += lin_int(pphir-phi0, one_by_pphidiff, f1, f2);

	return;
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
void CorrelationFunction::eiqxEdndp3(double ptr, double phir, double pyr, double * results, int loc_verb /*==0*/)
{
	double phi0, phi1, py0, py1;
	double f1, f2;

	int npphi_max = n_pphi_pts - 1;
	int npT_max = n_pT_pts - 1;
	int npY_max = n_pY_pts - 1;

	// locate pT interval
	int npt = 1;
	while ((ptr > SP_pT[npt]) &&
		(npt < npT_max)) ++npt;
	double pT0 = SP_pT[npt-1];
	double pT1 = SP_pT[npt];

	// locate pphi interval
	int nphi = 1, nphim1 = 0;
	if(phir < SP_pphi[0])                   //if angle is less than minimum angle grid point
	{
		phi0 = SP_pphi[npphi_max] - 2. * M_PI;
		phi1 = SP_pphi[0];
		nphi = 0;
		nphim1 = npphi_max;
	}
	else if(phir > SP_pphi[npphi_max])      //if angle is greater than maximum angle grid point
	{
		phi0 = SP_pphi[npphi_max];
		phi1 = SP_pphi[0] + 2. * M_PI;
		nphi = 0;
		nphim1 = npphi_max;
	}
	else                                            //if angle is within grid range
	{
		while ((phir > SP_pphi[nphi]) &&
			(nphi < npphi_max)) ++nphi;
		nphim1 = nphi - 1;
		phi0 = SP_pphi[nphim1];
		phi1 = SP_pphi[nphi];
	}

	// locate py interval
	long npy = 1, npym1 = 0;
	/*if(pyr < SP_Del_pY[0])			//if rapidity is less than minimum rapidity grid point
	{
		py0 = SP_Del_pY[0];
		py1 = SP_Del_pY[1];
		npy = 1;
		npym1 = 0;
	}
	else if(pyr > SP_Del_pY[npY_max])	//if rapidity is greater than maximum rapidity grid point
	{
		py0 = SP_Del_pY[npY_max];
		py1 = SP_Del_pY[npY_max];
		npy = npY_max;	//this just guarantees nearest neighbor if pY > pY_max
		npym1 = npY_max;
	}
	else						//if rapidity is within grid range
	{
		while ((pyr > SP_Del_pY[npy]) &&
				(npy < npY_max)) ++npy;
		npym1 = npy - 1;
		py0 = SP_Del_pY[npym1];
		py1 = SP_Del_pY[npy];
	}*/
	if(pyr > SP_Del_pY_max)	//if rapidity is greater than maximum rapidity grid point
	{
		py0 = SP_Del_pY_max;
		py1 = SP_Del_pY_max;
		npy = n_refinement_pts - 1;	//this just guarantees nearest neighbor if pY > pY_max
		npym1 = n_refinement_pts - 1;
	}
	else						//if rapidity is within grid range
	{
		//while ((pyr > SP_Del_pY[npy]) &&
		//		(npy < npY_max)) ++npy;
		long npym1 = floor((pyr-SP_Del_pY_min)/Delta_DpY);
		long npy = npym1 + 1;
		py0 = SP_Del_pY_min + (double)npym1 * Delta_DpY;
		py1 = SP_Del_pY_min + (double)npy * Delta_DpY;
	}

	if (pT0==pT1 || phi0==phi1)
	{
		cerr << "ERROR in eiqxEdndp3(): pT and/or pphi values equal!" << endl;
		exit(1);
	}

	double one_by_pTdiff = 1./(pT1 - pT0), one_by_pphidiff = 1./(phi1 - phi0), one_by_pYdiff = 1./(py1 - py0 + 1.e-100);
	double del_ptr_pt0 = ptr - pT0, del_phir_phi0 = phir - phi0, del_pyr_py0 = pyr - py0;

	// interpolate over pY first, then proceed as usual
	// set index for looping
	int qpt_cs_idx = 0;
	double f11_arr[qxnpts*qynpts*ntrig], f12_arr[qxnpts*qynpts*ntrig], f21_arr[qxnpts*qynpts*ntrig], f22_arr[qxnpts*qynpts*ntrig];
	double log_f11_arr[qxnpts*qynpts*ntrig], log_f12_arr[qxnpts*qynpts*ntrig], log_f21_arr[qxnpts*qynpts*ntrig], log_f22_arr[qxnpts*qynpts*ntrig];
	double sign_of_f11_arr[qxnpts*qynpts*ntrig], sign_of_f12_arr[qxnpts*qynpts*ntrig], sign_of_f21_arr[qxnpts*qynpts*ntrig], sign_of_f22_arr[qxnpts*qynpts*ntrig];

	// choose pt-pphi slice of resonance info arrays
	// flatten arrays, since these are quicker to process
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int itrig = 0; itrig < 2; ++itrig)
	{
		/*f11_arr[qpt_cs_idx] = Get_EdNd3p_moment(((npt-1)*n_pphi_pts+nphim1)*qxnpts + qpt_cs_idx, pyr);
		f12_arr[qpt_cs_idx] = Get_EdNd3p_moment(((npt-1)*n_pphi_pts+nphi)*qxnpts + qpt_cs_idx, pyr);
		f21_arr[qpt_cs_idx] = Get_EdNd3p_moment((npt*n_pphi_pts+nphim1)*qxnpts + qpt_cs_idx, pyr);
		f22_arr[qpt_cs_idx] = Get_EdNd3p_moment((npt*n_pphi_pts+nphi)*qxnpts + qpt_cs_idx, pyr);*/
		//cout << (((npt-1)*n_pphi_pts+nphim1)*qxnpts + qpt_cs_idx)*n_refinement_pts+npym1 << "   " << (((npt-1)*n_pphi_pts+nphim1)*qxnpts + qpt_cs_idx)*n_refinement_pts+npy << endl;
		f11_arr[qpt_cs_idx] = lin_int( del_pyr_py0, one_by_pYdiff, refined_resonance_grids[(((npt-1)*n_pphi_pts+nphim1)*qxnpts + qpt_cs_idx)*n_refinement_pts+npym1],
																	refined_resonance_grids[(((npt-1)*n_pphi_pts+nphim1)*qxnpts + qpt_cs_idx)*n_refinement_pts+npy]);
		f12_arr[qpt_cs_idx] = lin_int( del_pyr_py0, one_by_pYdiff, refined_resonance_grids[(((npt-1)*n_pphi_pts+nphi)*qxnpts + qpt_cs_idx)*n_refinement_pts+npym1],
																	refined_resonance_grids[(((npt-1)*n_pphi_pts+nphi)*qxnpts + qpt_cs_idx)*n_refinement_pts+npy]);
		f21_arr[qpt_cs_idx] = lin_int( del_pyr_py0, one_by_pYdiff, refined_resonance_grids[((npt*n_pphi_pts+nphim1)*qxnpts + qpt_cs_idx)*n_refinement_pts+npym1],
																	refined_resonance_grids[((npt*n_pphi_pts+nphim1)*qxnpts + qpt_cs_idx)*n_refinement_pts+npy]);
		f22_arr[qpt_cs_idx] = lin_int( del_pyr_py0, one_by_pYdiff, refined_resonance_grids[((npt*n_pphi_pts+nphi)*qxnpts + qpt_cs_idx)*n_refinement_pts+npym1],
																	refined_resonance_grids[((npt*n_pphi_pts+nphi)*qxnpts + qpt_cs_idx)*n_refinement_pts+npy]);

		log_f11_arr[qpt_cs_idx] = log(abs(f11_arr[qpt_cs_idx]) + 1.e-100);
		log_f12_arr[qpt_cs_idx] = log(abs(f12_arr[qpt_cs_idx]) + 1.e-100);
		log_f21_arr[qpt_cs_idx] = log(abs(f21_arr[qpt_cs_idx]) + 1.e-100);
		log_f22_arr[qpt_cs_idx] = log(abs(f22_arr[qpt_cs_idx]) + 1.e-100);

		sign_of_f11_arr[qpt_cs_idx] = sgn(f11_arr[qpt_cs_idx]);
		sign_of_f12_arr[qpt_cs_idx] = sgn(f12_arr[qpt_cs_idx]);
		sign_of_f21_arr[qpt_cs_idx] = sgn(f21_arr[qpt_cs_idx]);
		sign_of_f22_arr[qpt_cs_idx] = sgn(f22_arr[qpt_cs_idx]);

		++qpt_cs_idx;
	}


	// set index for looping
	qpt_cs_idx = 0;
	int qlist_idx = 0;

	if (ptr > PTCHANGE)                             // if pT interpolation point is larger than PTCHANGE (currently 1.0 GeV)
	{
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
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
			else if (sign_of_f11 * sign_of_f21 > 0) // if the two points have the same sign in the pT direction, interpolate logs
				f1 = sign_of_f11 * exp( lin_int(del_ptr_pt0, one_by_pTdiff, log_f11_arr[qpt_cs_idx], log_f21_arr[qpt_cs_idx]) );
			else                                    // otherwise, just interpolate original vals
				f1 = lin_int(del_ptr_pt0, one_by_pTdiff, f11_arr[qpt_cs_idx], f21_arr[qpt_cs_idx]);

			//*******************************************************************************************************************
			// set f2 next
			//*******************************************************************************************************************
			if ( ptr > pT1 && ( log_f22_arr[qpt_cs_idx] > log_f12_arr[qpt_cs_idx] || sign_of_f12 * sign_of_f22 < 0 ) )
					f2 = 0.0;
			else if (sign_of_f12 * sign_of_f22 > 0) // if the two points have the same sign in the pT direction, interpolate logs
					f2 = sign_of_f12 * exp( lin_int(del_ptr_pt0, one_by_pTdiff, log_f12_arr[qpt_cs_idx], log_f22_arr[qpt_cs_idx]) );
			else                                    // otherwise, just interpolate original vals
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
	        else if (sign_of_f11 * sign_of_f21 > 0) // if the two points have the same sign in the pT direction, interpolate logs
				f1 = sign_of_f11 * exp( lin_int(del_ptr_pt0, one_by_pTdiff, log_f11_arr[qpt_cs_idx+1], log_f21_arr[qpt_cs_idx+1]) );
	        else                                    // otherwise, just interpolate original vals
				f1 = lin_int(del_ptr_pt0, one_by_pTdiff, f11_arr[qpt_cs_idx+1], f21_arr[qpt_cs_idx+1]);

	        //*******************************************************************************************************************
	        // set f2 next
	        //*******************************************************************************************************************
	        if ( ptr > pT1 && ( log_f22_arr[qpt_cs_idx+1] > log_f12_arr[qpt_cs_idx+1] || sign_of_f12 * sign_of_f22 < 0 ) )
				f2 = 0.0;
	        else if (sign_of_f12 * sign_of_f22 > 0) // if the two points have the same sign in the pT direction, interpolate logs
				f2 = sign_of_f12 * exp( lin_int(del_ptr_pt0, one_by_pTdiff, log_f12_arr[qpt_cs_idx+1], log_f22_arr[qpt_cs_idx+1]) );
	        else                                    // otherwise, just interpolate original vals
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
		}   //end of all q-loops
	}
	else                                            // if pT is smaller than PTCHANGE, just use linear interpolation, no matter what
	{
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
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
		}       //end of all q-loops
	}

	return;
}

//End of file


///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
void CorrelationFunction::Load_decay_channel_info(int dc_idx, double K_T_local, double K_phi_local, double K_y_local)
{
	Mres = current_resonance_mass;
	Gamma = current_resonance_Gamma;
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
		double pstar_loc = sqrt( ((Mres+mass)*(Mres+mass) - s_loc)*((Mres-mass)*(Mres-mass) - s_loc) )/(2.0*Mres);
		double g_s_loc = g(s_loc);	//for n_body == 2, doesn't actually use s_loc since result is just a factor * delta(...); just returns factor
		double Estar_loc = sqrt(mass*mass + pstar_loc*pstar_loc);
		double psBmT = pstar_loc / mT;
		double DeltaY_loc = log(psBmT + sqrt(1.+psBmT*psBmT));

		VEC_n2_s_factor = br/(4.*M_PI*pstar_loc);	//==g_s_loc

		for(int iv = 0; iv < n_v_pts; ++iv)
		{
			double v_loc = v_pts[iv];
			double P_Y_loc = p_y + v_loc*DeltaY_loc;
			double mT_ch_P_Y_p_y = mT*cosh(v_loc*DeltaY_loc);
			double x2 = mT_ch_P_Y_p_y*mT_ch_P_Y_p_y - pT*pT;
			double MTbar_loc = Estar_loc*Mres*mT_ch_P_Y_p_y/x2;
			double DeltaMT_loc = Mres*pT*sqrt(Estar_loc*Estar_loc - x2)/x2;

			VEC_n2_P_Y[iv] = P_Y_loc;
			VEC_n2_v_factor[iv] = v_wts[iv]*DeltaY_loc/sqrt(x2);

			for(int izeta = 0; izeta < n_zeta_pts; ++izeta)
			{
				double zeta_loc = zeta_pts[izeta];
				double MT_loc = MTbar_loc + cos(zeta_loc)*DeltaMT_loc;
				double PT_loc = sqrt(MT_loc*MT_loc - Mres*Mres);
				double temp_cos_PPhi_tilde_loc = (mT*MT_loc*cosh(P_Y_loc-p_y) - Estar_loc*Mres)/(pT*PT_loc);
				//assume that PPhi_tilde is +ve in next step...
				double temp_sin_PPhi_tilde_loc = sqrt(1. - temp_cos_PPhi_tilde_loc*temp_cos_PPhi_tilde_loc);
				double PPhi_tilde_loc = place_in_range( atan2(temp_sin_PPhi_tilde_loc, temp_cos_PPhi_tilde_loc), Kphi_min, Kphi_max);

				VEC_n2_zeta_factor[NB2_indexer(iv, izeta)] = zeta_wts[izeta]*MT_loc;
				VEC_n2_PPhi_tilde[NB2_indexer(iv, izeta)] = place_in_range( K_phi_local + PPhi_tilde_loc, Kphi_min, Kphi_max);
				VEC_n2_PPhi_tildeFLIP[NB2_indexer(iv, izeta)] = place_in_range( K_phi_local - PPhi_tilde_loc, Kphi_min, Kphi_max);
				VEC_n2_PT[NB2_indexer(iv, izeta)] = PT_loc;
				//set P^+ components
				VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+0][0] = MT_loc * cosh(P_Y_loc);
				VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+0][1] = PT_loc * cos(K_phi_local + PPhi_tilde_loc);
				VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+0][2] = PT_loc * sin(K_phi_local + PPhi_tilde_loc);
				VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+0][3] = MT_loc * sinh(P_Y_loc);
				//set P^- components
				VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+1][0] = MT_loc * cosh(P_Y_loc);
				VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+1][1] = PT_loc * cos(K_phi_local - PPhi_tilde_loc);
				VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+1][2] = PT_loc * sin(K_phi_local - PPhi_tilde_loc);
				VEC_n2_Ppm[NB2_indexer(iv, izeta)*2+1][3] = MT_loc * sinh(P_Y_loc);
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

		// s-loop
		for (int is = 0; is < n_s_pts; ++is)
		{

			double s_loc = s_pts[is];
			double g_s_loc = g(s_loc);
			double pstar_loc = sqrt(((Mres+mass)*(Mres+mass) - s_loc)*((Mres-mass)*(Mres-mass) - s_loc))/(2.0*Mres);
			double Estar_loc = sqrt(mass*mass + pstar_loc*pstar_loc);
			double psBmT = pstar_loc / mT;
			double DeltaY_loc = log(psBmT + sqrt(1.+psBmT*psBmT));

			VEC_n3_s_factor[is] = s_wts[is]*g_s_loc;

			// v-loop
			for(int iv = 0; iv < n_v_pts; ++iv)
			{
				double v_loc = v_pts[iv];
				double P_Y_loc = p_y + v_loc*DeltaY_loc;
				double mT_ch_P_Y_p_y = mT*cosh(v_loc*DeltaY_loc);
				double x2 = mT_ch_P_Y_p_y*mT_ch_P_Y_p_y - pT*pT;
				double MTbar_loc = Estar_loc*Mres*mT_ch_P_Y_p_y/x2;
				double DeltaMT_loc = Mres*pT*sqrt(Estar_loc*Estar_loc - x2)/x2;

				VEC_n3_P_Y[is * n_v_pts + iv] = P_Y_loc;
				VEC_n3_v_factor[is * n_v_pts + iv] = v_wts[iv]*DeltaY_loc/sqrt(x2);

				// zeta-loop
				for(int izeta = 0; izeta < n_zeta_pts; ++izeta)
				{
					double zeta_loc = zeta_pts[izeta];
					double MT_loc = MTbar_loc + cos(zeta_loc)*DeltaMT_loc;
					double PT_loc = sqrt(MT_loc*MT_loc - Mres*Mres);
					double temp_cos_PPhi_tilde_loc = (mT*MT_loc*cosh(P_Y_loc-p_y) - Estar_loc*Mres)/(pT*PT_loc);
					//assume that PPhi_tilde is +ve in next step...
					double temp_sin_PPhi_tilde_loc = sqrt(1. - temp_cos_PPhi_tilde_loc*temp_cos_PPhi_tilde_loc);
					double PPhi_tilde_loc = place_in_range( atan2(temp_sin_PPhi_tilde_loc, temp_cos_PPhi_tilde_loc), Kphi_min, Kphi_max);

					VEC_n3_zeta_factor[NB3_indexer(is, iv, izeta)] = zeta_wts[izeta]*MT_loc;
					VEC_n3_PPhi_tilde[NB3_indexer(is, iv, izeta)] = place_in_range( K_phi_local + PPhi_tilde_loc, Kphi_min, Kphi_max);
					VEC_n3_PPhi_tildeFLIP[NB3_indexer(is, iv, izeta)] = place_in_range( K_phi_local - PPhi_tilde_loc, Kphi_min, Kphi_max);
					VEC_n3_PT[NB3_indexer(is, iv, izeta)] = PT_loc;
					//set P^+ components
					VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+0][0] = MT_loc * cosh(P_Y_loc);
					VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+0][1] = PT_loc * cos(K_phi_local + PPhi_tilde_loc);
					VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+0][2] = PT_loc * sin(K_phi_local + PPhi_tilde_loc);
					VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+0][3] = MT_loc * sinh(P_Y_loc);
					//set P^- components
					VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+1][0] = MT_loc * cosh(P_Y_loc);
					VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+1][1] = PT_loc * cos(K_phi_local - PPhi_tilde_loc);
					VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+1][2] = PT_loc * sin(K_phi_local - PPhi_tilde_loc);
					VEC_n3_Ppm[NB3_indexer(is, iv, izeta)*2+1][3] = MT_loc * sinh(P_Y_loc);
				}
			}
		}
	}

	return;
}

//End of file
