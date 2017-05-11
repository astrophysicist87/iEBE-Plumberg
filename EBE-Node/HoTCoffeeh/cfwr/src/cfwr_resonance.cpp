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

double * val11_arr, * val12_arr, * val21_arr, * val22_arr;

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
				chebyshev_a_cfs[idx][ipY] += exp(SP_Del_pY[kpY]) * chebTcfs[ipY * n_pY_pts + kpY] * current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,kpY,iqx,iqy,itrig)];
		}
		++idx;
	}

	return;
}

void CorrelationFunction::Refine_resonance_grids(int parent_resonance_particle_id)
{
	//Delta_DpY = (SP_Del_pY_max - SP_Del_pY_min) / (double)(n_refinement_pts - 1);
	Delta_DpY = (SP_Del_pY_max - 0.0) / (double)(n_refinement_pts - 1);

	long cfs_array_length = qxnpts * qynpts * ntrig;
	refined_resonance_grids = new double * [n_pT_pts * n_pphi_pts * n_refinement_pts];
	log_refined_grids = new double * [n_pT_pts * n_pphi_pts * n_refinement_pts];
	sgn_refined_grids = new double * [n_pT_pts * n_pphi_pts * n_refinement_pts];
	for (int ipTpphipY = 0; ipTpphipY < n_pT_pts * n_pphi_pts * n_refinement_pts; ++ipTpphipY)
	{
		refined_resonance_grids[ipTpphipY] = new double [cfs_array_length];
		log_refined_grids[ipTpphipY] = new double [cfs_array_length];
		sgn_refined_grids[ipTpphipY] = new double [cfs_array_length];
	}

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
			int tmp_index = (ipT * n_pphi_pts + ipphi)*n_refinement_pts + iii;
			double tmp_pY = SP_Del_pY_min + (double)iii * Delta_DpY;
			double tmp_result = exp(-tmp_pY) * gsl_cheb_eval (cs_accel_expEdNd3p, tmp_pY);
			refined_resonance_grids[tmp_index][(iqx * qynpts + iqy)*ntrig + itrig] = tmp_result;
			log_refined_grids[tmp_index][(iqx * qynpts + iqy)*ntrig + itrig] = log(abs(tmp_result)+1.e-100);
			sgn_refined_grids[tmp_index][(iqx * qynpts + iqy)*ntrig + itrig] = sgn(tmp_result);
		}
		++idx;
	}

	return;
}

void CorrelationFunction::Do_resonance_integrals(int parent_resonance_particle_id, int daughter_particle_id, int decay_channel, int iqt, int iqz)
{
	time_t rawtime;
  	struct tm * timeinfo;

	val11_arr = new double [qspace_cs_slice_length];
	val12_arr = new double [qspace_cs_slice_length];
	val21_arr = new double [qspace_cs_slice_length];
	val22_arr = new double [qspace_cs_slice_length];

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

	double loc_qz = qz_pts[iqz];
	double loc_qt = qt_pts[iqt];
	current_iqt = iqt;
	current_iqz = iqz;
	current_pY_shift = 0.5 * log(abs((loc_qt+loc_qz + 1.e-100)/(loc_qt-loc_qz + 1.e-100)));
	//current_pY_shift = 0.0;

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
					double Del_PKY = VEC_n2_P_Y[iv] - current_pY_shift;
					double PKphi = VEC_n2_PPhi_tilde[NB2_indexer(iv,izeta)];

					for (int tempidx = 0; tempidx <= 1; ++tempidx)
					{
						if (tempidx != 0)
							PKphi = VEC_n2_PPhi_tildeFLIP[NB2_indexer(iv,izeta)];		//also takes Pp --> Pm
						currentPpm = VEC_n2_Ppm[NB2_indexer(iv,izeta)*2 + tempidx];
						//spectra
						if ( iqt == (qtnpts - 1)/2 && iqz == (qznpts - 1)/2 )
							Edndp3(PKT, PKphi, &Csum);
						//space-time moments
						if (!IGNORE_LONG_LIVED_RESONANCES || Gamma >= hbarC / max_lifetime)
						{
							//if (iv==0 && izeta==0 && tempidx==0)
							//	cout << "Doing eiqxEdndp3(nb=2, parent_pid = " << parent_resonance_particle_id
							//			<< "): IGNORE_LONG_LIVED_RESONANCES = " << IGNORE_LONG_LIVED_RESONANCES << "::: "
							//			<< Gamma << "   " << hbarC / max_lifetime << "   " << max_lifetime << endl;
							Set_val_arrays(PKT, PKphi, Del_PKY);
							eiqxEdndp3(PKT, PKphi, Del_PKY, Csum_vec, local_verbose);
						}
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
			for (int itrig = 0; itrig < ntrig; ++itrig)
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
						double Del_PKY = VEC_n3_P_Y[is*n_v_pts+iv] - current_pY_shift;
						double PKphi = VEC_n3_PPhi_tilde[NB3_indexer(is,iv,izeta)];

						for (int tempidx = 0; tempidx <= 1; ++tempidx)
						{
							if (tempidx != 0)
								PKphi = VEC_n3_PPhi_tildeFLIP[NB3_indexer(is,iv,izeta)];		//also takes Pp --> Pm
							currentPpm = VEC_n3_Ppm[NB3_indexer(is,iv,izeta)*2+tempidx];
							//spectra
							if ( iqt == (qtnpts - 1)/2 && iqz == (qznpts - 1)/2 )
								Edndp3(PKT, PKphi, &Csum);
							//space-time moments
							if (!IGNORE_LONG_LIVED_RESONANCES || Gamma >= hbarC / max_lifetime)
							{
								//if (is==0 && iv==0 && izeta==0 && tempidx==0)
								//	cout << "Doing eiqxEdndp3(nb=3, parent_pid = " << parent_resonance_particle_id
								//			<< "): IGNORE_LONG_LIVED_RESONANCES = " << IGNORE_LONG_LIVED_RESONANCES << "::: "
								//			<< Gamma << "   " << hbarC / max_lifetime << "   " << max_lifetime << endl;
								Set_val_arrays(PKT, PKphi, Del_PKY);
								eiqxEdndp3(PKT, PKphi, Del_PKY, Csum_vec, local_verbose);
							}
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
			for (int itrig = 0; itrig < ntrig; ++itrig)
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
	for (int ipTpphipY = 0; ipTpphipY < n_pT_pts * n_pphi_pts * n_refinement_pts; ++ipTpphipY)
	{
		delete [] refined_resonance_grids[ipTpphipY];
		delete [] log_refined_grids[ipTpphipY];
		delete [] sgn_refined_grids[ipTpphipY];
	}
	delete [] refined_resonance_grids;
	delete [] log_refined_grids;
	delete [] sgn_refined_grids;

	delete [] val11_arr;
	delete [] val12_arr;
	delete [] val21_arr;
	delete [] val22_arr;

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
void CorrelationFunction::eiqxEdndp3(double ptr, double phir, double spyr, double * results, int loc_verb /*==0*/)
{
	double phi0, phi1, py0, py1;
	double val11, val12, val21, val22;	//store intermediate results of pT interpolation
	double val1, val2;					//store intermediate results of pphi interpolation
	
	double pyr = abs(spyr);				//used for checking
	double parity_factor = 1.0;
	if (abs(qt_pts[current_iqt]) < abs(qz_pts[current_iqz]))
		parity_factor = -1.0;		//used to get correct sign of imaginary part of spectra
									//with Del_pY --> -Del_pY for |qt| < |qz|
									//real part of spectra always symmetric about Del_pY == 0
	double reflection_factor = 1.0;
	
	bool pY_out_of_range = false;

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
	if(phir < SP_pphi[0])                   //if angle is less than minimum angle grid point
	{
		phi0 = SP_pphi[npphi_max] - 2. * M_PI;
		phi1 = SP_pphi[0];
		nphi = 0;
		nphim1 = npphi_max;
		pY_out_of_range = true;
	}
	else if(phir > SP_pphi[npphi_max])      //if angle is greater than maximum angle grid point
	{
		phi0 = SP_pphi[npphi_max];
		phi1 = SP_pphi[0] + 2. * M_PI;
		nphi = 0;
		nphim1 = npphi_max;
		pY_out_of_range = true;
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
	if(pyr > SP_Del_pY_max)	//if rapidity is greater than maximum rapidity grid point
	{
		py0 = SP_Del_pY_max;
		py1 = SP_Del_pY_max;
		npy = n_refinement_pts - 1;	//this just guarantees nearest neighbor if pY > pY_max
		npym1 = n_refinement_pts - 1;
	}
	else if(pyr < SP_Del_pY_min)	//else if rapidity is less than minimum rapidity grid point
	{
		py0 = SP_Del_pY_min;
		py1 = SP_Del_pY_min;
		npy = 0;	//this just guarantees nearest neighbor if pY > pY_max
		npym1 = 0;
	}
	else						//if rapidity is within grid range
	{
		npym1 = floor((pyr-SP_Del_pY_min)/Delta_DpY);
		npy = npym1 + 1;
		py0 = SP_Del_pY_min + (double)npym1 * Delta_DpY;
		py1 = SP_Del_pY_min + (double)npy * Delta_DpY;
	}

	if (pT0==pT1 || phi0==phi1)
	{
		cerr << "ERROR in eiqxEdndp3(): pT and/or pphi values equal!" << endl;
		exit(1);
	}

	bool ptr_greater_than_pT1 = (ptr > pT1);

	double one_by_pTdiff = 1./(pT1 - pT0), one_by_pphidiff = 1./(phi1 - phi0), one_by_pYdiff = 1./(py1 - py0 + 1.e-10);
	double del_ptr_pt0 = ptr - pT0, del_phir_phi0 = phir - phi0, del_pyr_py0 = pyr - py0;

	// set index for looping
	int qpt_cs_idx = 0;
	int qlist_idx = 0;

	if ( ptr_greater_than_pT1 )                             // if pT interpolation point is outside of grid
	{
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		{
			double arg = one_by_Gamma_Mres * dot_four_vectors(qlist[qlist_idx], currentPpm);
			double akr = 1./(1.+arg*arg);
			double aki = arg/(1.+arg*arg);

			for (int iCS = 0; iCS < 2; ++iCS)	//cosine (rapidity-even), then sine (rapidity-odd)
			{
				/////////////////////////////////////////////////////////////////
				// DO REAL PART
				/////////////////////////////////////////////////////////////////
				// interpolate over pT values
				/////////////////////////////////////////////////////////////////
				val11 = val11_arr[qpt_cs_idx];
				val21 = val21_arr[qpt_cs_idx];
				val12 = val12_arr[qpt_cs_idx];
				val22 = val22_arr[qpt_cs_idx];

				//////////////////////////////////////
				// interpolate val11 and val12 over the pphi direction to get val1
				// similarly for val21 and val22 --> val2
				//////////////////////////////////////
				val1 = lin_int(del_phir_phi0, one_by_pphidiff, val11, val21);
				val2 = lin_int(del_phir_phi0, one_by_pphidiff, val12, val22);

				// finally, get the interpolated value
				double Zkr = lin_int(del_pyr_py0, one_by_pYdiff, val1, val2);

			    /////////////////////////////////////////////////////////////////
			    // DO IMAGINARY PART
			    /////////////////////////////////////////////////////////////////
			    // interpolate over pT values
			    /////////////////////////////////////////////////////////////////
				val11 = val11_arr[qpt_cs_idx+1];
				val21 = val21_arr[qpt_cs_idx+1];
				val12 = val12_arr[qpt_cs_idx+1];
				val22 = val22_arr[qpt_cs_idx+1];

				//////////////////////////////////////
				// interpolate val11 and val12 over the pphi direction to get val1
				// similarly for val21 and val22 --> val2
				//////////////////////////////////////
				val1 = lin_int(del_phir_phi0, one_by_pphidiff, val11, val21);
				val2 = lin_int(del_phir_phi0, one_by_pphidiff, val12, val22);

				// finally, get the interpolated value
				double Zki = reflection_factor * lin_int(del_pyr_py0, one_by_pYdiff, val1, val2);

			    /////////////////////////////////////////////////////
			    // Finally, update results vectors appropriately
			    /////////////////////////////////////////////////////
			    //--> update the real part of weighted daughter spectra
			    results[qpt_cs_idx] += akr*Zkr-aki*Zki;
			    //--> update the imaginary part of weighted daughter spectra
			    results[qpt_cs_idx+1] += akr*Zki+aki*Zkr;

				//needed to exploit symmetries of sine component
				reflection_factor *= parity_factor;
			    qpt_cs_idx += 2;
			}
			qlist_idx++;
		}   //end of all q-loops
	}
	else if (ptr > PTCHANGE)                             // if pT interpolation point is larger than PTCHANGE (currently 1.0 GeV)
	{
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		{
			double arg = one_by_Gamma_Mres * dot_four_vectors(qlist[qlist_idx], currentPpm);
			double akr = 1./(1.+arg*arg);
			double aki = arg/(1.+arg*arg);

			for (int iCS = 0; iCS < 2; ++iCS)	//cosine (rapidity-even), then sine (rapidity-odd)
			{
				/////////////////////////////////////////////////////////////////
				// DO REAL PART
				/////////////////////////////////////////////////////////////////
				// interpolate over pT values
				/////////////////////////////////////////////////////////////////
				val11 = val11_arr[qpt_cs_idx];
				val21 = val21_arr[qpt_cs_idx];
				val12 = val12_arr[qpt_cs_idx];
				val22 = val22_arr[qpt_cs_idx];

				//////////////////////////////////////
				// interpolate val11 and val12 over the pphi direction to get val1
				// similarly for val21 and val22 --> val2
				//////////////////////////////////////
				val1 = lin_int(del_phir_phi0, one_by_pphidiff, val11, val21);
				val2 = lin_int(del_phir_phi0, one_by_pphidiff, val12, val22);

				// finally, get the interpolated value
				double Zkr = lin_int(del_pyr_py0, one_by_pYdiff, val1, val2);

			    /////////////////////////////////////////////////////////////////
			    // DO IMAGINARY PART
			    /////////////////////////////////////////////////////////////////
			    // interpolate over pT values
			    /////////////////////////////////////////////////////////////////
				val11 = val11_arr[qpt_cs_idx+1];
				val21 = val21_arr[qpt_cs_idx+1];
				val12 = val12_arr[qpt_cs_idx+1];
				val22 = val22_arr[qpt_cs_idx+1];

				//////////////////////////////////////
				// interpolate val11 and val12 over the pphi direction to get val1
				// similarly for val21 and val22 --> val2
				//////////////////////////////////////
				val1 = lin_int(del_phir_phi0, one_by_pphidiff, val11, val21);
				val2 = lin_int(del_phir_phi0, one_by_pphidiff, val12, val22);

				// finally, get the interpolated value
				double Zki = reflection_factor * lin_int(del_pyr_py0, one_by_pYdiff, val1, val2);

			    /////////////////////////////////////////////////////
			    // Finally, update results vectors appropriately
			    /////////////////////////////////////////////////////
			    //--> update the real part of weighted daughter spectra
			    results[qpt_cs_idx] += akr*Zkr-aki*Zki;
			    //--> update the imaginary part of weighted daughter spectra
			    results[qpt_cs_idx+1] += akr*Zki+aki*Zkr;

				//needed to exploit symmetries of sine component
				reflection_factor *= parity_factor;
			    qpt_cs_idx += 2;
			}
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

			for (int iCS = 0; iCS < 2; ++iCS)
			{
				/////////////////////////////////////////////////////////////////
				// DO REAL PART
				/////////////////////////////////////////////////////////////////
				// interpolate over pT values
				/////////////////////////////////////////////////////////////////
				val11 = val11_arr[qpt_cs_idx];
				val21 = val21_arr[qpt_cs_idx];
				val12 = val12_arr[qpt_cs_idx];
				val22 = val22_arr[qpt_cs_idx];

				val1 = lin_int(del_phir_phi0, one_by_pphidiff, val11, val21);
				val2 = lin_int(del_phir_phi0, one_by_pphidiff, val12, val22);

				double Zkr = lin_int(del_pyr_py0, one_by_pYdiff, val1, val2);

				/////////////////////////////////////////////////////////////////
				// DO IMAGINARY PART
				/////////////////////////////////////////////////////////////////
				// interpolate over pT values
				/////////////////////////////////////////////////////////////////
				val11 = val11_arr[qpt_cs_idx+1];
				val21 = val21_arr[qpt_cs_idx+1];
				val12 = val12_arr[qpt_cs_idx+1];
				val22 = val22_arr[qpt_cs_idx+1];

				val1 = lin_int(del_phir_phi0, one_by_pphidiff, val11, val21);
				val2 = lin_int(del_phir_phi0, one_by_pphidiff, val12, val22);

				double Zki = reflection_factor * lin_int(del_pyr_py0, one_by_pYdiff, val1, val2);

				/////////////////////////////////////////////////////
				// Finally, update results vectors appropriately
				/////////////////////////////////////////////////////
				//--> update the real part of weighted daughter spectra
				results[qpt_cs_idx] += akr*Zkr-aki*Zki;
				//--> update the imaginary part of weighted daughter spectra
				results[qpt_cs_idx+1] += akr*Zki+aki*Zkr;

				//needed to exploit symmetries of sine component
				reflection_factor *= parity_factor;
				qpt_cs_idx += 2;
			}
			qlist_idx++;
		}       //end of all q-loops
	}

	/*qpt_cs_idx = 0;
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
		cout << "SANITY CHECK: " << qt_pts[current_iqt] << "   " << qx_pts[iqx] << "   " << qy_pts[iqy] << "   " << qz_pts[current_iqz] << "   "
				<< results[qpt_cs_idx++] << "   " << results[qpt_cs_idx++] << "   " << results[qpt_cs_idx++] << "   " << results[qpt_cs_idx++] << endl;

if (1) exit(8);*/

	return;
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
void CorrelationFunction::Set_val_arrays(double ptr, double phir, double spyr)
{
	double phi0, phi1, py0, py1;
	double val11, val12, val21, val22;	//store intermediate results of pT interpolation
	double val1, val2;					//store intermediate results of pphi interpolation
	
	double pyr = abs(spyr);				//used for checking
	
	bool pY_out_of_range = false;

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
	if(phir < SP_pphi[0])                   //if angle is less than minimum angle grid point
	{
		phi0 = SP_pphi[npphi_max] - 2. * M_PI;
		phi1 = SP_pphi[0];
		nphi = 0;
		nphim1 = npphi_max;
		pY_out_of_range = true;
	}
	else if(phir > SP_pphi[npphi_max])      //if angle is greater than maximum angle grid point
	{
		phi0 = SP_pphi[npphi_max];
		phi1 = SP_pphi[0] + 2. * M_PI;
		nphi = 0;
		nphim1 = npphi_max;
		pY_out_of_range = true;
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
	if(pyr > SP_Del_pY_max)	//if rapidity is greater than maximum rapidity grid point
	{
		py0 = SP_Del_pY_max;
		py1 = SP_Del_pY_max;
		npy = n_refinement_pts - 1;	//this just guarantees nearest neighbor if pY > pY_max
		npym1 = n_refinement_pts - 1;
	}
	else if(pyr < SP_Del_pY_min)	//else if rapidity is less than minimum rapidity grid point
	{
		py0 = SP_Del_pY_min;
		py1 = SP_Del_pY_min;
		npy = 0;	//this just guarantees nearest neighbor if pY > pY_max
		npym1 = 0;
	}
	else						//if rapidity is within grid range
	{
		npym1 = floor((pyr-SP_Del_pY_min)/Delta_DpY);
		npy = npym1 + 1;
		py0 = SP_Del_pY_min + (double)npym1 * Delta_DpY;
		py1 = SP_Del_pY_min + (double)npy * Delta_DpY;
	}

	if (pT0==pT1 || phi0==phi1)
	{
		cerr << "ERROR in eiqxEdndp3(): pT and/or pphi values equal!" << endl;
		exit(1);
	}

	bool ptr_greater_than_pT1 = (ptr > pT1);

	double one_by_pTdiff = 1./(pT1 - pT0), one_by_pphidiff = 1./(phi1 - phi0), one_by_pYdiff = 1./(py1 - py0 + 1.e-10);
	double del_ptr_pt0 = ptr - pT0, del_phir_phi0 = phir - phi0, del_pyr_py0 = pyr - py0;

	double * f111_arr = refined_resonance_grids[((npt-1)*n_pphi_pts+nphim1)*n_refinement_pts + npym1];
	double * f112_arr = refined_resonance_grids[((npt-1)*n_pphi_pts+nphim1)*n_refinement_pts + npy];
	double * f121_arr = refined_resonance_grids[((npt-1)*n_pphi_pts+nphi)*n_refinement_pts + npym1];
	double * f122_arr = refined_resonance_grids[((npt-1)*n_pphi_pts+nphi)*n_refinement_pts + npy];
	double * f211_arr = refined_resonance_grids[(npt*n_pphi_pts+nphim1)*n_refinement_pts + npym1];
	double * f212_arr = refined_resonance_grids[(npt*n_pphi_pts+nphim1)*n_refinement_pts + npy];
	double * f221_arr = refined_resonance_grids[(npt*n_pphi_pts+nphi)*n_refinement_pts + npym1];
	double * f222_arr = refined_resonance_grids[(npt*n_pphi_pts+nphi)*n_refinement_pts + npy];

	// set index for looping
	int qpt_cs_idx = 0;

	if ( ptr_greater_than_pT1 )                             // if pT interpolation point is outside of grid
	{

		double * log_f111_arr = log_refined_grids[((npt-1)*n_pphi_pts+nphim1)*n_refinement_pts + npym1];
		double * log_f112_arr = log_refined_grids[((npt-1)*n_pphi_pts+nphim1)*n_refinement_pts + npy];
		double * log_f121_arr = log_refined_grids[((npt-1)*n_pphi_pts+nphi)*n_refinement_pts + npym1];
		double * log_f122_arr = log_refined_grids[((npt-1)*n_pphi_pts+nphi)*n_refinement_pts + npy];
		double * log_f211_arr = log_refined_grids[(npt*n_pphi_pts+nphim1)*n_refinement_pts + npym1];
		double * log_f212_arr = log_refined_grids[(npt*n_pphi_pts+nphim1)*n_refinement_pts + npy];
		double * log_f221_arr = log_refined_grids[(npt*n_pphi_pts+nphi)*n_refinement_pts + npym1];
		double * log_f222_arr = log_refined_grids[(npt*n_pphi_pts+nphi)*n_refinement_pts + npy];

		double * sign_of_f111_arr = sgn_refined_grids[((npt-1)*n_pphi_pts+nphim1)*n_refinement_pts + npym1];
		double * sign_of_f112_arr = sgn_refined_grids[((npt-1)*n_pphi_pts+nphim1)*n_refinement_pts + npy];
		double * sign_of_f121_arr = sgn_refined_grids[((npt-1)*n_pphi_pts+nphi)*n_refinement_pts + npym1];
		double * sign_of_f122_arr = sgn_refined_grids[((npt-1)*n_pphi_pts+nphi)*n_refinement_pts + npy];
		double * sign_of_f211_arr = sgn_refined_grids[(npt*n_pphi_pts+nphim1)*n_refinement_pts + npym1];
		double * sign_of_f212_arr = sgn_refined_grids[(npt*n_pphi_pts+nphim1)*n_refinement_pts + npy];
		double * sign_of_f221_arr = sgn_refined_grids[(npt*n_pphi_pts+nphi)*n_refinement_pts + npym1];
		double * sign_of_f222_arr = sgn_refined_grids[(npt*n_pphi_pts+nphi)*n_refinement_pts + npy];

		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int iCS = 0; iCS < 2; ++iCS)	//cos/sin
		for (int iRI = 0; iRI < 2; ++iRI)	//real/imag
		{
			double sign_of_f111 = sign_of_f111_arr[qpt_cs_idx];
			double sign_of_f112 = sign_of_f112_arr[qpt_cs_idx];
			double sign_of_f121 = sign_of_f121_arr[qpt_cs_idx];
			double sign_of_f122 = sign_of_f122_arr[qpt_cs_idx];
			double sign_of_f211 = sign_of_f211_arr[qpt_cs_idx];
			double sign_of_f212 = sign_of_f212_arr[qpt_cs_idx];
			double sign_of_f221 = sign_of_f221_arr[qpt_cs_idx];
			double sign_of_f222 = sign_of_f222_arr[qpt_cs_idx];

			// set val11
			if ( sign_of_f111 * sign_of_f211 > 0 && log_f211_arr[qpt_cs_idx] < log_f111_arr[qpt_cs_idx] )
				val11_arr[qpt_cs_idx] = sign_of_f111 * exp( lin_int(del_ptr_pt0, one_by_pTdiff, log_f111_arr[qpt_cs_idx], log_f211_arr[qpt_cs_idx]) );
			else
				val11_arr[qpt_cs_idx] = 0.0;

			// set val21
			if ( sign_of_f121 * sign_of_f221 > 0 && log_f221_arr[qpt_cs_idx] < log_f121_arr[qpt_cs_idx] )
				val21_arr[qpt_cs_idx] = sign_of_f121 * exp( lin_int(del_ptr_pt0, one_by_pTdiff, log_f121_arr[qpt_cs_idx], log_f221_arr[qpt_cs_idx]) );
			else
				val21_arr[qpt_cs_idx] = 0.0;

			// set val12
			if ( sign_of_f112 * sign_of_f212 > 0 && log_f212_arr[qpt_cs_idx] < log_f112_arr[qpt_cs_idx] )
				val12_arr[qpt_cs_idx] = sign_of_f112 * exp( lin_int(del_ptr_pt0, one_by_pTdiff, log_f112_arr[qpt_cs_idx], log_f212_arr[qpt_cs_idx]) );
			else
				val12_arr[qpt_cs_idx] = 0.0;

			// set val22
			if ( sign_of_f122 * sign_of_f222 > 0 && log_f222_arr[qpt_cs_idx] < log_f122_arr[qpt_cs_idx] )
				val22_arr[qpt_cs_idx] = sign_of_f122 * exp( lin_int(del_ptr_pt0, one_by_pTdiff, log_f122_arr[qpt_cs_idx], log_f222_arr[qpt_cs_idx]) );
			else
				val22_arr[qpt_cs_idx] = 0.0;
		    ++qpt_cs_idx;
		}
	}
	else if (ptr > PTCHANGE)                             // if pT interpolation point is larger than PTCHANGE (currently 1.0 GeV)
	{
		double * log_f111_arr = log_refined_grids[((npt-1)*n_pphi_pts+nphim1)*n_refinement_pts + npym1];
		double * log_f112_arr = log_refined_grids[((npt-1)*n_pphi_pts+nphim1)*n_refinement_pts + npy];
		double * log_f121_arr = log_refined_grids[((npt-1)*n_pphi_pts+nphi)*n_refinement_pts + npym1];
		double * log_f122_arr = log_refined_grids[((npt-1)*n_pphi_pts+nphi)*n_refinement_pts + npy];
		double * log_f211_arr = log_refined_grids[(npt*n_pphi_pts+nphim1)*n_refinement_pts + npym1];
		double * log_f212_arr = log_refined_grids[(npt*n_pphi_pts+nphim1)*n_refinement_pts + npy];
		double * log_f221_arr = log_refined_grids[(npt*n_pphi_pts+nphi)*n_refinement_pts + npym1];
		double * log_f222_arr = log_refined_grids[(npt*n_pphi_pts+nphi)*n_refinement_pts + npy];

		double * sign_of_f111_arr = sgn_refined_grids[((npt-1)*n_pphi_pts+nphim1)*n_refinement_pts + npym1];
		double * sign_of_f112_arr = sgn_refined_grids[((npt-1)*n_pphi_pts+nphim1)*n_refinement_pts + npy];
		double * sign_of_f121_arr = sgn_refined_grids[((npt-1)*n_pphi_pts+nphi)*n_refinement_pts + npym1];
		double * sign_of_f122_arr = sgn_refined_grids[((npt-1)*n_pphi_pts+nphi)*n_refinement_pts + npy];
		double * sign_of_f211_arr = sgn_refined_grids[(npt*n_pphi_pts+nphim1)*n_refinement_pts + npym1];
		double * sign_of_f212_arr = sgn_refined_grids[(npt*n_pphi_pts+nphim1)*n_refinement_pts + npy];
		double * sign_of_f221_arr = sgn_refined_grids[(npt*n_pphi_pts+nphi)*n_refinement_pts + npym1];
		double * sign_of_f222_arr = sgn_refined_grids[(npt*n_pphi_pts+nphi)*n_refinement_pts + npy];

		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int iCS = 0; iCS < 2; ++iCS)	//cos/sin
		for (int iRI = 0; iRI < 2; ++iRI)	//real/imag
		{
			double sign_of_f111 = sign_of_f111_arr[qpt_cs_idx];
			double sign_of_f112 = sign_of_f112_arr[qpt_cs_idx];
			double sign_of_f121 = sign_of_f121_arr[qpt_cs_idx];
			double sign_of_f122 = sign_of_f122_arr[qpt_cs_idx];
			double sign_of_f211 = sign_of_f211_arr[qpt_cs_idx];
			double sign_of_f212 = sign_of_f212_arr[qpt_cs_idx];
			double sign_of_f221 = sign_of_f221_arr[qpt_cs_idx];
			double sign_of_f222 = sign_of_f222_arr[qpt_cs_idx];

			// set val11
			if (sign_of_f111 * sign_of_f211 > 0) // if the two points have the same sign in the pT direction, interpolate logs
				val11_arr[qpt_cs_idx] = sign_of_f111 * exp( lin_int(del_ptr_pt0, one_by_pTdiff, log_f111_arr[qpt_cs_idx], log_f211_arr[qpt_cs_idx]) );
			else                                    // otherwise, just interpolate original vals
				val11_arr[qpt_cs_idx] = lin_int(del_ptr_pt0, one_by_pTdiff, f111_arr[qpt_cs_idx], f211_arr[qpt_cs_idx]);

			// set val21
			if (sign_of_f121 * sign_of_f221 > 0) // if the two points have the same sign in the pT direction, interpolate logs
				val21_arr[qpt_cs_idx] = sign_of_f121 * exp( lin_int(del_ptr_pt0, one_by_pTdiff, log_f121_arr[qpt_cs_idx], log_f221_arr[qpt_cs_idx]) );
			else                                    // otherwise, just interpolate original vals
				val21_arr[qpt_cs_idx] = lin_int(del_ptr_pt0, one_by_pTdiff, f121_arr[qpt_cs_idx], f221_arr[qpt_cs_idx]);

			// set val12
			if (sign_of_f112 * sign_of_f212 > 0) // if the two points have the same sign in the pT direction, interpolate logs
				val12_arr[qpt_cs_idx] = sign_of_f112 * exp( lin_int(del_ptr_pt0, one_by_pTdiff, log_f112_arr[qpt_cs_idx], log_f212_arr[qpt_cs_idx]) );
			else                                    // otherwise, just interpolate original vals
				val12_arr[qpt_cs_idx] = lin_int(del_ptr_pt0, one_by_pTdiff, f112_arr[qpt_cs_idx], f212_arr[qpt_cs_idx]);

			// set val22
			if (sign_of_f122 * sign_of_f222 > 0) // if the two points have the same sign in the pT direction, interpolate logs
				val22_arr[qpt_cs_idx] = sign_of_f122 * exp( lin_int(del_ptr_pt0, one_by_pTdiff, log_f122_arr[qpt_cs_idx], log_f222_arr[qpt_cs_idx]) );
			else                                    // otherwise, just interpolate original vals
				val22_arr[qpt_cs_idx] = lin_int(del_ptr_pt0, one_by_pTdiff, f122_arr[qpt_cs_idx], f222_arr[qpt_cs_idx]);
			++qpt_cs_idx;
		}
	}
	else                                            // if pT is smaller than PTCHANGE, just use linear interpolation, no matter what
	{
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int iCS = 0; iCS < 2; ++iCS)	//cos/sin
		for (int iRI = 0; iRI < 2; ++iRI)	//real/imag
		{
			val11_arr[qpt_cs_idx] = lin_int(del_ptr_pt0, one_by_pTdiff, f111_arr[qpt_cs_idx], f211_arr[qpt_cs_idx]);
			val21_arr[qpt_cs_idx] = lin_int(del_ptr_pt0, one_by_pTdiff, f121_arr[qpt_cs_idx], f221_arr[qpt_cs_idx]);
			val12_arr[qpt_cs_idx] = lin_int(del_ptr_pt0, one_by_pTdiff, f112_arr[qpt_cs_idx], f212_arr[qpt_cs_idx]);
			val22_arr[qpt_cs_idx] = lin_int(del_ptr_pt0, one_by_pTdiff, f122_arr[qpt_cs_idx], f222_arr[qpt_cs_idx]);
			++qpt_cs_idx;
		}
	}

	return;
}



//End of file
