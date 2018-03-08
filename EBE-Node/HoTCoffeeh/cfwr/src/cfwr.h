#ifndef CFWR_H
#define CFWR_H

#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<fstream>
#include<vector>
#include<set>
#include<queue>
#include <fenv.h>

#include "H5Cpp.h"

#include "readindata.h"
#include "parameters.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"
#include "chebyshev.h"
#include "ParameterReader.h"

using namespace std;

typedef struct
{
	int resonance_particle_id;		// keeps track of current resonance's index in all_particles array
	int resonance_idx;			// keeps track of current resonance's index in chosen_resonances vector
	int nbody;
	int resonance_sign;
	double resonance_mass;
	double resonance_mu;
	double resonance_gspin;
	double resonance_Gamma;
	double resonance_total_br;
	double resonance_direct_br;
	double * resonance_decay_masses;
	int * resonance_decay_monvals;
	double * resonance_decay_Gammas;
	string resonance_name;
	bool include_channel;
} decay_info;

struct Correlationfunction3D_data
{
	size_t data_length;
	double * q_o;
	double * q_s;
	double * q_l;
	double * y;
	double * sigma;
};

int Fittarget_correlfun3D_f (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr);
int Fittarget_correlfun3D_df (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr);
int Fittarget_correlfun3D_fdf (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr);
int Fittarget_correlfun3D_f_withlambda (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr);
int Fittarget_correlfun3D_df_withlambda (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr);
int Fittarget_correlfun3D_fdf_withlambda (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr);

class CorrelationFunction
{
	private:
		ParameterReader * paraRdr;

		int USE_PLANE_PSI_ORDER;
		int INCLUDE_DELTA_F;
		int GROUPING_PARTICLES;
		double PARTICLE_DIFF_TOLERANCE;
		//int USE_LAMBDA;
		//int USE_LOG_FIT;
		int CALCULATE_CF_MODE;
		int USE_EXTRAPOLATION;
		int IGNORE_LONG_LIVED_RESONANCES;
		int FIT_WITH_PROJECTED_CFVALS;
		int FLESH_OUT_CF;
		int n_order;
		double tol;
		int flagneg;
		double max_lifetime;

		//header info
		int n_pT_pts, n_pphi_pts, n_pY_pts, nKT, nKphi;
		int qtnpts, qxnpts, qynpts, qznpts;
		double KT_min, KT_max;
		double init_qt, init_qx, init_qy, init_qz;
		double delta_qt, delta_qx, delta_qy, delta_qz;

		//particle information 
		string particle_name;
		double particle_mass;
		int particle_monval;
		int particle_id;     //particle id
		double particle_sign;   //+/- 1 for Fermi/Bose statistics for baryon/meson
		double particle_gspin;  //particle degeneracy 
		double particle_mu;
		double current_total_resonance_percentage, previous_total_resonance_percentage;
		particle_info * all_particles;
		vector<int> chosen_resonances;
		bool thermal_pions_only;
		int Nparticle, NchosenParticle;
		int target_particle_id;		//the particle whose spectra (with resonance contributions) you want to compute
		int current_level_of_output;
		int qspace_cs_slice_length;
		int full_FO_length;
		int FO_length;
		int n_alpha_points;
		//int nFO_cutoff;
		//int number_of_percentage_markers;
		double qtmax, qx_max, qy_max, qz_max;
		int new_nqxpts, new_nqypts, new_nqzpts;
		double q_space_CF_cutoff;		// when correlator falls below this value,
							//	set correlator to zero for any q-points further away from q-origin than that
		double ** current_q_space_cutoff;	// point in q-space at which cutoff of CF begins (depends on pT and pphi)
				
		//array to hold previous and current resonance info
		decay_info * decay_channels;
		int current_parent_resonance;
		int current_decay_channel_idx, current_resonance_particle_id, previous_resonance_idx, current_resonance_idx, current_reso_nbody;
		double current_resonance_mu, current_resonance_mass, current_resonance_Gamma, current_m2_Gamma, current_m3_Gamma;
		double current_resonance_total_br, current_resonance_direct_br, current_daughter_mass, current_daughter_Gamma;
		double * current_resonance_decay_masses, * P_eval, * alpha_mu;
		int previous_decay_channel_idx, previous_resonance_particle_id, previous_reso_nbody;
		double previous_resonance_mu, previous_resonance_mass, previous_resonance_Gamma, previous_m2_Gamma, previous_m3_Gamma;
		double previous_resonance_total_br, previous_resonance_direct_br, previous_daughter_mass, previous_daughter_Gamma;
		double * previous_resonance_decay_masses;
		
		//arrays to hold results of resonance phase-space integrations
		double * current_dN_dypTdpTdphi_moments;
		double ** current_daughters_dN_dypTdpTdphi_moments;
		double * thermal_target_dN_dypTdpTdphi_moments;
		double * full_target_dN_dypTdpTdphi_moments;
		double * thermal_target_Yeq0_moments;
		double * full_target_Yeq0_moments;

		// needed these to avoid too many trigonometric evaluations
		double ** oscx, ** oscy;

		// skip cells with negative FOcell_density
		vector<vector<int> > FOcells_to_include;
	
		//needed for resonance calculations
		//	kinematic info
		double pstar, Estar, Yp, Ym, DeltaY, MTbar, DeltaMT, MTp, MTm, Qfunc;
		//	pair momentum info, currently assumes pT != 0
		double p_y, pT, pphi, mT, mass, ch_p_y, sh_p_y;
		//	resonance momentum info
		double P_Y, PT, PPhi, MT, Mres, PPhip, PPhim, m2, m3, Gamma, br, m2Gamma, m3Gamma, one_by_Gamma_Mres;
		double * Pp, * Pm, * currentPpm;

		//SP momentum arrays for interpolation grid
		double * SP_pT, * SP_pphi;
		double * SP_pT_wts, * SP_pphi_wts;
		double * sin_SP_pphi, * cos_SP_pphi;
		double ** SP_p0, ** SP_pz;
		double * plane_angle;
		double adjusted_SP_Del_pY_minimum;

		//Freeze-out surface information
		FO_surf* FOsurf_ptr;
		double Tdec, Edec, Pdec, muRES, signRES, gRES, S_prefactor;
	
		//single particle spectra for plane angle determination
		double SP_p_y, mean_pT;

		//pair momentum
		double K_y, ch_K_y, sh_K_y;
		double current_K_phi, cos_cKphi, sin_cKphi;
		double beta_perp, beta_l;
		double * K_T, * K_phi, * K_phi_weight;
		    
		//momentum rapidity grid
		double * eta_s, * eta_s_weight, * ch_eta_s, * sh_eta_s;
		double * base_Del_eta_s, * base_Del_eta_s_weight;
		double * SP_Del_pY, * ch_SP_pY, * sh_SP_pY;
		double * chebTcfs;
		double ** chebyshev_a_cfs, ** refined_resonance_grids, ** log_refined_grids, ** sgn_refined_grids;	//for resonance interpolation
		double ** exp_table_11, ** exp_table_21, ** exp_table_12, ** exp_table_22;
		bool * grids_calculated;

		//points and weights for resonance integrals
		double v_min, v_max, zeta_min, zeta_max, s_min, s_max;
		double * zeta_pts, * v_pts, * s_pts;
		double * zeta_wts, * v_wts, * s_wts;

		//some arrays to save unnecessary multiple calculations for resonances
		//	use these for n_body = 2
		double VEC_n2_spt, VEC_n2_g_s, VEC_n2_s_factor;
		double * VEC_n2_P_Y, * VEC_n2_PPhi_tilde, * VEC_n2_PPhi_tildeFLIP, * VEC_n2_PT;
		double * VEC_n2_v_factor, * VEC_n2_zeta_factor;
		double ** VEC_n2_Ppm;
		//	use these for n_body = 3
		double * VEC_n3_s_factor, * VEC_n3_g_s, * VEC_n3_P_Y, * VEC_n3_v_factor;
		double * VEC_n3_PPhi_tilde, * VEC_n3_PPhi_tildeFLIP, * VEC_n3_PT, * VEC_n3_zeta_factor;
		double ** VEC_n3_Ppm;
		double * ssum_vec, * vsum_vec, * zetasum_vec, * Csum_vec;
		
		double *** spectra, *** thermal_spectra, *** log_spectra, *** sign_spectra;

		//use this to hold full spectra (after various q-space symmetries have been exploited)
		double * reflected_moments;
		
		// relative momentum information
		double * qo_pts, * qs_pts, * ql_pts, * q_pts, * q_axes, * qt_pts, * qx_pts, * qy_pts, * qz_pts;
		double * q_out, * q_side, * q_long;
		int q1npts, q2npts, q3npts;		//123 indexing allows these to refer to either q[osl]npts or q[xyz]npts
		double * q1_pts, * q2_pts, * q3_pts;
		int iqt0, iqx0, iqy0, iqz0, ipY0;
		vector<vector<int> > sorted_q_pts_list;
		double ** qlist, * current_qlist_slice;
		vector<vector<int> > q_axes_and_rays;
		
		//store correlation functions
		//double *** Correl_3D;
		double ***** CFvals, ***** thermalCFvals, ***** crosstermCFvals, ***** resonancesCFvals;
		double *** fleshed_out_CF, *** fleshed_out_thermal, *** fleshed_out_crossterm, *** fleshed_out_resonances;
		double *** Correl_3D_err;
		double ** lambda_Correl, ** lambda_Correl_err;
		double ** lambda_QM;
		int *** correlator_minus_one_cutoff_norms;
		double * qx_fleshed_out_pts, * qy_fleshed_out_pts, * qz_fleshed_out_pts;

		//HBT radii coefficients
		double ** R2_side_GF, ** R2_out_GF, ** R2_long_GF, ** R2_outside_GF, ** R2_sidelong_GF, ** R2_outlong_GF;
		double ** R2_side_GF_C, ** R2_out_GF_C, ** R2_long_GF_C, ** R2_outside_GF_C, ** R2_sidelong_GF_C, ** R2_outlong_GF_C;
		double ** R2_side_GF_S, ** R2_out_GF_S, ** R2_long_GF_S, ** R2_outside_GF_S, ** R2_sidelong_GF_S, ** R2_outlong_GF_S;
		double ** R2_side_err, ** R2_out_err, ** R2_long_err, ** R2_outside_err, ** R2_sidelong_err, ** R2_outlong_err;

		double ** R2_side_QM, ** R2_out_QM, ** R2_long_QM, ** R2_outside_QM, ** R2_sidelong_QM, ** R2_outlong_QM;
		double ** R2_side_QM_C, ** R2_out_QM_C, ** R2_long_QM_C, ** R2_outside_QM_C, ** R2_sidelong_QM_C, ** R2_outlong_QM_C;
		double ** R2_side_QM_S, ** R2_out_QM_S, ** R2_long_QM_S, ** R2_outside_QM_S, ** R2_sidelong_QM_S, ** R2_outlong_QM_S;

		double ** res_sign_info, ** res_log_info, ** res_moments_info;
		double ** spec_sign_info, ** spec_log_info, ** spec_vals_info;
		
		//miscellaneous
		ofstream * global_out_stream_ptr;
		string path;
		string no_df_stem;
		int n_resonance, n_decay_channels;
		int n_body;
		double global_plane_psi;
		set<int> daughter_resonance_indices;

		int current_ipT, current_ipphi, current_ipY, current_iqt, current_iqz;
		double current_pY_shift;

		//some private methods		
		bool particles_are_the_same(int idx1, int idx2);
		bool Search_for_similar_particle(int dc_idx, int * result);

	public:
		//library of inline functions
		inline int indexer(const int ipt, const int ipphi, const int iqt, const int iqx, const int iqy, const int iqz, const int itrig);
		inline int indexer(const int ipt, const int ipphi, const int ipY, const int iqt, const int iqx, const int iqy, const int iqz, const int itrig);
		inline int indexer2(const int ipt, const int ipphi, const int ipY, const int iqt, const int iqx, const int iqy, const int iqz);
		inline int indexer3(const int ipt, const int ipphi, const int isurf);
		inline int mom_indexer(const int ipt, const int ipphi, const int ipY);
		inline int fixQTQZ_indexer(const int ipT, const int ipphi, const int ipY, const int iqx, const int iqy, const int itrig);
		inline int indexer4(const int ipt, const int ipphi, const int iqx, const int iqy);
		inline int FM_indexer(const int ipY, const int iqt, const int iqx, const int iqy, const int iqz);
		inline int HDF_indexer(const int ir, const int iqt, const int iqz);
		inline int NB2_indexer(const int iv, const int izeta);
		inline int NB3_indexer(const int is, const int iv, const int izeta);
		inline int qT_trig_indexer(const int iqx, const int iqy, const int iCSlong, const int iCStrans);
		inline void addElementToQueue(priority_queue<pair<double, size_t> >& p, pair<double, size_t> elem, size_t max_size);
		inline void set_to_zero(double * array, size_t arraylength);
		inline double dot_four_vectors(double * a, double * b);

		void Fourier_transform_emission_function(int iqt, int iqz);
		void Compute_phase_space_integrals(int iqt, int iqz);
		void Update_sourcefunction(particle_info* particle, int FOarray_length, int particle_idx);
		bool fexists(const char *filename);

		void Reset_FOcells_array();
		void Dump_FOcells(int local_pid);
		void Load_FOcells(int local_pid);

		// HDF routines
		void Initialize_HDF_resonance_array();
		void Reset_HDF_resonance_array();
		void Close_HDF_resonance_array();
		void Initialize_HDF_target_thermal_array();
		void Initialize_HDF_target_full_array();
		// resonances
		int Access_resonance_in_HDF_array(int local_pid, int iqt, int iqz, int access_mode, double * resonance_array_to_fill, bool verbose = false);
		int Administrate_resonance_HDF_array(int administration_mode);
		int Copy_chunk(int current_resonance_index, int reso_idx_to_be_copied);
		// Bessel coefficients
		int Access_besselcoeffs_in_HDF_array(int ipY, int access_mode, double * besselcoeffs_array_to_fill, int particle_mode = 0);
		int Administrate_besselcoeffs_HDF_array(int administration_mode, int particle_mode = 0);
		// target thermal moments...
		int Administrate_target_thermal_HDF_array(int administration_mode);
		int Access_target_thermal_in_HDF_array(int iqt, int iqz, int access_mode, double * target_full_array_to_fill, bool verbose = false);
		int Administrate_target_full_HDF_array(int administration_mode);
		int Access_target_full_in_HDF_array(int iqt, int iqz, int access_mode, double * target_full_array_to_fill, bool verbose = false);

		void Set_dN_dypTdpTdphi_moments(int local_pid, int iqt, int iqz);
		void Set_all_Bessel_grids(int iqt, int iqz, int particle_mode = 0);
		void Set_Y_eq_0_Bessel_grids(int iqt, int iqz, double * BC_chunk);
		void Set_target_moments(int iqt, int iqz);
		void Set_thermal_target_moments(int iqt, int iqz);
		void Set_full_target_moments(int iqt, int iqz);
		void Set_giant_arrays(int iqt, int iqx, int iqy, int iqz);
		void Cal_dN_dypTdpTdphi_no_weights(int local_pid);
		void Cal_dN_dypTdpTdphi_with_weights(
					int local_pid, int ipY, int iqt, int iqz,
					double * BC_chunk, int local_part_mode);
		void Cal_dN_dypTdpTdphi_with_weights_Yeq0_alternate(int iqt, int iqz);
		void Cal_dN_dypTdpTdphi_no_weights_Yeq0_alternate();
		void Cal_dN_dypTdpTdphi_with_weights_function_approx(
					int local_pid, double pT, double pphi, double pY,
					double qt, double qx, double qy, double qz,
					double * cosqx_dN_dypTdpTdphi, double * sinqx_dN_dypTdpTdphi);
		void Cal_dN_dypTdpTdphi_with_weights_function_approx(
					int local_pid, double pT, double pphi, double p_Y,
					double qt, double qx, double qy, double qz,
					double * cosLcosT_dN_dypTdpTdphi, double * cosLsinT_dN_dypTdpTdphi,
					double * sinLcosT_dN_dypTdpTdphi, double * sinLsinT_dN_dypTdpTdphi);
		void Cal_dN_dypTdpTdphi_with_weights_function_etas_integ(
					int local_pid, double pT, double pphi, double p_Y,
					double qt, double qx, double qy, double qz,
					double * cosLcosT_dN_dypTdpTdphi, double * cosLsinT_dN_dypTdpTdphi,
					double * sinLcosT_dN_dypTdpTdphi, double * sinLsinT_dN_dypTdpTdphi,
					double use_Boltzmann_approx = 0.0);
		void Cal_dN_dypTdpTdphi_with_weights_function_and_decay_etas_integ(
					int local_pid, double pT, double pphi, double Del_p_Y,
					double qt, double qx, double qy, double qz,
					double * cosLcosT_dN_dypTdpTdphi, double * cosLsinT_dN_dypTdpTdphi,
					double * sinLcosT_dN_dypTdpTdphi, double * sinLsinT_dN_dypTdpTdphi,
					double * res_RE, double * res_IM );
		void Cal_dN_dypTdpTdphi_no_weights_toy(int local_pid);
		void Cal_dN_dypTdpTdphi_with_weights_toy(
					int local_pid, int iqt, int iqz, int ipY,
					double * moments_to_update);
		void Cal_dN_dypTdpTdphi_with_weights_toy_func(
					int local_pid, double pT, double pphi, double pY,
					double qt, double qx, double qy, double qz,
					double * mom_CC, double * mom_CS,
					double * mom_SC, double * mom_SS);
		void Do_resonance_integrals(int parent_resonance_particle_id, int daughter_particle_id, int decay_channel, int iqt, int iqz);
		void Clear_and_set_exp_table_nb2();
		void Clear_and_set_exp_table_nb3();
		void Tabulate_resonance_Chebyshev_coefficients(int parent_resonance_particle_id);
		void Refine_resonance_grids(int parent_resonance_particle_id);
		void Set_current_daughter_info(int dc_idx, int daughter_idx);
		void Set_current_particle_info(int dc_idx);
		void Set_target_pphiavgd_CFs();
		bool Do_this_decay_channel(int dc_idx);
		bool Do_this_daughter_particle(int dc_idx, int daughter_idx, int * daughter_resonance_pid);
		void Get_spacetime_moments(int dc_idx, int iqt, int iqz);
		void Recycle_spacetime_moments();
		void Load_resonance_and_daughter_spectra(int local_pid, int iqt, int iqz);
		void Update_daughter_spectra(int local_pid, int iqt, int iqz);
		void Set_spectra_logs_and_signs(int local_pid);
		void Allocate_decay_channel_info();
		void Load_decay_channel_info_nb2(int dc_idx, double K_T_local, double K_phi_local, double K_y_local);
		void Load_decay_channel_info_nb3(int dc_idx, double K_T_local, double K_phi_local, double K_y_local);
		void Delete_decay_channel_info();
		int Set_daughter_list(int parent_resonance_index);

		void Fill_out_pts(double * pointsarray, int numpoints, double max_val, int spacing_type);

		void Reflect_in_qz_and_qt();

		//miscellaneous
		void Set_path(string path_in);
		void Set_use_delta_f();
		void Set_FOsurf_ptr(FO_surf* FOsurf_ptr_in, int FO_length_in);
		void Set_eiqx_matrices();

		double get_Q();
		double g(double s);
		double place_in_range(double phi, double min, double max);
		void Get_current_decay_string(int dc_idx, string * decay_string);
		int lookup_resonance_idx_from_particle_id(int particle_id);
		int list_daughters(int parent_resonance_index, set<int> * daughter_resonance_indices_ptr, particle_info * particle, int Nparticle);
		void eiqxEdndp3(double ptr, double phir, double pyr, double * results, int loc_verb = 0);
		void Set_val_arrays(double ptr, double phir, double spyr);
		void Edndp3(double ptr, double phir, double * result, int loc_verb = 0);
		void Set_correlation_function_q_pts();
		void Set_q_points();
		void Set_qlist(int iqt, int iqz);
		void Set_sorted_q_pts_list();
		void Get_q_points(double qo, double qs, double ql, double KT, double Kphi, double * qgridpts);
		void Allocate_resonance_running_sum_vectors();
		void Delete_resonance_running_sum_vectors();
		void Zero_resonance_running_sum_vector(double * vec);
		void Setup_temp_arrays(double ***** local_temp_moments, double ******* temp_moments_array);
		void Teardown_temp_arrays(double ***** local_temp_moments, double ******* temp_moments_array);
		void Setup_current_daughters_dN_dypTdpTdphi_moments(int n_daughter);
		void Cleanup_current_daughters_dN_dypTdpTdphi_moments(int n_daughter);
		void Allocate_osc_arrays(int FOarray_length);
		void Delete_osc_arrays();
		void R2_Fourier_transform(int ipt, double plane_psi, int mode);

		// Gaussian fit / correlation function routines
		void Allocate_CFvals();
		void Delete_CFvals();
		void Allocate_fleshed_out_CF();
		void Delete_fleshed_out_CF();
		void Flesh_out_CF(int ipt, int ipphi, double sample_scale = 1.0);
		double interpolate_CF(double *** current_C_slice, double qx0, double qy0, double qz0, int ipt, int thermal_or_resonances);
		double interpolate_qi(double q0, double qi0, double qi1, double f1, double f2, bool use_linear);
		void Get_GF_HBTradii();
		double get_CF(int ipt, int ipphi, int iqt, int iqx, int iqy, int iqz, bool return_projected_value);
		void get_CF_terms(double * totalresult, double * thermalresult, double * crosstermresult, double * resonanceresult,
									int ipt, int ipphi, int iqt, int iqx, int iqy, int iqz, bool return_projected_value);
		void Compute_correlationfunction(double * totalresult, double * thermalresult, double * crosstermresult, double * resonanceresult,
										int ipt, int ipphi, int iqx, int iqy, int iqz, double qt_interp, int interp_flag = 0, bool project_CF = true);
		void Cal_correlationfunction(bool project_CF = true);
		void find_minimum_chisq_correlationfunction_full(double *** Correl_3D, int ipt, int ipphi, bool fleshing_out_CF = true);
		void gsl_polynomial_fit(const vector<double> &data_x, const vector<double> &data_y, double * results, const int order, double * chisq, const int n);

		// input and output function prototypes
		void Output_total_target_dN_dypTdpTdphi();
		void Output_total_target_eiqx_dN_dypTdpTdphi(double current_fraction = -1.0);
		void Readin_total_target_eiqx_dN_dypTdpTdphi();
		void Output_total_eiqx_dN_dypTdpTdphi(int local_pid);
		void Output_total_eiqx_dN_dypTdpTdphi(int local_pid, int iqt, int iqz);
		void Output_thermal_target_eiqx_dN_dypTdpTdphi(int iqt, int iqz);
		void Output_results(int mode);
		void Readin_results(int mode);
		void Read_in_all_dN_dypTdpTdphi();
		void Output_chosen_resonances();
		void Output_resonance_fraction();
		void Read_in_correlationfunction();
		void Output_correlationfunction(bool project_CF = true);
		void Output_lambdas();
		void Output_fleshed_out_correlationfunction(int ipt, int ipphi, bool project_CF = true);
		void Dump_spectra_array(string output_filename, double *** array_to_dump);
		void Load_spectra_array(string output_filename, double *** array_to_read);

		//parameters that the user is free to define
		double plumberg_test_variable;
		bool use_delta_f;
		bool append_output;
		int n_events;
		vector<int> osr;
		int initial_event, currentfolderindex;
		bool read_in_all_dN_dypTdpTdphi, output_all_dN_dypTdpTdphi;
		double fraction_of_resonances;

		// need to hold giant array stuff
		//thermal target array
		H5::DataSpace * tta_dataspace, * tta_memspace;
		H5::H5File * tta_file;
		H5::DataSet * tta_dataset;
		//full target array
		H5::DataSpace * tfa_dataspace, * tfa_memspace;
		H5::H5File * tfa_file;
		H5::DataSet * tfa_dataset;
		//all resonance array
		H5::DataSpace * resonance_dataspace, * resonance_memspace;
		H5::H5File * resonance_file;
		H5::DataSet * resonance_dataset;
		//bessel coefficients array
		H5::DataSpace * besselcoeffs_dataspace, * besselcoeffs_memspace;
		H5::H5File * besselcoeffs_file;
		H5::DataSet * besselcoeffs_dataset;

		CorrelationFunction(ParameterReader * paraRdr_in, particle_info* particle, particle_info* all_particles_in, int Nparticle,
				vector<int> chosen_resonances, int particle_idx, ofstream& myout);
		~CorrelationFunction();

		//misc
		double S_x_p( int local_pid, int isurf, double eta_s,
						double pT, double pphi, double pY );
		int print_fit_state_3D (size_t iteration, gsl_multifit_fdfsolver * solver_ptr);
		void Fit_Correlationfunction3D_withlambda(double *** Correl_3D, int ipt, int ipphi, bool fleshing_out_CF = true);
		int print_fit_state_3D_withlambda (size_t iteration, gsl_multifit_fdfsolver * solver_ptr);
		inline double get_fit_results(int i, gsl_multifit_fdfsolver * solver_ptr);
		inline double get_fit_err (int i, gsl_matrix * covariance_ptr);


};

#endif
