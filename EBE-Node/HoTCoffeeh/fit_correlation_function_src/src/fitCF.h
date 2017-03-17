#ifndef FITCF_H
#define FITCF_H

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

#include<gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>            // gsl random number generators
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>

#include "readindata.h"
#include "parameters.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"
#include "chebyshev.h"

using namespace std;

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

class FitCF
{
	private:
		//header info
		int n_interp_pT_pts, n_interp_pphi_pts;
		int qtnpts, qxnpts, qynpts, qznpts;
		double init_qt, init_qx, init_qy, init_qz;

		particle_info * particle;
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
		particle_info * chosen_particle;
		vector<int> chosen_resonances, chosen_events;
		bool thermal_pions_only;
		int Nparticle, NchosenParticle;
		int target_particle_id;		//the particle whose spectra (with resonance contributions) you want to compute
		int current_level_of_output;
		int qspace_cs_slice_length;

		int nEvents;	//number of events in ensemble; default is 1

		double q_space_CF_cutoff;		// when correlator falls below this value,
							//	set correlator to zero for any q-points further away from q-origin than that
		double ** current_q_space_cutoff;	// point in q-space at which cutoff of CF begins (depends on pT and pphi)
						
		//arrays to hold results of resonance phase-space integrations
		//double ******* thermal_target_dN_dypTdpTdphi_moments;
		//double ******* full_target_dN_dypTdpTdphi_moments;
		//double ******* current_dN_dypTdpTdphi_moments;
		double * thermal_target_dN_dypTdpTdphi_moments;
		double * full_target_dN_dypTdpTdphi_moments;

		double ***** target_pphiavgd_CFs;
		double ***** target_pphivar_CFs;

		//SP momentum arrays for interpolation grid
		double * SPinterp_pT, * SPinterp_pphi;
		double * SPinterp_pT_wts, * SPinterp_pphi_wts;
		double * sin_SPinterp_pphi, * cos_SPinterp_pphi;
		double ** SPinterp_p0, ** SPinterp_pz;

		//pair momentum
		double K_y, ch_K_y, sh_K_y;
		double current_K_phi, cos_cKphi, sin_cKphi;
		double beta_perp, beta_l;
		double * K_T, * K_phi, * K_phi_weight;
		    
		//spatial rapidity grid
		double * eta_s, * ch_eta_s, * sh_eta_s, * eta_s_weight;

		Chebyshev * approx_R2s, * approx_R2o, * approx_R2l, * approx_R2os, * approx_R2sl, * approx_R2ol;
		double * flat_spectra;
		double ***** tmp_moments_real;
		double ***** tmp_moments_imag;
		
		double *** spectra, *** abs_spectra, *** thermal_spectra, *** log_spectra, *** sign_spectra;
		
		// relative momentum information
		double * qo_pts, * qs_pts, * ql_pts, * q_pts, * q_axes, * qt_pts, * qx_pts, * qy_pts, * qz_pts;
		double * q_out, * q_side, * q_long;
		int q1npts, q2npts, q3npts;		//123 indexing allows these to refer to either q[osl]npts or q[xyz]npts
		double * q1_pts, * q2_pts, * q3_pts;
		int iqt0, iqx0, iqy0, iqz0;
		vector<vector<int> > sorted_q_pts_list;
		double ** qlist, * current_qlist_slice;
		vector<vector<int> > q_axes_and_rays;
		
		//store correlation functions
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
		
		//miscellaneous
		ofstream * global_out_stream_ptr;
		int global_folderindex;
		string global_path;
		string global_runfolder;
		string global_resultsfolder_stem;
		string no_df_stem;

	public:
		inline int indexer(const int ipt, const int ipphi, const int iqt, const int iqx, const int iqy, const int iqz, const int itrig);
		inline int indexer2(const int iKT, const int iKPHI, const int iqt, const int iqx, const int iqy, const int iqz, const int itrig);
		inline double lin_int(double x_m_x1, double one_by_x2_m_x1, double f1, double f2);

		void replace_parentheses(std::string & tempstring);
		bool fexists(const char *filename);

		void Regulate_CF(int ipt, int iqt, int iqx, int iqy, int iqz, double * CF, double * projCF);
		void Regulate_CF_Hampel(int ipt, int iqx, int iqy, int iqz,
												double * pphi_CF_slice, double * pphi_CF_slice_term1, double * pphi_CF_slice_term2, double * pphi_CF_slice_term3);
		void Regulate_CF_Hampel_v2(int ipt, int iqx, int iqy, int iqz,
												double * pphi_CF_slice, double * pphi_CF_slice_term1, double * pphi_CF_slice_term2, double * pphi_CF_slice_term3);

		void Fill_out_pts(double * pointsarray, int numpoints, double max_val, int spacing_type);

		//miscellaneous
		void Set_path(string path);
		void Set_resultsfolder_stem(string usrdef_stem);
		void Set_runfolder(string runfolder);
		void Set_use_delta_f(bool usrdef_usedeltaf);
		double place_in_range(double phi, double min, double max);
		int lookup_resonance_idx_from_particle_id(int particle_id);

		void Average_total_target_eiqx_dN_dypTdpTdphi();
		void Set_correlation_function_q_pts();
		void Set_q_points();
		void Set_sorted_q_pts_list();
		void Get_q_points(double qo, double qs, double ql, double KT, double Kphi, double * qgridpts);
		void test_interpolator();
		void R2_Fourier_transform(int ipt, double plane_psi, int mode);
		double Extrapolate_Gaussian_1D(double q0, double qi0, double qi1, double f1, double f2);
		double Extrapolate_Gaussian_2D(double * q0, double * qi0, double * qi1, double (*vals) [2]);
		double Extrapolate_Gaussian_3D(double * q0, double * qi0, double * qi1, double (*vals) [2][2]);

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
		void get_CF(double * totalresult, double * thermalresult, double * crosstermresult, double * resonanceresult,
									int ipt, int ipphi, int iqt, int iqx, int iqy, int iqz, bool return_projected_value);
		void Compute_correlationfunction(double * totalresult, double * thermalresult, double * crosstermresult, double * resonanceresult,
										int ipt, int ipphi, int iqx, int iqy, int iqz, double qt_interp, int interp_flag = 0);
		void Cal_correlationfunction();
		void Fit_Correlationfunction3D(double *** Correl_3D, int ipt, int ipphi, bool fleshing_out_CF = true);
		int print_fit_state_3D (size_t iteration, gsl_multifit_fdfsolver * solver_ptr);
		void Fit_Correlationfunction3D_withlambda(double *** Correl_3D, int ipt, int ipphi, bool fleshing_out_CF = true);
		int print_fit_state_3D_withlambda (size_t iteration, gsl_multifit_fdfsolver * solver_ptr);
		//int Read_correlationfunction(int iKT, int iKphi);
		inline double get_fit_results(int i, gsl_multifit_fdfsolver * solver_ptr);
		inline double get_fit_err (int i, gsl_matrix * covariance_ptr);
		double gsl_polynomial_fit(const vector<double> &data_x, const vector<double> &data_y, const int order, double & chisq, bool verbose = false);
		double best_fit_rational_function(vector<double> & xdata, vector<double> & ydata, int n, int m, double x, bool & error_report);
		void Get_total_target_eiqx_dN_dypTdpTdphi_on_pair_momentum_grid(double * eiqx_EdNd3p_in, double * eiqx_EdNd3p_out, double * KT, double * KPHI);

		// input and output function prototypes
		void Readin_total_target_eiqx_dN_dypTdpTdphi(int folderindex);
		void Readin_total_target_eiqx_dN_dypTdpTdphi(string filename);
		void Output_total_target_eiqx_dN_dypTdpTdphi();
		void Readin_total_target_eiqx_dN_dypTdpTdphi_evavg();
		void Output_results(int mode);
		void Readin_results(int mode);
		void Read_in_all_dN_dypTdpTdphi();
		void Readin_resonance_fraction(int folderindex);
		void Output_correlationfunction();
		void Output_fleshed_out_correlationfunction(int ipt, int ipphi);
		void Dump_spectra_array(string output_filename, double *** array_to_dump);
		void Load_spectra_array(string output_filename, double *** array_to_read);
		void Output_lambdas();

		//parameters that the user is free to define
		double plumberg_test_variable;
		bool use_delta_f;
		bool append_output;
		int n_events;
		vector<int> osr;
		int initial_event, currentfolderindex;
		bool read_in_all_dN_dypTdpTdphi, output_all_dN_dypTdpTdphi;
		double fraction_of_resonances;
		double * SPinterp_pT_public;

		FitCF(particle_info* all_particles_in, int Nparticle, int particle_idx, vector<int> chosen_events, ofstream& myout,
					const int n_interp_pT_pts_in, const int n_interp_pphi_pts_in, const int qtnpts_in, const int qxnpts_in, const int qynpts_in, const int qznpts_in);
		~FitCF();

};

#endif
