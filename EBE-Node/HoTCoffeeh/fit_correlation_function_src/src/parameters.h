#ifndef PARAMETERS_H
#define PARAMETERS_H

#include<string>
#include<sstream>

using namespace std;

#define USE_OLD_INTERP					true		// duh
#define VERBOSE 						3			// specifies level of output - 0 is lowest (no output)
#define DEBUG							false		// flag for output of debugging statements
#define USE_LAMBDA						true		// fit correlation function with adjustable intercept parameter
#define IGNORE_LONG_LIVED_RESONANCES	true		// particularly, whether to include eta or eta' in spectra calculations
													// true means C(q=0) ~ 1 + \lambda
#define QT_POINTS_SPACING				1			// 0 - uniform from -qmax to +qmax
													// 1 - Chebyshev nodes from -qmax to +qmax
#define QX_POINTS_SPACING				0
#define QY_POINTS_SPACING				0
#define QZ_POINTS_SPACING				0
#define Q_AXES_AND_RAYS_ONLY			false		// true - only do points along q-axes (only works for odd points right now)
													// false - do full grid
#define FIT_WITH_PROJECTED_CFVALS		true		// as opposed to unprojected CFvals...
#define FLESH_OUT_CF					true		// refines grid via interpolation before fitting
#define REGULATE_CF						false		// true (false) means (don't) try to catch spurious values of projected
													// or regular CF and replace them with median value in that window
#define SLICE_OF_FLESH_ONLY				true		// full correlation function (fleshed out) is typically a HUGE file (~10GB),
													// so use this to output q-slices only (MUCH) smaller
#define THERMAL_ONLY					false		// duh



const int ntrig = 2;			// for cos or sin

const double hbarC = 0.197327053;		//GeV*fm
const double hbarC3 = 0.00768351405;
const double hbarCm1 = 5.067728853;
const double twopi = 2.*M_PI;
const double MeVToGeV = 0.001;

//particle information
const int Maxparticle=400;            //size of array for storage of the particles
const int Maxdecaychannel=13;
const int Maxdecaypart=5;

//spatial rapidity information
const int eta_s_npts = 15;
const double eta_s_i = 0.0;
const double eta_s_f = 4.0;

//extrapolation information
const int polynomial_fit_order = 4;
const int rational_function_numerator_order = 3;
const int rational_function_denominator_order = 4;
const int UDPMsize = 15;
const int UDPMTsize = 10;
static double usr_def_pc_markers[UDPMsize] = {
					0.00, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82,
					0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90
				};
static double usr_def_pc_markers_thinned[UDPMTsize] = { 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90 };

//phase space integral info
const int s_npts = 12;
const int v_npts = 12;
const int zeta_npts = 12;

//relative momentum information
const int qonpts = 11;
const int qsnpts = 11;
const int qlnpts = 11;
const int qnpts = 1;
const double delta_q = 0.005;
const double init_q = 0.0;

const int new_nqpts = 51;
const int new_nqxpts = new_nqpts;
const int new_nqypts = new_nqpts;
const int new_nqzpts = new_nqpts;

//all direction-specific q points information here
//const int qtnpts = 13;
//const int qxnpts = 7;
//const int qynpts = 7;
//const int qznpts = 7;
//try to make max. sqrt(q dot q) ~ 0.025 GeV or so
const double delta_qt = 0.00625;
const double delta_qx = 0.025;
const double delta_qy = 0.025;
const double delta_qz = 0.0125;
//const double delta_qx = 0.003;
//const double delta_qy = 0.003;
//const double delta_qz = 0.003;
//const double init_qt = -0.5*double(qtnpts-1)*delta_qt;
//const double init_qx = -0.5*double(qxnpts-1)*delta_qx;
//const double init_qy = -0.5*double(qynpts-1)*delta_qy;
//const double init_qz = -0.5*double(qznpts-1)*delta_qz;

//single particle spectra info
const int n_SP_pT = 15;
const int n_SP_pphi = 48;
const double SP_pT_min = 0.0;
const double SP_pT_max = 3.0;

//parameters for interpolation grid
//  - polar
//const int n_interp_pT_pts = 15;
//const int n_interp_pphi_pts = 12;
const double interp_pT_min = 0.0;
const double interp_pphi_min = 0.0;
const double interp_pT_max = 4.0;
const double interp_pphi_max = 2.*M_PI;
//const double Del2_pT = (interp_pT_max - interp_pT_min) / (double)(n_interp_pT_pts-1);
//const double Del2_pphi = (interp_pphi_max - interp_pphi_min) / (double)(n_interp_pphi_pts-1);
const double Del_pT = 0.01;
const double Del2_pphi = 0.01;

//correlation function info
const int corrfuncdim = 1;
const bool lambdaflag = USE_LAMBDA;
const double correlator_minus_one_cutoff = 0.0;		//zero means all calculations happen as usual

//pair momentum info
const int n_localp_T = 101;
const double localp_T_min = 0.01;
const double localp_T_max = 1.01;
const int n_localp_phi = 48;
const double localp_phi_min = 0.0;
const double localp_phi_max = 2*M_PI;

const int n_order = 1;

const size_t fit_max_iterations = 1000;  // stop at this point if not converged 
const double fit_tolerance = 1e-6;

//misc. resonance info
const double max_lifetime = 100.;	// fm/c

namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}

#endif
