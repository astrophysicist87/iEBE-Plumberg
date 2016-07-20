#ifndef PARAMETERS_H
#define PARAMETERS_H

#include<string>
#include<sstream>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>            // gsl random number generators
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>

#include "H5Cpp.h"

using namespace std;

#define VERBOSE 			3		// specifies level of output - 0 is lowest (no output)
#define QT_POINTS_SPACING		1		// 0 - uniform from -qmax to +qmax
							// 1 - Chebyshev nodes from -qmax to +qmax
#define QX_POINTS_SPACING		0		// same
#define QY_POINTS_SPACING		0		// same
#define QZ_POINTS_SPACING		0		// same

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

//HDF information
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
const int UDPMsize = 15;
static double usr_def_pc_markers[UDPMsize] = {
					0.00, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82,
					0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90
				};

//phase space integral info
const int s_npts = 12;
const int v_npts = 12;
const int zeta_npts = 12;

//relative momentum information
const double delta_qt = 0.00625;
const double delta_qx = 0.0125;
const double delta_qy = 0.0125;
const double delta_qz = 0.015;
const int new_nqpts = 51;	//for fleshing out

//single particle spectra info
const int n_SP_pT = 15;
const int n_SP_pphi = 48;
const double SP_pT_min = 0.0;
const double SP_pT_max = 3.0;

//parameters for interpolation grid
const double interp_pT_min = 0.0;
const double interp_pphi_min = 0.0;
const double interp_pT_max = 4.0;
const double interp_pphi_max = 2.*M_PI;

//pair momentum info
const double Kphi_min = 0.0;
const double Kphi_max = 2*M_PI;

const size_t fit_max_iterations = 1000;  // stop at this point if not converged 
const double fit_tolerance = 1e-6;

#endif
