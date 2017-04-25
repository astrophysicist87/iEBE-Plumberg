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
#include <gsl/gsl_math.h>
#include <gsl/gsl_chebyshev.h>

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

// General information
const int ntrig = 2;			// for cos or sin
const double hbarC = 0.197327053;		//GeV*fm
const double hbarC3 = 0.00768351405;
const double hbarCm1 = 5.067728853;
const double twopi = 2.0*M_PI;
const double MeVToGeV = 0.001;

// Particle information
const int Maxparticle=400;            //size of array for storage of the particles
const int Maxdecaychannel=13;
const int Maxdecaypart=5;

// Relative momentum information
const double delta_qt = 0.00625;
const double delta_qx = 0.0125;
const double delta_qy = 0.0125;
const double delta_qz = 0.015;
const int new_nqpts = 51;	//for fleshing out

// Single particle spectra info
const double SP_pT_min = 0.0;
const double SP_pphi_min = 0.0;
const double SP_Del_pY_min = 0.0;
const double SP_pT_max = 4.0;
const double SP_pphi_max = 2.0*M_PI;
const double SP_Del_pY_max = 4.0;

// Pair momentum info
const double Kphi_min = 0.0;
const double Kphi_max = 2.0*M_PI;

// Phase-space integral information
const int n_zeta_pts = 12;
const int n_v_pts = 12;
const int n_s_pts = 12;

// Fitting information
const size_t fit_max_iterations = 1000;  // stop at this point if not converged 
const double fit_tolerance = 1e-6;

#endif
