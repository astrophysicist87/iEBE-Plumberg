#ifndef PARAMETERS_H
#define PARAMETERS_H

#include<string>
#include<sstream>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

using namespace std;

#define VERBOSE 			1		// specifies level of output - 0 is lowest (no output)

// General information
const double hbarC=0.197327053;  	// GeV*fm
const double twopi = 2.*M_PI;
const double MeVToGeV = 0.001;

// Particle information
const int Maxparticle=400;			// size of array for storage of the particles
const int Maxdecaychannel=13;
const int Maxdecaypart=5;

// Spatial rapidity information
const int eta_s_npts = 15;
const double eta_s_i = 0.0;
const double eta_s_f = 4.0;

// Single particle spectra info
const double SP_pT_min = 0.0;
const double SP_pphi_min = 0.0;
const double SP_pT_max = 4.0;
const double SP_pphi_max = 2.*M_PI;

// Pair momentum info
const double Kphi_min = 0.0;
const double Kphi_max = 2*M_PI;

// Phase-space integral information
const int n_zeta_pts = 12;
const int n_v_pts = 12;
const int n_s_pts = 12;

#endif
