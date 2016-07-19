#ifndef PARAMETERS_H
#define PARAMETERS_H

#include<string>
#include<sstream>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>// gsl random number generators
#include <gsl/gsl_randist.h>  // gsl random number distributions
#include <gsl/gsl_vector.h>   // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>     // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting

using namespace std;

#define VERBOSE 			1		// specifies level of output - 0 is lowest (no output)

const double hbarC=0.197327053;  //GeV*fm
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

//single particle spectra info
const int n_SP_pT = 15;
const int n_SP_pphi = 48;
const double SP_pT_min = 0.0;
const double SP_pT_max = 3.0;

//parameters for interpolation grid
const int n_interp_pT_pts = 15;
const int n_interp_pphi_pts = 48;
const double interp_pT_min = 0.0;
const double interp_pphi_min = 0.0;
const double interp_pT_max = 4.0;
const double interp_pphi_max = 2.*M_PI;
const double Del2_pT = (interp_pT_max - interp_pT_min) / (double)(n_interp_pT_pts-1);
const double Del2_pphi = (interp_pphi_max - interp_pphi_min) / (double)(n_interp_pphi_pts-1);

//pair momentum info
const int n_localp_T = 101;
const double localp_T_min = 0.01;
const double localp_T_max = 1.01;
const int n_localp_phi = 48;
const double localp_phi_min = 0.0;
const double localp_phi_max = 2.0 * M_PI;

#endif
