#ifndef EMISS_FN_FIT_H
#define EMISS_FN_FIT_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>
#include <string>
#include <time.h>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <complex>
#include <limits>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>            // gsl random number generators
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting

using namespace std;

typedef struct
{
   double rval, phival, etaval, tauval;
   double Emiss_fn, error;
}emiss_fcn_data;

struct data
{
	size_t my_n;
	double * y;
	double * sigma;
};

struct phys_params
{
	double M_perp, T_0, eta_0, Y_rapidity;
	double Rad, eta_f, tau_f, K_perp, Phi_K;
	double Delta_tau, Delta_eta, S0;

	double v_3_bar, eps_3_bar, psi_3_bar;
};
vector<double>* emiss_fn_fit(double KT, double Kphi, double KY, vector<emiss_fcn_data>* passed_ptr);

#endif
