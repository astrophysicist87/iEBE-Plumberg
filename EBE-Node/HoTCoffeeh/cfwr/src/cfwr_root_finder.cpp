#include<cmath>
#include<vector>
#include<cstdlib>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "cfwr.h"
#include "gauss_quadrature.h"
#include "bessel.h"

using namespace std;

struct chebyshev_params
{
	double * coefficients;
};

gsl_cheb_series *cs_accel_expFunc;

inline double chebyshev_evaluator(double pY_local, void * params)
{
	struct chebyshev_params *p = (struct chebyshev_params *) params;
	cs_accel_expFunc->c = p->coefficients;
	return ( exp(-pY_local) * gsl_cheb_eval (cs_accel_expFunc, pY_local) );
}

double CorrelationFunction::root_finder(double * chebyshev_cfs, double a, double b)
{
	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r = 0;
	double x_lo = estimate_pY_shift - 0.1, x_hi = estimate_pY_shift + 0.1;	//try this for now
	gsl_function F;
	struct chebyshev_params params;
	params.coefficients = chebyshev_cfs;
	cs_accel_expFunc->a = a;
	cs_accel_expFunc->b = b;

	F.function = &chebyshev_evaluator;
	F.params = &params;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set (s, &F, x_lo, x_hi);

	do
	{
		iter++;
		status = gsl_root_fsolver_iterate (s);
		r = gsl_root_fsolver_root (s);
		x_lo = gsl_root_fsolver_x_lower (s);
		x_hi = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);
	}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free (s);

	return (r);
}

//End of file
