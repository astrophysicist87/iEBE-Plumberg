#ifndef CHEBYSHEV_LIBRARY_H
#define CHEBYSHEV_LIBRARY_H

using namespace std;

namespace csf
{
	// n = 0
	inline double T0(double x)
	{
		return 1.0;
	}

	// n = 1
	inline double T1(double x)
	{
		return x;
	}

	// n = 2
	inline double T2(double x)
	{
		return (2.0 * x*x) - 1.0;
	}

	/*
	 *	Tn(x)
	 */
	inline double Tfun(unsigned int n, double x)
	{
		if (n == 0)
		{
			return T0(x);
		}
		else if (n == 1)
		{
			return T1(x);
		}
		else if (n == 2)
		{
			return T2(x);
		}

		/* We could simply do this:
		    return (2.0 * x * Tn(n - 1, x)) - Tn(n - 2, x) ;
		   but it could be slow for large n */
 
		double tnm1(T2(x));
		double tnm2(T1(x));
		double tn(tnm1);
	
	    for (unsigned int l = 3; l <= n; l++)
	    { 
			tn = (2.0 * x * tnm1) - tnm2;
			tnm2 = tnm1;
			tnm1 = tn;
	    }

		return tn;
	}
}

#endif

//End of file
