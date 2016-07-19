#include <cmath>
#include <iostream>
#include <stdlib.h>

#include "chebyshev_library.h"
#include "chebyshev.h"

using namespace std;
using namespace csf;

Chebyshev::Chebyshev(double * fpts_in, int * numbers_of_points_in, int * orders_in, double * lower_limits_in, double * upper_limits_in, int dimension_in)
{
	dimension = dimension_in;
	modes = new int [dimension];
	for (int idim = 0; idim < dimension; ++idim)
		modes[idim] = 0;	//by default, assume standard Chebyshev interpolation in every direction

	numbers_of_points = new int [dimension];
	orders = new int [dimension];
	lower_limits = new double [dimension];
	upper_limits = new double [dimension];

	for (int idim = 0; idim < dimension; ++idim)
	{
		numbers_of_points[idim] = numbers_of_points_in[idim];
		orders[idim] = orders_in[idim];
		if (1+orders[idim] != numbers_of_points[idim])	//only order of interpolation currently supported
		{
			cerr << "Warning: falling back to orders[" << idim << "] = numbers_of_points[" << idim << "] - 1 = " << numbers_of_points[idim] - 1 << endl;
			orders[idim] = numbers_of_points[idim] - 1;
		}
		lower_limits[idim] = lower_limits_in[idim];
		upper_limits[idim] = upper_limits_in[idim];
	}

	set_total_coeffs_length();
	set_total_fpts_length();

	fpts = new double [total_fpts_length];
	for (int ifpt = 0; ifpt < total_fpts_length; ++ifpt)
		fpts[ifpt] = fpts_in[ifpt];

	coeffs = new double [total_coeffs_length];
	for (int ic = 0; ic < total_coeffs_length; ++ic)
		coeffs[ic] = 0.0;

	//once everything is set, go ahead and do appropriate calculations so that
	//approximating polynomial is ready to evaluate
	get_Chebyshev_coefficients();

	return;
}

Chebyshev::Chebyshev(double * fpts_in, int * numbers_of_points_in, int * orders_in, double * lower_limits_in, double * upper_limits_in, int dimension_in, int * modes_in)
{
	dimension = dimension_in;
	modes = new int [dimension];
	for (int idim = 0; idim < dimension; ++idim)
		modes[idim] = modes_in[idim];

	numbers_of_points = new int [dimension];
	orders = new int [dimension];
	lower_limits = new double [dimension];
	upper_limits = new double [dimension];

	for (int idim = 0; idim < dimension; ++idim)
	{
		numbers_of_points[idim] = numbers_of_points_in[idim];
		orders[idim] = orders_in[idim];
		if (1+orders[idim] != numbers_of_points[idim])	//only order of interpolation currently supported
		{
			cerr << "Warning: falling back to orders[" << idim << "] = numbers_of_points[" << idim << "] - 1 = " << numbers_of_points[idim] - 1 << endl;
			orders[idim] = numbers_of_points[idim] - 1;
		}
		lower_limits[idim] = lower_limits_in[idim];
		upper_limits[idim] = upper_limits_in[idim];
	}

	set_total_coeffs_length();
	set_total_fpts_length();

	fpts = new double [total_fpts_length];
	for (int ifpt = 0; ifpt < total_fpts_length; ++ifpt)
		fpts[ifpt] = fpts_in[ifpt];

	coeffs = new double [total_coeffs_length];
	for (int ic = 0; ic < total_coeffs_length; ++ic)
		coeffs[ic] = 0.0;

	//once everything is set, go ahead and do appropriate calculations so that
	//approximating polynomial is ready to evaluate
	get_Chebyshev_coefficients();

	return;
}

Chebyshev::~Chebyshev()
{
	delete [] fpts;
	delete [] coeffs;
	delete [] numbers_of_points;
	delete [] orders;
	delete [] lower_limits;
	delete [] upper_limits;

	for (int ic = 0; ic < total_coeffs_length; ++ic)
		delete [] coeffs_indices[ic];

	for (int ifpt = 0; ifpt < total_fpts_length; ++ifpt)
		delete [] fpts_indices[ifpt];

	delete [] coeffs_indices;
	delete [] fpts_indices;

	return;
}

inline double Chebyshev::dot(double * x, double * y, int length)
{
    double sum = 0.0;
    for (int i = 0; i < length; ++i)
        sum += x[i]*y[i];
    return (sum);
}

inline double semi_infinite_scaling(double x, double L)
{
	return ( L * ( 1. + x ) / ( 1. - x ) );
}

inline double semi_infinite_inverse_scaling(double x, double L)
{
	return ( ( x - L ) / ( x + L ) );
}

void Chebyshev::set_nodes(int number_of_points, double * nodes)
{
    for (int i = 1; i <= number_of_points; ++i)
        nodes[i-1] = -cos(M_PI * (2.*i-1.) / (2.*number_of_points));

    return;
}

void Chebyshev::get_nodes(double a, double b, int number_of_points, double * nodes, double * adjnodes)
{
    double halfwidth = 0.5*(b-a);

    set_nodes(number_of_points, nodes);

    for (int i = 1; i <= number_of_points; ++i)
       adjnodes[i-1] =  (nodes[i-1] + 1.0)*halfwidth + a;

    return;
}

void Chebyshev::get_Chebyshev_points(int number_of_points, int order, double * nodes, double ** Tpts)
{
    for (int j = 0; j <= order; ++j)
    for (int i = 1; i <= number_of_points; ++i)
        Tpts[j][i-1] = csf::Tfun(j, nodes[i-1]);

    return;
}

void Chebyshev::get_indices(int n, int total_dims_size, int * dims, int * indices, int lookup_index)
{
    // let dims be {11,5,5,7}
    // then lookup_index == i4 + 7*i3 + 7*5*i2 + 7*5*5*i1
    int tmp = total_dims_size;
    int resid = lookup_index;
    for (int id = 0; id < n; ++id)
    {
        tmp /= dims[id];
        indices[id] = resid / tmp;
        resid = resid % tmp;
    }
    return;
}

void Chebyshev::set_coeffs_indices(int total_coeffs_length)
{
    coeffs_indices = new int * [total_coeffs_length];
	for (int ic = 0; ic < total_coeffs_length; ++ic)
    {
        coeffs_indices[ic] = new int [dimension];
        get_indices(dimension, total_coeffs_length, numbers_of_points, coeffs_indices[ic], ic);
    }
    return;
}

void Chebyshev::set_ifpt_indices(int total_fpts_length)
{
    fpts_indices = new int * [total_fpts_length];
    for (int ifpt = 0; ifpt < total_fpts_length; ++ifpt)
    {
        fpts_indices[ifpt] = new int [dimension];
        get_indices(dimension, total_fpts_length, numbers_of_points, fpts_indices[ifpt], ifpt);
    }
    return;
}

void Chebyshev::set_total_coeffs_length()
{
    total_coeffs_length = 1;
    for (int idim = 0; idim < dimension; ++idim)
        total_coeffs_length *= (1 + orders[idim]);
	return;
}

void Chebyshev::set_total_fpts_length()
{
    total_fpts_length = 1;
    for (int idim = 0; idim < dimension; ++idim)
        total_fpts_length *= numbers_of_points[idim];
	return;
}

/*USAGE: debugger(__LINE__, __FILE__)*/
void inline debugger(int cln, const char* cfn)
{
	cerr << "You made it to " << cfn << ":" << cln << "!" << endl;
	return;
}

void Chebyshev::get_Chebyshev_coefficients()
{
    double ** nodes = new double * [dimension];		//x_k, k=1..m
    double ** adjnodes = new double * [dimension];	//z_k
    double *** Tpts = new double ** [dimension];

	//for each dimension, set number of nodes and order of interpolation in that dimension
    for (int idim = 0; idim < dimension; ++idim)
    {
        int current_number_of_points = numbers_of_points[idim];
        int current_order = orders[idim];
        nodes[idim] = new double [current_number_of_points];
        adjnodes[idim] = new double [current_number_of_points];
    
        // this is for the Chebyshev points
        Tpts[idim] = new double * [current_order+1];
        for (int io = 0; io <= current_order; ++io)
            Tpts[idim][io] = new double [current_number_of_points];

        get_nodes(lower_limits[idim], upper_limits[idim], current_number_of_points, nodes[idim], adjnodes[idim]);
    
        get_Chebyshev_points(current_number_of_points, current_order, nodes[idim], Tpts[idim]);
    }

    set_coeffs_indices(total_coeffs_length);
    set_ifpt_indices(total_fpts_length);

	//compute the coefficients (a_i) here...
    for (int ic = 0; ic < total_coeffs_length; ++ic)
    {
        double num = 0.0;
        for (int ifpt = 0; ifpt < total_fpts_length; ++ifpt)
        {
            double numprod = 1.0;
            for (int idim = 0; idim < dimension; ++idim)
                numprod *= Tpts[idim][coeffs_indices[ic][idim]][fpts_indices[ifpt][idim]];
            num += fpts[ifpt] * numprod;
        }
        double den = 1.0;
        for (int idim = 0; idim < dimension; ++idim)
            den *= dot( Tpts[idim][coeffs_indices[ic][idim]], Tpts[idim][coeffs_indices[ic][idim]], numbers_of_points[idim]);
        coeffs[ic] = num / den;
    }

    for (int idim = 0; idim < dimension; ++idim)
	{
	    delete [] nodes[idim];
	    delete [] adjnodes[idim];
	    for (int io = 0; io <= orders[idim]; ++io)
	        delete [] Tpts[idim][io];
	    delete [] Tpts[idim];
	}
    delete [] nodes;
    delete [] adjnodes;
    delete [] Tpts;

    return;
}

double Chebyshev::eval(double * p)
{

	//p is point at which to evaluate approximating polynomial
	double * unadj_p = new double [dimension];	//the unadjusted point
	for (int idim = 0; idim < dimension; ++idim)
	{
		switch (modes[idim])
		{
			case 0:		//standard Chebyshev interpolation on finite interval
				unadj_p[idim] = 2.0 * ((p[idim] - lower_limits[idim]) / (upper_limits[idim] - lower_limits[idim])) - 1.0;
				break;
			case 1:		//semi-infinite Chebyshev interpolation: scale set by upper_limits[idim]
				unadj_p[idim] = semi_infinite_inverse_scaling(p[idim], upper_limits[idim]);
				break;
			default:
				cerr << "Chebyshev(): invalid mode choice of modes[dim = " << idim << "] == " << modes[idim] << "!" << endl;
				exit(1);
				break;
		}
	}

	double result = 0.0;
	for (int ic = 0; ic < total_coeffs_length; ++ic)
	{
		double prod = 1.0;
		for (int idim = 0; idim < dimension; ++idim)
			prod *= csf::Tfun(coeffs_indices[ic][idim], unadj_p[idim]);
		result += coeffs[ic] * prod;
	}

	delete [] unadj_p;
	return (result);
}

void Chebyshev::eval(double ** points, int npoints, double * results)
{
	//double ** points (npoints x dimension) - points at which to evaluate approximating polynomial
	//int npoints - obvious
	//double * results - approximated results stored here
	for (int ip = 0; ip < npoints; ++ip)
		results[ip] = eval(points[ip]);
	return;
}

//End of file
