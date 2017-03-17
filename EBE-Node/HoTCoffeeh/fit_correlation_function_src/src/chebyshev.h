#ifndef CHEBYSHEV_H
#define CHEBYSHEV_H

#include <cmath>
#include "chebyshev_library.h"

using namespace std;

class Chebyshev
{
	private:
		int * modes;
		//double ** Tpts;
        double * lower_limits, * upper_limits, * coeffs, * fpts;
        int dimension, total_coeffs_length, total_fpts_length, * numbers_of_points, * orders;
        //double ** nodes, ** adjusted_nodes;
        int ** coeffs_indices, ** fpts_indices;

		inline double dot(double * x, double * y, int length);
		void set_nodes(int number_of_points, double * nodes);
		//void get_nodes(double a, double b, int number_of_points, double * nodes, double * adjnodes);
		void get_Chebyshev_points(int number_of_points, int order, double * nodes, double ** Tpts);
		void set_total_coeffs_length();
		void set_total_fpts_length();
		void get_indices(int n, int total_dims_size, int * dims, int * indices, int lookup_index);
		void set_coeffs_indices(int total_coeffs_length);
		void set_ifpt_indices(int total_fpts_length);
		void get_Chebyshev_coefficients();

	public:
		Chebyshev(double * fpts_in, int * numbers_of_points_in, int * orders_in, double * lower_limits_in, double * upper_limits_in, int dimension_in);
		Chebyshev(double * fpts_in, int * numbers_of_points_in, int * orders_in, double * lower_limits_in, double * upper_limits_in, int dimension_in, int * modes_in);
		~Chebyshev();

		void get_nodes(double a, double b, int number_of_points, double * nodes, double * adjnodes);

		double eval(double * p);
		void eval(double ** points, int npoints, double * results);
};

#endif
//End of file
