#ifndef HBTEG_H
#define HBTEG_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>

#include "ParameterReader.h"
#include "Arsenal.h"
#include "EventRecord.h"
#include "ParticleRecord.h"

using namespace std;

class HBT_event_generator
{
	private:
		ParameterReader * paraRdr;

		//header info
		int n_pT_pts, n_pphi_pts, n_pY_pts;
		int n_KT_pts, n_Kphi_pts, n_KL_pts;

		int n_pT_bins, n_pphi_bins, n_pY_bins;
		int n_KT_bins, n_Kphi_bins, n_KL_bins;

		int n_qo_pts, n_qs_pts, n_ql_pts;
		int n_qo_bins, n_qs_bins, n_ql_bins;

		double pT_min, pT_max, pphi_min, pphi_max, pY_min, pY_max;
		double KT_min, KT_max, Kphi_min, Kphi_max, KL_min, KL_max;

		double init_qo, init_qs, init_ql;
		double delta_qo, delta_qs, delta_ql;

		double pT_bin_width, pphi_bin_width, pY_bin_width;
		double KT_bin_width, Kphi_bin_width, KL_bin_width;

		vector<string> all_file_names;
		vector<EventRecord> allEvents;

		vector<double> pT_pts, pphi_pts, pY_pts;
		vector<double> KT_pts, Kphi_pts, KL_pts;

		vector<double> dN_pTdpTdpphidpY;
		
		/*
		//store correlation functions
		//double *** Correl_3D;
		double ***** CFvals, ***** thermalCFvals, ***** crosstermCFvals, ***** resonancesCFvals;
		double *** fleshed_out_CF, *** fleshed_out_thermal, *** fleshed_out_crossterm, *** fleshed_out_resonances;
		double *** Correl_3D_err;
		double ** lambda_Correl, ** lambda_Correl_err;
		double ** lambda_QM;
		int *** correlator_minus_one_cutoff_norms;
		double * qx_fleshed_out_pts, * qy_fleshed_out_pts, * qz_fleshed_out_pts;

		//HBT radii coefficients
		double ** R2_side_GF, ** R2_out_GF, ** R2_long_GF, ** R2_outside_GF, ** R2_sidelong_GF, ** R2_outlong_GF;
		double ** R2_side_GF_C, ** R2_out_GF_C, ** R2_long_GF_C, ** R2_outside_GF_C, ** R2_sidelong_GF_C, ** R2_outlong_GF_C;
		double ** R2_side_GF_S, ** R2_out_GF_S, ** R2_long_GF_S, ** R2_outside_GF_S, ** R2_sidelong_GF_S, ** R2_outlong_GF_S;
		double ** R2_side_err, ** R2_out_err, ** R2_long_err, ** R2_outside_err, ** R2_sidelong_err, ** R2_outlong_err;

		double ** R2_side_QM, ** R2_out_QM, ** R2_long_QM, ** R2_outside_QM, ** R2_sidelong_QM, ** R2_outlong_QM;
		double ** R2_side_QM_C, ** R2_out_QM_C, ** R2_long_QM_C, ** R2_outside_QM_C, ** R2_sidelong_QM_C, ** R2_outlong_QM_C;
		double ** R2_side_QM_S, ** R2_out_QM_S, ** R2_long_QM_S, ** R2_outside_QM_S, ** R2_sidelong_QM_S, ** R2_outlong_QM_S;
		*/
		
		//miscellaneous
		string path;
		ostream & out;
		ostream & err;

		/*inline int bin( double datapoint, double * bin_limits );

		double bin_function( vector<double> mom1,
								vector<double> mom2,
								int mode,
								void * bin_function_parameters );
		double bin_function_mode_1( vector<double> mom1,
								vector<double> mom2,
								void * bin_function_parameters );*/

		//double bin_function_mode_2( vector<double> mom1,
		//						vector<double> mom2,
		//						void * bin_function_parameters );


		vector<double> numerator;



	public:

		// Constructors, destructors, and initializers
		HBT_event_generator( ParameterReader * paraRdr_in,
								const vector<EventRecord> & allEvents_in,
								ostream & out_stream = std::cout,
								ostream & err_stream = std::cerr )
								:
								out(out_stream),
								err(err_stream)
								{ initialize_all( paraRdr_in, allEvents_in ); };

		void initialize_all(ParameterReader * paraRdr_in,
								const vector<EventRecord> & allEvents_in);

		~HBT_event_generator();

		// Library functions
		//inline int bin_function( double datapoint, const vector<double> & points );
		//inline int indexer(int ipT, int ipphi, int ipY);

		//inline double get_q0(double m, double qo, double qs, double ql, double KT, double KL);

		inline int bin_function( double datapoint, const vector<double> & points )
		{
			//out << "Here..." << endl;

			int result = (int)( ( datapoint - points[0] )
							* double( points.size()-1 )
							/ ( points[points.size()-1] - points[0] ) );

			//out << "...to here." << endl;

			// Assume uniform bin-widths for now
			return ( result );
		}

		inline int indexer(int ipT, int ipphi, int ipY)
		{
			return ( ( ipT * n_pphi_bins + ipphi ) * 1 + ipY );
		}


		inline double get_q0(double m, double qo, double qs, double ql, double KT, double KL)
		{
			double xi2 = m*m + KT*KT + KL*KL + 0.25*(qo*qo + qs*qs + ql*ql);

			return ( sqrt(xi2 + qo*KT + ql*KL) - sqrt(xi2 - qo*KT - ql*KL) );
		}


		// Functions to compute single-particle spectra
		void Compute_spectra();
		// Sub-methods
		void Compute_dN_pTdpTdpphidpY();
		//
		void Compute_dN_pTdpTdpphi();
		void Compute_dN_2pipTdpTdpY();
		void Compute_dN_dpphidpY();
		//
		void Compute_dN_2pipTdpT();
		void Compute_dN_dpphi();
		void Compute_dN_2pidpY();
		// total multiplicity
		//void Compute_N();







		//void Compute_correlation_function();

		//void Compute_numerator();
		//void Compute_denominator();








		// Input/output
		//void Output_correlation_function();

		//std::ostream out, err;

};

#endif
