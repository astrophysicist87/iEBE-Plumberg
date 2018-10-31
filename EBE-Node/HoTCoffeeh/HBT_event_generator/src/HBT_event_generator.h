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
#include "EventRecord.h"
#include "ParticleRecord.h"

using namespace std;

class HBT_event_generator
{
	private:
		ParameterReader * paraRdr;

		//header info
		int n_pT_pts, n_pphi_pts, n_pY_pts, nKT, nKphi;
		int qonpts, qsnpts, qlnpts;
		int qonbins, qsnbins, qlnbins;
		double KT_min, KT_max;
		double init_qo, init_qs, init_ql;
		double delta_qo, delta_qs, delta_ql;

		vector<string> all_file_names;
		
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
		ostream out, err;
		string path;

		inline int bin( double datapoint, double * bin_limits );

		double bin_function( vector<double> mom1,
								vector<double> mom2,
								int mode,
								void * bin_function_parameters );
		double bin_function_mode_1( vector<double> mom1,
								vector<double> mom2,
								void * bin_function_parameters );

		//double bin_function_mode_2( vector<double> mom1,
		//						vector<double> mom2,
		//						void * bin_function_parameters );


	public:

		HBT_event_generator( ParameterReader * paraRdr_in,
								vector<EventRecord> &allEvents,
								ostream& out_stream = std::cout,
								ostream& err_stream = std::cerr );

		~HBT_event_generator();

		// Main procedure
		void Compute_correlation_function();

		// Auxiliary procedures
		void Compute_numerator();
		void Compute_denominator();

		// Input/output
		void Output_correlation_function();

};

#endif
