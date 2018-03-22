#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>

#include "cfwr.h"
#include "cfwr_lib.h"
#include "Stopwatch.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"
#include "chebyshev.h"

using namespace std;

void CorrelationFunction::Get_GF_HBTradii()
{
	*global_out_stream_ptr << "--> Getting HBT radii by Gaussian fit method" << endl;

	if (FLESH_OUT_CF)
		Allocate_fleshed_out_CF();

	for (int ipt = 0; ipt < n_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	{
		//if (ipt > 0 or ipphi > 0)
		//	continue;
		*global_out_stream_ptr << "   --> Doing pT = " << SP_pT[ipt] << ", pphi = " << SP_pphi[ipphi] << "..." << endl;
		
		//determine whether to use fleshed out / projected CFvals
		double *** CF_for_fitting = CFvals[ipt][ipphi];
		if (FLESH_OUT_CF)
		{
			Flesh_out_CF(ipt, ipphi);
			CF_for_fitting = fleshed_out_CF;
		}

		//finally, do fits
		//Fit_Correlationfunction3D_withlambda( CF_for_fitting, ipt, ipphi, FLESH_OUT_CF );
		find_minimum_chisq_correlationfunction_full( CF_for_fitting, ipt, ipphi, FLESH_OUT_CF );
	}

	if (FLESH_OUT_CF)
		Delete_fleshed_out_CF();

	return;
}

//to save time, don't bother making grid large enough to interpolate to OSL
//this function leaves open the option of just interpolating over the qt-direction
void CorrelationFunction::Cal_correlationfunction(bool project_CF /*==true*/)
{
	string projection_stem = ( project_CF ) ? "projected" : "unprojected" ;
	*global_out_stream_ptr << "Calculating the " << projection_stem << " correlation function..." << endl;

	//exploits convenient symmetries to get full correlation function
	//as need for remainder of calculation
	Reflect_in_qz_and_qt();

	// Can't interpolate if there's only one point in qt-direction!
	if (qtnpts == 1)
		return;

	// chooses the qo, qs, ql (or qx, qy, ql) points at which to evaluate correlation function,
	// and allocates the array to hold correlation function values
	Set_correlation_function_q_pts();
	Allocate_CFvals();

	double * q_interp = new double [4];

	// Then compute full correlation function
	*global_out_stream_ptr << "Computing the full, " << projection_stem << " correlator in XYZ coordinates..." << endl;
	Stopwatch sw;
	sw.Start();
	for (int ipt = 0; ipt < n_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	{
//if (ipt > 0 || ipphi > 0) continue;
		Get_q_points(qx_pts[iqx], qy_pts[iqy], qz_pts[iqz], SP_pT[ipt], SP_pphi[ipphi], q_interp);
		//returns only projected value automatically if appropriate options are specified!
		double tmp1 = 0.0, tmp2 = 0.0, tmp2a = 0.0, tmp3 = 0.0;
		Compute_correlationfunction(&tmp1, &tmp2, &tmp2a, &tmp3, ipt, ipphi, iqx, iqy, iqz, q_interp[0], 0, project_CF);
		CFvals[ipt][ipphi][iqx][iqy][iqz] = 1.0 + tmp1;		//C == Ct + Cct + Cr + 1
		thermalCFvals[ipt][ipphi][iqx][iqy][iqz] = tmp2;	//Ct
		crosstermCFvals[ipt][ipphi][iqx][iqy][iqz] = tmp2a;	//Cct
		resonancesCFvals[ipt][ipphi][iqx][iqy][iqz] = tmp3;	//Cr
	}
	sw.Stop();
	*global_out_stream_ptr << "Finished computing " << projection_stem << " correlator in " << sw.printTime() << " seconds." << endl;

	//output the un-regulated correlation function to separate file for debugging purposes
	Output_correlationfunction(project_CF);

	//added for checking...
	if (FLESH_OUT_CF)
	{
		Allocate_fleshed_out_CF();

		for (int ipt = 0; ipt < n_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			//if (ipt > 0 or ipphi > 0)
			//	continue;
			Flesh_out_CF(ipt, ipphi);
			Output_fleshed_out_correlationfunction(ipt, ipphi, project_CF);
		}
	
		Delete_fleshed_out_CF();
	}

	delete [] q_interp;

	return;
}

void CorrelationFunction::Compute_correlationfunction(double * totalresult, double * thermalresult, double * CTresult, double * resonanceresult,
										int ipt, int ipphi, int iqx, int iqy, int iqz, double qt_interp, int interp_flag /*==0*/, bool project_CF /*==true*/)
{
	int qidx = binarySearch(qt_pts, qtnpts, qt_interp);
	double q_min = -qtmax;	//the range in qt expected to be covered by
	double q_max = qtmax;	//Chebyshev interpolation (for any qt spacing choice)

	bool q_point_is_outside_grid = ( ( qidx < 0 || qidx >= qtnpts ) && ( qt_interp < q_min || qt_interp > q_max ) );

	if (qxnpts==1 && qynpts==1 && qznpts==1)
	{
		double tmpC = 0.0, tmpCt = 0.0, tmpCct = 0.0, tmpCr = 0.0;
		*global_out_stream_ptr << "Warning: qt_interp point essentially zero!" << endl
								<< "\t qt_interp = " << qt_interp << " out of {q_min, q_max} = {" << q_min << ", " << q_max << "};" << endl
								<< "\t " << ipt << "   " << ipphi << "   " << iqx << "   " << iqy << "   " << iqz << endl;
		get_CF_terms(&tmpC, &tmpCt, &tmpCct, &tmpCr, ipt, ipphi, iqt0, iqx0, iqy0, iqz0, project_CF && !thermal_pions_only);
		*totalresult = tmpC;
		*thermalresult = tmpCt;
		*CTresult = tmpCct;
		*resonanceresult = tmpCr;
	}
	else if (!q_point_is_outside_grid)
	{
		//assumes qt-grid has already been computed at (adjusted) Chebyshev nodes!!!
		if (QT_POINTS_SPACING == 1 && interp_flag == 0)
		{
			double C_at_q[qtnpts], Ct_at_q[qtnpts], Cct_at_q[qtnpts], Cr_at_q[qtnpts];	//C - 1
			double tmpC = 0.0, tmpCt = 0.0, tmpCct = 0.0, tmpCr = 0.0;
	
			// set CF values along qt-slice for interpolation
			for (int iqtidx = 0; iqtidx < qtnpts; ++iqtidx)
			{
				//return C - 1!!!
				get_CF_terms(&tmpC, &tmpCt, &tmpCct, &tmpCr, ipt, ipphi, iqtidx, iqx, iqy, iqz, project_CF && !thermal_pions_only);
				C_at_q[iqtidx] = tmpC;
				Ct_at_q[iqtidx] = tmpCt;
				Cct_at_q[iqtidx] = tmpCct;
				Cr_at_q[iqtidx] = tmpCr;
			}

			//set up Chebyshev calculation
			int npts_loc[1] = { qtnpts };
			int os[1] = { qtnpts - 1 };
			double lls[1] = { q_min };
			double uls[1] = { q_max };
			double point[1] = { qt_interp };
			int dim_loc = 1;

			Chebyshev cf(C_at_q, npts_loc, os, lls, uls, dim_loc);
			Chebyshev cft(Ct_at_q, npts_loc, os, lls, uls, dim_loc);
			Chebyshev cfct(Cct_at_q, npts_loc, os, lls, uls, dim_loc);
			Chebyshev cfr(Cr_at_q, npts_loc, os, lls, uls, dim_loc);

			*totalresult = cf.eval(point);
			*thermalresult = cft.eval(point);
			*CTresult = cfct.eval(point);
			*resonanceresult = cfr.eval(point);

			if (*thermalresult < 0.0 || *resonanceresult < 0.0 )
			{
				cerr << "WARNING: " << qt_interp << "   " << ipt << "   " << ipphi << "   " << iqx << "   " << iqy << "   " << iqz << "   "
						<< *totalresult << "   " << *thermalresult << "   " << *CTresult << "   " << *resonanceresult << endl;
				for (int iqtidx = 0; iqtidx < qtnpts; ++iqtidx)
				{
					get_CF_terms(&tmpC, &tmpCt, &tmpCct, &tmpCr, ipt, ipphi, iqtidx, iqx, iqy, iqz, project_CF && !thermal_pions_only);
					C_at_q[iqtidx] = tmpC;
					Ct_at_q[iqtidx] = tmpCt;
					Cct_at_q[iqtidx] = tmpCct;
					Cr_at_q[iqtidx] = tmpCr;
					cerr << "Check CF terms: " << qt_pts[iqtidx] << "   " << ipt << "   " << ipphi << "   " << iqx << "   " << iqy << "   " << iqz << "   "
							<< tmpC << "   " << tmpCt << "   " << tmpCct << "   " << tmpCr << endl;
				}
			}
		}
		else if (QT_POINTS_SPACING == 2 && interp_flag == 0)	//use a different distribution of qt points for better accuracy
		{
			const int n = (qtnpts + 1) / 2;
			double tmptan = tan(M_PI / (4.0*n));
			double q_min_local = ( qt_interp <= 0.0 ) ? q_min : q_min*tmptan*tmptan ;
			double q_max_local = ( qt_interp <= 0.0 ) ? q_max*tmptan*tmptan : q_max ;

			double C_at_q[n], Ct_at_q[n], Cct_at_q[n], Cr_at_q[n];	//C - 1
			double tmpC = 0.0, tmpCt = 0.0, tmpCct = 0.0, tmpCr = 0.0;
	
			const int iqt_i = ( qt_interp <= 0.0 ) ? 0 : n-1 ;
			const int iqt_f = ( qt_interp <= 0.0 ) ? n-1 : qtnpts-1 ;
			const int shift = ( qt_interp <= 0.0 ) ? 0 : n-1 ;

			// set CF values along qt-slice for interpolation
			for (int iqtidx = iqt_i; iqtidx <= iqt_f; ++iqtidx)
			{
				//return C - 1!!!
				get_CF_terms(&tmpC, &tmpCt, &tmpCct, &tmpCr, ipt, ipphi, iqtidx, iqx, iqy, iqz, project_CF && !thermal_pions_only);
				C_at_q[iqtidx-shift] = tmpC;
				Ct_at_q[iqtidx-shift] = tmpCt;
				Cct_at_q[iqtidx-shift] = tmpCct;
				Cr_at_q[iqtidx-shift] = tmpCr;
			}

			//set up Chebyshev calculation
			int npts_loc[1] = { n };
			int os[1] = { n - 1 };
			double lls[1] = { q_min_local };
			double uls[1] = { q_max_local };
			double point[1] = { qt_interp };
			int dim_loc = 1;

			Chebyshev cf(C_at_q, npts_loc, os, lls, uls, dim_loc);
			Chebyshev cft(Ct_at_q, npts_loc, os, lls, uls, dim_loc);
			Chebyshev cfct(Cct_at_q, npts_loc, os, lls, uls, dim_loc);
			Chebyshev cfr(Cr_at_q, npts_loc, os, lls, uls, dim_loc);
			*totalresult = cf.eval(point);
			*thermalresult = cft.eval(point);
			*CTresult = cfct.eval(point);
			*resonanceresult = cfr.eval(point);

			
			//if (*thermalresult < 0.0 || *resonanceresult < 0.0 )
			//{
				cerr << "WARNING(OKAY): " << qt_interp << "   " << ipt << "   " << ipphi << "   " << iqx << "   " << iqy << "   " << iqz << "   "
						<< *totalresult << "   " << *thermalresult << "   " << *CTresult << "   " << *resonanceresult << endl;
				for (int iqtidx = iqt_i; iqtidx <= iqt_f; ++iqtidx)
				{
					get_CF_terms(&tmpC, &tmpCt, &tmpCct, &tmpCr, ipt, ipphi, iqtidx, iqx, iqy, iqz, project_CF && !thermal_pions_only);
					C_at_q[iqtidx-shift] = tmpC;
					Ct_at_q[iqtidx-shift] = tmpCt;
					Cct_at_q[iqtidx-shift] = tmpCct;
					Cr_at_q[iqtidx-shift] = tmpCr;
					cerr << "Check CF terms: " << qt_pts[iqtidx] << "   " << ipt << "   " << ipphi << "   " << iqx << "   " << iqy << "   " << iqz << "   "
							<< tmpC << "   " << tmpCt << "   " << tmpCct << "   " << tmpCr << endl;
				}
			//}
			
		}
		else if (QT_POINTS_SPACING == 3 && interp_flag == 0)	//use a different distribution of qt points for better accuracy
		{
			const int n = qtnpts;

			double tmptan = 25.0*tan(M_PI/(2.0*n));
			double x0_loc = 0.0;
			double L_loc = tmptan * q_min;

			double C_at_q[n], Ct_at_q[n], Cct_at_q[n], Cr_at_q[n];	//C - 1
			double tmpC = 0.0, tmpCt = 0.0, tmpCct = 0.0, tmpCr = 0.0;
	
			// set CF values along qt-slice for interpolation
			for (int iqtidx = 0; iqtidx < qtnpts; ++iqtidx)
			{
				//return C - 1!!!
				get_CF_terms(&tmpC, &tmpCt, &tmpCct, &tmpCr, ipt, ipphi, iqtidx, iqx, iqy, iqz, project_CF && !thermal_pions_only);
				C_at_q[iqtidx] = tmpC;
				Ct_at_q[iqtidx] = tmpCt;
				Cct_at_q[iqtidx] = tmpCct;
				Cr_at_q[iqtidx] = tmpCr;
			}

			//set up Chebyshev calculation
			int npts_loc[1] = { n };
			int os[1] = { n - 1 };
			double lls[1] = { x0_loc };
			double uls[1] = { L_loc };
			int modes[1] = { 2 };
			double point[1] = { qt_interp };
			int dim_loc = 1;

			Chebyshev cf(C_at_q, npts_loc, os, lls, uls, dim_loc, modes);
			Chebyshev cft(Ct_at_q, npts_loc, os, lls, uls, dim_loc, modes);
			Chebyshev cfct(Cct_at_q, npts_loc, os, lls, uls, dim_loc, modes);
			Chebyshev cfr(Cr_at_q, npts_loc, os, lls, uls, dim_loc, modes);
			*totalresult = cf.eval(point);
			*thermalresult = cft.eval(point);
			*CTresult = cfct.eval(point);
			*resonanceresult = cfr.eval(point);

			cerr << "WARNING: " << qt_interp << "   " << ipt << "   " << ipphi << "   " << iqx << "   " << iqy << "   " << iqz << "   "
					<< *totalresult << "   " << *thermalresult << "   " << *CTresult << "   " << *resonanceresult << endl;
			for (int iqtidx = 0; iqtidx < qtnpts; ++iqtidx)
			{
				get_CF_terms(&tmpC, &tmpCt, &tmpCct, &tmpCr, ipt, ipphi, iqtidx, iqx, iqy, iqz, project_CF && !thermal_pions_only);
				C_at_q[iqtidx] = tmpC;
				Ct_at_q[iqtidx] = tmpCt;
				Cct_at_q[iqtidx] = tmpCct;
				Cr_at_q[iqtidx] = tmpCr;
				cerr << "Check CF terms: " << qt_pts[iqtidx] << "   " << ipt << "   " << ipphi << "   " << iqx << "   " << iqy << "   " << iqz << "   "
						<< tmpC << "   " << tmpCt << "   " << tmpCct << "   " << tmpCr << endl;
			}
		}
		else if (QT_POINTS_SPACING == 4 && interp_flag == 0)	//use a different distribution of qt points for better accuracy
		{
			const int n = (qtnpts + 1) / 2;

			double tmpsin1 = sin(M_PI / (4.0*n));
			double tmpsin2 = sin(3.0*M_PI / (4.0*n));
			double tmpcos1 = cos(M_PI / (4.0*n));
			double tmpcos2 = cos(M_PI / (2.0*n));
			double tmpcos3 = cos(3.0*M_PI / (2.0*n));

			double prefactor_x0 = -2.0 * tmpsin1*tmpsin1*tmpsin2*tmpsin2 / (tmpcos2+tmpcos3);
			double prefactor_L = 2.0 * tmpcos1*tmpcos1*tmpsin2*tmpsin2 / (tmpcos2+tmpcos3);
			double x0_loc = ( qt_interp <= 0.0 ) ? prefactor_x0 * q_min : prefactor_x0 * q_max ;
			double L_loc = ( qt_interp <= 0.0 ) ? prefactor_L * q_min : prefactor_L * q_max ;

			double C_at_q[n], Ct_at_q[n], Cct_at_q[n], Cr_at_q[n];	//C - 1
			double tmpC = 0.0, tmpCt = 0.0, tmpCct = 0.0, tmpCr = 0.0;
	
			const int iqt_i = ( qt_interp <= 0.0 ) ? 0 : n-1 ;
			const int iqt_f = ( qt_interp <= 0.0 ) ? n-1 : qtnpts-1 ;
			const int shift = ( qt_interp <= 0.0 ) ? 0 : n-1 ;

			// set CF values along qt-slice for interpolation
			for (int iqtidx = iqt_i; iqtidx <= iqt_f; ++iqtidx)
			{
				//return C - 1!!!
				get_CF_terms(&tmpC, &tmpCt, &tmpCct, &tmpCr, ipt, ipphi, iqtidx, iqx, iqy, iqz, project_CF && !thermal_pions_only);
				C_at_q[iqtidx-shift] = tmpC;
				Ct_at_q[iqtidx-shift] = tmpCt;
				Cct_at_q[iqtidx-shift] = tmpCct;
				Cr_at_q[iqtidx-shift] = tmpCr;
			}

			//set up Chebyshev calculation
			int npts_loc[1] = { n };
			int os[1] = { n - 1 };
			double lls[1] = { x0_loc };
			double uls[1] = { L_loc };
			int modes[1] = { 1 };
			double point[1] = { qt_interp };
			int dim_loc = 1;

			Chebyshev cf(C_at_q, npts_loc, os, lls, uls, dim_loc, modes);
			Chebyshev cft(Ct_at_q, npts_loc, os, lls, uls, dim_loc, modes);
			Chebyshev cfct(Cct_at_q, npts_loc, os, lls, uls, dim_loc, modes);
			Chebyshev cfr(Cr_at_q, npts_loc, os, lls, uls, dim_loc, modes);
			*totalresult = cf.eval(point);
			*thermalresult = cft.eval(point);
			*CTresult = cfct.eval(point);
			*resonanceresult = cfr.eval(point);

			cerr << "WARNING(OKAY): " << qt_interp << "   " << ipt << "   " << ipphi << "   " << iqx << "   " << iqy << "   " << iqz << "   "
					<< *totalresult << "   " << *thermalresult << "   " << *CTresult << "   " << *resonanceresult << endl;
			for (int iqtidx = iqt_i; iqtidx <= iqt_f; ++iqtidx)
			{
				get_CF_terms(&tmpC, &tmpCt, &tmpCct, &tmpCr, ipt, ipphi, iqtidx, iqx, iqy, iqz, project_CF && !thermal_pions_only);
				C_at_q[iqtidx-shift] = tmpC;
				Ct_at_q[iqtidx-shift] = tmpCt;
				Cct_at_q[iqtidx-shift] = tmpCct;
				Cr_at_q[iqtidx-shift] = tmpCr;
				cerr << "Check CF terms: " << qt_pts[iqtidx] << "   " << ipt << "   " << ipphi << "   " << iqx << "   " << iqy << "   " << iqz << "   "
						<< tmpC << "   " << tmpCt << "   " << tmpCct << "   " << tmpCr << endl;
			}
		}
		else	//if not using Chebyshev nodes in qt-direction, just use straight-up linear(0) or cubic(1) interpolation
		{
			double C_at_q[qtnpts], Ct_at_q[qtnpts], Cct_at_q[qtnpts], Cr_at_q[qtnpts];	//C - 1
			double tmpC = 0.0, tmpCt = 0.0, tmpCct = 0.0, tmpCr = 0.0;
debugger(__LINE__, __FILE__);

			// set CF values along qt-slice for interpolation
			for (int iqtidx = 0; iqtidx < qtnpts; ++iqtidx)
			{
				//return C - 1!!!
				get_CF_terms(&tmpC, &tmpCt, &tmpCct, &tmpCr, ipt, ipphi, iqtidx, iqx, iqy, iqz, project_CF && !thermal_pions_only);
				C_at_q[iqtidx] = tmpC;
				Ct_at_q[iqtidx] = tmpCt;
				Cct_at_q[iqtidx] = tmpCct;
				Cr_at_q[iqtidx] = tmpCr;
			}
debugger(__LINE__, __FILE__);
			*totalresult = interpolate1D(qt_pts, C_at_q, qt_interp, qtnpts, 1, false);
			*thermalresult = interpolate1D(qt_pts, Ct_at_q, qt_interp, qtnpts, 1, false);
			*CTresult = interpolate1D(qt_pts, Cct_at_q, qt_interp, qtnpts, 1, false);
			*resonanceresult = interpolate1D(qt_pts, Cr_at_q, qt_interp, qtnpts, 1, false);
		}
	}
	else
	{
		*global_out_stream_ptr << "Warning: qt_interp point was outside of computed grid!" << endl
								<< "\t qt_interp = " << qt_interp << " out of {q_min, q_max} = {" << q_min << ", " << q_max << "};" << endl
								<< "\t " << ipt << "   " << ipphi << "   " << iqx << "   " << iqy << "   " << iqz << endl;
		*totalresult = 0.0;
	}
	return;
}


//******************************************************************
// Routines for refining correlation function grid via interpolation
//******************************************************************

void CorrelationFunction::Flesh_out_CF(int ipt, int ipphi, double sample_scale /*==1.0*/)
{
	//declare needed quantities here
	double qxmin = 0.9999*sample_scale*qx_pts[0], qxmax = 0.9999*sample_scale*qx_pts[qxnpts-1];
	double qymin = 0.9999*sample_scale*qy_pts[0], qymax = 0.9999*sample_scale*qy_pts[qynpts-1];
	double qzmin = 0.9999*sample_scale*qz_pts[0], qzmax = 0.9999*sample_scale*qz_pts[qznpts-1];

	new_nqxpts = ( qxnpts > 1 ) ? new_nqpts : 1;
	new_nqypts = ( qynpts > 1 ) ? new_nqpts : 1;
	new_nqzpts = ( qznpts > 1 ) ? new_nqpts : 1;

	double new_Del_qx = (qxmax - qxmin)/(double(new_nqxpts-1)+1.e-100);
	double new_Del_qy = (qymax - qymin)/(double(new_nqypts-1)+1.e-100);
	double new_Del_qz = (qzmax - qzmin)/(double(new_nqzpts-1)+1.e-100);

	double *** current_C_slice = CFvals[ipt][ipphi];

	//cout << "(qxmin, qxmax, new_Del_qx) = (" << qxmin << ", " << qxmax << ", " << new_Del_qx << ")" << endl;
	//cout << "(qymin, qymax, new_Del_qy) = (" << qymin << ", " << qymax << ", " << new_Del_qy << ")" << endl;
	//cout << "(qzmin, qzmax, new_Del_qz) = (" << qzmin << ", " << qzmax << ", " << new_Del_qz << ")" << endl;
	
	for (int iqx = 0; iqx < new_nqpts; ++iqx)
	for (int iqy = 0; iqy < new_nqpts; ++iqy)
	for (int iqz = 0; iqz < new_nqpts; ++iqz)
	{
		double qx0 = qxmin + double(iqx) * new_Del_qx;
		double qy0 = qymin + double(iqy) * new_Del_qy;
		double qz0 = qzmin + double(iqz) * new_Del_qz;

		qx_fleshed_out_pts[iqx] = qx0;
		qy_fleshed_out_pts[iqy] = qy0;
		qz_fleshed_out_pts[iqz] = qz0;

		fleshed_out_thermal[iqx][iqy][iqz] = interpolate_CF(thermalCFvals[ipt][ipphi], qx0, qy0, qz0, ipt, 0);
		fleshed_out_crossterm[iqx][iqy][iqz] = interpolate_CF(crosstermCFvals[ipt][ipphi], qx0, qy0, qz0, ipt, 1);
		fleshed_out_resonances[iqx][iqy][iqz] = interpolate_CF(resonancesCFvals[ipt][ipphi], qx0, qy0, qz0, ipt, 2);
		fleshed_out_CF[iqx][iqy][iqz] = 1.0 + fleshed_out_thermal[iqx][iqy][iqz] + fleshed_out_crossterm[iqx][iqy][iqz] + fleshed_out_resonances[iqx][iqy][iqz];
	}

	//set up Chebyshev calculation
	/*
	const int n = (qxnpts + 1) / 2;
	double tmptan = tan(M_PI / (4.0*n));
	double qx_min_local_neg = -qx_max, qx_min_local_pos = -qx_max*tmptan*tmptan ;
	double qx_max_local_neg = qx_max*tmptan*tmptan, qx_max_local_pos = qx_max ;

	int npts_loc[1] = { n };
	int os[1] = { n - 1 };
	double lls[1] = { qx_min_local_neg };
	double uls[1] = { qx_max_local_neg };
	int dim_loc = 1;

	double C_at_neg_q[n];	//C - 1
	for (int iqd = 0; iqd < n; iqd++)
		C_at_neg_q[iqd] = current_C_slice[iqd][0][0];
	Chebyshev cfNEG(C_at_neg_q, npts_loc, os, lls, uls, dim_loc);

	double C_at_pos_q[n];	//C - 1
	lls[0] = qx_min_local_pos;
	uls[0] = qx_max_local_pos;
	for (int iqd = n-1; iqd < qxnpts; iqd++)
		C_at_pos_q[iqd-(n-1)] = current_C_slice[iqd][0][0];
	Chebyshev cfPOS(C_at_pos_q, npts_loc, os, lls, uls, dim_loc);

	for (int iqx = 0; iqx < (new_nqxpts+1)/2; ++iqx)
	{
		double qx0 = qxmin + double(iqx) * new_Del_qx;

		qx_fleshed_out_pts[iqx] = qx0;
		double point[1] = { qx0 };

		fleshed_out_CF[iqx][0][0] = cfNEG.eval(point);
	}

	for (int iqx = (new_nqxpts+1)/2; iqx < new_nqxpts; ++iqx)
	{
		double qx0 = qxmin + double(iqx) * new_Del_qx;

		qx_fleshed_out_pts[iqx] = qx0;
		double point[1] = { qx0 };

		const int iqt_i = ( qx0 <= 0.0 ) ? 0 : n-1 ;
		const int iqt_f = ( qx0 <= 0.0 ) ? n-1 : qxnpts-1 ;
		const int shift = ( qx0 <= 0.0 ) ? 0 : n-1 ;

		fleshed_out_CF[iqx][0][0] = cfPOS.eval(point);
	}
	*/

	return;
}

double CorrelationFunction::interpolate_CF(double *** current_C_slice, double qx0, double qy0, double qz0, int ipt, int thermal_or_resonances)
{
	//int thermal_or_resonances - 	0: interpolate thermal CF (assuming basically Gaussian)
	//								1: interpolate resonance contributions (linear near origin, exponential further out)

	int iqx0_loc = binarySearch(qx_pts, qxnpts, qx0, true, true);
	int iqy0_loc = binarySearch(qy_pts, qynpts, qy0, true, true);
	int iqz0_loc = binarySearch(qz_pts, qznpts, qz0, true, true);

	if (iqx0_loc == -1 || iqy0_loc == -1 || iqz0_loc == -1)
	{
		cerr << "Interpolation failed: exiting!" << endl;
		exit(1);
	}

	double fx0y0z0 = current_C_slice[iqx0_loc][iqy0_loc][iqz0_loc];
	double fx0y0z1 = current_C_slice[iqx0_loc][iqy0_loc][iqz0_loc+1];
	double fx0y1z0 = current_C_slice[iqx0_loc][iqy0_loc+1][iqz0_loc];
	double fx0y1z1 = current_C_slice[iqx0_loc][iqy0_loc+1][iqz0_loc+1];
	double fx1y0z0 = current_C_slice[iqx0_loc+1][iqy0_loc][iqz0_loc];
	double fx1y0z1 = current_C_slice[iqx0_loc+1][iqy0_loc][iqz0_loc+1];
	double fx1y1z0 = current_C_slice[iqx0_loc+1][iqy0_loc+1][iqz0_loc];
	double fx1y1z1 = current_C_slice[iqx0_loc+1][iqy0_loc+1][iqz0_loc+1];

	double fx0 = current_C_slice[iqx0_loc][0][0];
	double fx1 = current_C_slice[iqx0_loc+1][0][0];

	//interpolate here...
	double qx_0 = qx_pts[iqx0_loc];
	double qx_1 = qx_pts[iqx0_loc+1];
	double qy_0 = qy_pts[iqy0_loc];
	double qy_1 = qy_pts[iqy0_loc+1];
	double qz_0 = qz_pts[iqz0_loc];
	double qz_1 = qz_pts[iqz0_loc+1];

	double fxiyizi = 0.0;
	if (thermal_or_resonances == 0)
	{
		double sqx02 = sgn(qx0)*qx0*qx0;
		double sqy02 = sgn(qy0)*qy0*qy0;
		double sqz02 = sgn(qz0)*qz0*qz0;
		double sqx_02 = sgn(qx_0)*qx_0*qx_0;
		double sqy_02 = sgn(qy_0)*qy_0*qy_0;
		double sqz_02 = sgn(qz_0)*qz_0*qz_0;
		double sqx_12 = sgn(qx_1)*qx_1*qx_1;
		double sqy_12 = sgn(qy_1)*qy_1*qy_1;
		double sqz_12 = sgn(qz_1)*qz_1*qz_1;

		//interpolate over qz-points first
		///////////////////////
		//interpolate over each pair of qz-points
		double fx0y0zi = interpolate_qi(sqz02, sqz_02, sqz_12, fx0y0z0, fx0y0z1, false);
		double fx0y1zi = interpolate_qi(sqz02, sqz_02, sqz_12, fx0y1z0, fx0y1z1, false);    
		double fx1y0zi = interpolate_qi(sqz02, sqz_02, sqz_12, fx1y0z0, fx1y0z1, false);
		double fx1y1zi = interpolate_qi(sqz02, sqz_02, sqz_12, fx1y1z0, fx1y1z1, false);
		///////////////////////
	
		//interpolate over qy-points next
		double fx0yizi = interpolate_qi(sqy02, sqy_02, sqy_12, fx0y0zi, fx0y1zi, false);
		double fx1yizi = interpolate_qi(sqy02, sqy_02, sqy_12, fx1y0zi, fx1y1zi, false);
	
		//finally, interpolate over qx-points
		fxiyizi = interpolate_qi(sqx02, sqx_02, sqx_12, fx0yizi, fx1yizi, false);
	}
	else
	{
		fxiyizi = interpolate3D(qx_pts, qy_pts, qz_pts, current_C_slice, qx0, qy0, qz0, qxnpts, qynpts, qznpts, 1, true);
	}

	return (fxiyizi);
}

double CorrelationFunction::interpolate_qi(double q0, double qi0, double qi1, double f1, double f2, bool use_linear)
{
	double if1 = f1;
	double if2 = f2;
	bool use_log = (!use_linear) && (f1 > 0.0) && (f2 > 0.0);
	if (use_log)
	{
		if1 = log(f1);
		if2 = log(f2);
	}
	double tmp_result = lin_int(q0 - qi0, 1./(qi1-qi0), if1, if2);
	if (use_log)
		tmp_result = exp(tmp_result);

	return (tmp_result);
}


double CorrelationFunction::get_CF(int ipt, int ipphi, int iqt, int iqx, int iqy, int iqz, bool return_projected_value)
{	//pY==0
	double nonFTd_spectra = spectra[target_particle_id][ipt][ipphi];
	double cos_transf_spectra = full_target_Yeq0_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,0)]
									+ full_target_Yeq0_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,3)];		//add real components
	double sin_transf_spectra = full_target_Yeq0_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,1)]
									+ full_target_Yeq0_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,2)];		//add imaginary components

	if (return_projected_value)
	{
		//with no resonances
		double nonFTd_tspectra = thermal_spectra[target_particle_id][ipt][ipphi];
		double cos_transf_tspectra = thermal_target_Yeq0_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,0)]
										+ thermal_target_Yeq0_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,3)];	//add real components
		double sin_transf_tspectra = thermal_target_Yeq0_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,1)]
										+ thermal_target_Yeq0_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,2)];	//add imaginary components

		double projected_nonFTd_spectra = nonFTd_tspectra + (nonFTd_spectra - nonFTd_tspectra) / fraction_of_resonances;
		double projected_cos_transf_spectra = cos_transf_tspectra + (cos_transf_spectra - cos_transf_tspectra) / fraction_of_resonances;
		double projected_sin_transf_spectra = sin_transf_tspectra + (sin_transf_spectra - sin_transf_tspectra) / fraction_of_resonances;

		//projected full value of correlation function by projecting moments *first*
		double projected_num = projected_cos_transf_spectra*projected_cos_transf_spectra + projected_sin_transf_spectra*projected_sin_transf_spectra;
		double projected_den = projected_nonFTd_spectra*projected_nonFTd_spectra;
		return (1. + projected_num / projected_den);
	}
	else
	{
		double num = cos_transf_spectra*cos_transf_spectra + sin_transf_spectra*sin_transf_spectra;
		double den = nonFTd_spectra*nonFTd_spectra;
		return (1. + num / den);
	}
}

void CorrelationFunction::get_CF_terms(double * totalresult, double * thermalresult, double * crosstermresult, double * resonanceresult,
									int ipt, int ipphi, int iqt, int iqx, int iqy, int iqz, bool return_projected_value)
{
	//thermal
	double nonFTd_tspectra = thermal_spectra[target_particle_id][ipt][ipphi];
	double cos_transf_tspectra = thermal_target_Yeq0_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,0)]
									+ thermal_target_Yeq0_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,3)];	//add real components
	double sin_transf_tspectra = thermal_target_Yeq0_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,1)]
									+ thermal_target_Yeq0_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,2)];	//add imaginary components
	//total
	double nonFTd_spectra = spectra[target_particle_id][ipt][ipphi];
	double cos_transf_spectra = full_target_Yeq0_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,0)]
									+ full_target_Yeq0_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,3)];		//add real components
	double sin_transf_spectra = full_target_Yeq0_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,1)]
									+ full_target_Yeq0_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,2)];		//add imaginary components

	if (return_projected_value)
	{
		nonFTd_spectra = nonFTd_tspectra + (nonFTd_spectra - nonFTd_tspectra) / fraction_of_resonances;
		cos_transf_spectra = cos_transf_tspectra + (cos_transf_spectra - cos_transf_tspectra) / fraction_of_resonances;
		sin_transf_spectra = sin_transf_tspectra + (sin_transf_spectra - sin_transf_tspectra) / fraction_of_resonances;
	}

	//non-thermal
	double NT_spectra = nonFTd_spectra - nonFTd_tspectra;
	double cosNT_spectra = cos_transf_spectra - cos_transf_tspectra;
	double sinNT_spectra = sin_transf_spectra - sin_transf_tspectra;

	double num = cos_transf_tspectra*cos_transf_tspectra + sin_transf_tspectra*sin_transf_tspectra;
	double den = nonFTd_spectra*nonFTd_spectra;
	*thermalresult = num / den;

	num = 2.0 * (cosNT_spectra*cos_transf_tspectra+sinNT_spectra*sin_transf_tspectra);			//cross term
	*crosstermresult = num / den;

	num = cosNT_spectra*cosNT_spectra+sinNT_spectra*sinNT_spectra;								//resonances only
	*resonanceresult = num / den;

	num = cos_transf_spectra*cos_transf_spectra + sin_transf_spectra*sin_transf_spectra;
	*totalresult = num / den;

	return;
}

//**************************************************************
// Gaussian fit routine below
//**************************************************************
void CorrelationFunction::find_minimum_chisq_correlationfunction_full(
		double *** Correl_3D, int ipt, int ipphi, bool fleshing_out_CF /*== true*/)
{
	double * q1pts = qx_pts;
	double * q2pts = qy_pts;
	double * q3pts = qz_pts;
	if (fleshing_out_CF)
	{
		q1npts = new_nqxpts;
		q2npts = new_nqypts;
		q3npts = new_nqzpts;
		q1pts = qx_fleshed_out_pts;
		q2pts = qy_fleshed_out_pts;
		q3pts = qz_fleshed_out_pts;
	}
	const size_t data_length = q1npts*q2npts*q3npts;  // # of points

    double lambda, R_o, R_s, R_l, R_os;
    int dim = 5;
    int s_gsl;

    double *V = new double [dim];
    double *qweight = new double [dim];
    double **T = new double* [dim];
    for(int i = 0; i < dim; i++)
    {
        V[i] = 0.0;
        T[i] = new double [dim];
        for(int j = 0; j < dim; j++)
            T[i][j] = 0.0;
    }

    gsl_matrix * T_gsl = gsl_matrix_alloc (dim, dim);
    gsl_matrix * T_inverse_gsl = gsl_matrix_alloc (dim, dim);
    gsl_permutation * perm = gsl_permutation_alloc (dim);

	double ckp = cos_SP_pphi[ipphi], skp = sin_SP_pphi[ipphi];
	double CF_err = 1.e-3;
	for (int i = 0; i < q1npts; i++)
	for (int j = 0; j < q2npts; j++)
	for (int k = 0; k < q3npts; k++)
    {
        double q_out_local = q1pts[i] * ckp + q2pts[j] * skp;
        double q_side_local = -q1pts[i] * skp + q2pts[j] * ckp;
        double q_long_local = q3pts[k];
        double correl_local = Correl_3D[i][j][k]-1;
		//*global_out_stream_ptr << "\t\t" << q1pts[i] << "   " << q2pts[i] << "   " << q3pts[i] << "   " << q_out_local << "   " << q_side_local << "   " << q_long_local << "   " << correl_local << endl;
        if(correl_local < 1e-15) continue;
		//if (i==(q1npts-1)/2 && j==(q2npts-1)/2 && k==(q3npts-1)/2)
		//	Correlfun3D_data.sigma[idx] = 1.e10;	//ignore central point
        double sigma_k_prime = CF_err/correl_local;
            
        double inv_sigma_k_prime_sq = 1./(sigma_k_prime*sigma_k_prime);
        double log_correl_over_sigma_sq = log(correl_local)*inv_sigma_k_prime_sq;

        qweight[0] = - 1.0;
        qweight[1] = q_out_local*q_out_local;
        qweight[2] = q_side_local*q_side_local;
        qweight[3] = q_long_local*q_long_local;
        qweight[4] = q_out_local*q_side_local;

        for(int ij = 0; ij < dim; ij++)
        {
            V[ij] += qweight[ij]*log_correl_over_sigma_sq;
            T[0][ij] += qweight[ij]*inv_sigma_k_prime_sq;
        }

        for(int ij = 1; ij < dim; ij++)
            T[ij][0] = T[0][ij];
            

        for(int ij = 1; ij < dim; ij++)
        {
            for(int lm = 1; lm < dim; lm++)
                T[ij][lm] += -qweight[ij]*qweight[lm]*inv_sigma_k_prime_sq;
        }
    }
    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
            gsl_matrix_set(T_gsl, i, j, T[i][j]);

    // Make LU decomposition of matrix T_gsl
    gsl_linalg_LU_decomp (T_gsl, perm, &s_gsl);
    // Invert the matrix m
    gsl_linalg_LU_invert (T_gsl, perm, T_inverse_gsl);

    double **T_inverse = new double* [dim];
    for(int i = 0; i < dim; i++)
    {
        T_inverse[i] = new double [dim];
        for(int j = 0; j < dim; j++)
            T_inverse[i][j] = gsl_matrix_get(T_inverse_gsl, i, j);
    }
    double *results = new double [dim];
    for(int i = 0; i < dim; i++)
    {
        results[i] = 0.0;
        for(int j = 0; j < dim; j++)
            results[i] += T_inverse[i][j]*V[j];
    }

	lambda_Correl[ipt][ipphi] = exp(results[0]);
	lambda_Correl_err[ipt][ipphi] = 0.0;
	R2_out_GF[ipt][ipphi] = results[1]*hbarC*hbarC;
	R2_side_GF[ipt][ipphi] = results[2]*hbarC*hbarC;
	R2_long_GF[ipt][ipphi] = results[3]*hbarC*hbarC;
	R2_outside_GF[ipt][ipphi] = results[4]*hbarC*hbarC;
	R2_out_err[ipt][ipphi] = 0.0;
	R2_side_err[ipt][ipphi] = 0.0;
	R2_long_err[ipt][ipphi] = 0.0;
	R2_outside_err[ipt][ipphi] = 0.0;


    double chi_sq = 0.0;
	for (int i = 0; i < q1npts; i++)
	for (int j = 0; j < q2npts; j++)
	for (int k = 0; k < q3npts; k++)
    {
        double q_out_local = q1pts[i] * ckp + q2pts[j] * skp;
        double q_side_local = -q1pts[i] * skp + q2pts[j] * ckp;
        double q_long_local = q3pts[k];
        double correl_local = Correl_3D[i][j][k]-1;
        if(correl_local < 1e-15) continue;
        double sigma_k_prime = CF_err/correl_local;

        chi_sq += pow((log(correl_local) - results[0] 
                       + results[1]*q_out_local*q_out_local 
                       + results[2]*q_side_local*q_side_local
                       + results[3]*q_long_local*q_long_local
                       + results[4]*q_out_local*q_side_local), 2)
                  /sigma_k_prime/sigma_k_prime;
    }
    //cout << "chi_sq/d.o.f = " << chi_sq/(qnpts - dim) << endl;
    //chi_sq_per_dof = chi_sq/(qnpts - dim);

    // clean up
    gsl_matrix_free (T_gsl);
    gsl_matrix_free (T_inverse_gsl);
    gsl_permutation_free (perm);

    delete [] qweight;
    delete [] V;
    for(int i = 0; i < dim; i++)
    {
        delete [] T[i];
        delete [] T_inverse[i];
    }
    delete [] T;
    delete [] T_inverse;
    delete [] results;
}

void CorrelationFunction::Fit_Correlationfunction3D_withlambda(double *** Correl_3D, int ipt, int ipphi, bool fleshing_out_CF /*== true*/)
{
	double * q1pts = qx_pts;
	double * q2pts = qy_pts;
	double * q3pts = qz_pts;
	if (fleshing_out_CF)
	{
		q1npts = new_nqpts;
		q2npts = new_nqpts;
		q3npts = new_nqpts;
		q1pts = qx_fleshed_out_pts;
		q2pts = qy_fleshed_out_pts;
		q3pts = qz_fleshed_out_pts;
	}
	const size_t data_length = q1npts*q2npts*q3npts;  // # of points
	const size_t n_para = 5;  // # of parameters

	// allocate space for a covariance matrix of size p by p
	gsl_matrix *covariance_ptr = gsl_matrix_alloc (n_para, n_para);

	// allocate and setup for generating gaussian distibuted random numbers
	gsl_rng_env_setup ();
	const gsl_rng_type *type = gsl_rng_default;
	gsl_rng *rng_ptr = gsl_rng_alloc (type);

	//set up test data
	struct Correlationfunction3D_data Correlfun3D_data;
	Correlfun3D_data.data_length = data_length;
	Correlfun3D_data.q_o = new double [data_length];
	Correlfun3D_data.q_s = new double [data_length];
	Correlfun3D_data.q_l = new double [data_length];
	Correlfun3D_data.y = new double [data_length];
	Correlfun3D_data.sigma = new double [data_length];

	int idx = 0;
	double ckp = cos_SP_pphi[ipphi], skp = sin_SP_pphi[ipphi];
	for (int i = 0; i < q1npts; i++)
	for (int j = 0; j < q2npts; j++)
	for (int k = 0; k < q3npts; k++)
	{
		Correlfun3D_data.q_o[idx] = q1pts[i] * ckp + q2pts[j] * skp;
		Correlfun3D_data.q_s[idx] = -q1pts[i] * skp + q2pts[j] * ckp;
		Correlfun3D_data.q_l[idx] = q3pts[k];
		Correlfun3D_data.y[idx] = Correl_3D[i][j][k];
		//Correlfun3D_data.sigma[idx] = Correl_3D_err[i][j][k];
		Correlfun3D_data.sigma[idx] = 1.e-3;
//if (i==(q1npts-1)/2 && j==(q2npts-1)/2 && k==(q3npts-1)/2)
//	Correlfun3D_data.sigma[idx] = 1.e10;	//ignore central point
		idx++;
	}
	double para_init[n_para] = { 1.0, 1.0, 1.0, 1.0, 1.0 };  // initial guesses of parameters

	gsl_vector_view xvec_ptr = gsl_vector_view_array (para_init, n_para);
  
	// set up the function to be fit 
	gsl_multifit_function_fdf target_func;
	target_func.f = &Fittarget_correlfun3D_f_withlambda;        // the function of residuals
	target_func.df = &Fittarget_correlfun3D_df_withlambda;      // the gradient of this function
	target_func.fdf = &Fittarget_correlfun3D_fdf_withlambda;    // combined function and gradient
	target_func.n = data_length;              // number of points in the data set
	target_func.p = n_para;              // number of parameters in the fit function
	target_func.params = &Correlfun3D_data;  // structure with the data and error bars

	const gsl_multifit_fdfsolver_type *type_ptr = gsl_multifit_fdfsolver_lmsder;
	gsl_multifit_fdfsolver *solver_ptr = gsl_multifit_fdfsolver_alloc (type_ptr, data_length, n_para);
	gsl_multifit_fdfsolver_set (solver_ptr, &target_func, &xvec_ptr.vector);

	size_t iteration = 0;         // initialize iteration counter
	print_fit_state_3D_withlambda (iteration, solver_ptr);
	int status;  		// return value from gsl function calls (e.g., error)
	do
	{
		iteration++;
      
		// perform a single iteration of the fitting routine
		status = gsl_multifit_fdfsolver_iterate (solver_ptr);

		// print out the status of the fit
		//if (VERBOSE > 2)
			cout << "status = " << gsl_strerror (status) << endl;

		// customized routine to print out current parameters
		//if (VERBOSE > 2)
			print_fit_state_3D_withlambda (iteration, solver_ptr);

		if (status)    // check for a nonzero status code
		{
			break;  // this should only happen if an error code is returned 
		}

		// test for convergence with an absolute and relative error (see manual)
		status = gsl_multifit_test_delta (solver_ptr->dx, solver_ptr->x, fit_tolerance, fit_tolerance);
	}
	while (status == GSL_CONTINUE && iteration < fit_max_iterations);

	// calculate the covariance matrix of the best-fit parameters
	gsl_multifit_covar (solver_ptr->J, 0.0, covariance_ptr);

	// print out the covariance matrix using the gsl function (not elegant!)
	if (VERBOSE > 2)
		cout << endl << "Covariance matrix: " << endl;
	if (VERBOSE > 2)
		gsl_matrix_fprintf (stdout, covariance_ptr, "%g");

	cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
	cout.precision (5);		                // # of digits in doubles

	int width = 7;		// setw width for output
	/*cout << endl << "Best fit results:" << endl;
	cout << "lambda      = " << setw (width) << get_fit_results (0, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (0, covariance_ptr) << endl;

	cout << "R2o = " << setw (width) << get_fit_results (1, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (1, covariance_ptr) << endl;

	cout << "R2s      = " << setw (width) << get_fit_results (2, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (2, covariance_ptr) << endl;

	cout << "R2l      = " << setw (width) << get_fit_results (3, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (3, covariance_ptr) << endl;
  
	cout << "R2os      = " << setw (width) << get_fit_results (4, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (4, covariance_ptr) << endl;
    
	cout << "status = " << gsl_strerror (status) << endl;
	cout << "--------------------------------------------------------------------" << endl;*/

	double chi = gsl_blas_dnrm2(solver_ptr->f);
	double dof = data_length - n_para;
	double c = GSL_MAX_DBL(1, chi/sqrt(dof));

	lambda_Correl[ipt][ipphi] = get_fit_results(0, solver_ptr);
	lambda_Correl_err[ipt][ipphi] = c*get_fit_err(0, covariance_ptr);
	R2_out_GF[ipt][ipphi] = fabs(get_fit_results(1, solver_ptr))*hbarC*hbarC;
	R2_side_GF[ipt][ipphi] = fabs(get_fit_results(2, solver_ptr))*hbarC*hbarC;
	R2_long_GF[ipt][ipphi] = fabs(get_fit_results(3, solver_ptr))*hbarC*hbarC;
	R2_outside_GF[ipt][ipphi] = get_fit_results(4, solver_ptr)*hbarC*hbarC;
	R2_out_err[ipt][ipphi] = c*get_fit_err(1, covariance_ptr)*hbarC*hbarC;
	R2_side_err[ipt][ipphi] = c*get_fit_err(2, covariance_ptr)*hbarC*hbarC;
	R2_long_err[ipt][ipphi] = c*get_fit_err(3, covariance_ptr)*hbarC*hbarC;
	R2_outside_err[ipt][ipphi] = c*get_fit_err(4, covariance_ptr)*hbarC*hbarC;

	//if (VERBOSE > 1) 
	//{
		cout << "final results: " << endl;
		cout << scientific << setw(10) << setprecision(5) 
			<< "chisq/dof = " << chi*chi/dof << endl;
		cout << scientific << setw(10) << setprecision(5) 
			<< " lambda = " << lambda_Correl[ipt][ipphi] << " +/- " << lambda_Correl_err[ipt][ipphi] << endl;
		cout << " R2_out[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << R2_out_GF[ipt][ipphi] << " +/- " << R2_out_err[ipt][ipphi] << endl;
		cout << " R2_side[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << R2_side_GF[ipt][ipphi] << " +/- " << R2_side_err[ipt][ipphi] << endl;
		cout << " R2_long[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << R2_long_GF[ipt][ipphi] << " +/- " << R2_long_err[ipt][ipphi] << endl;
		cout << " R2_outside[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << R2_outside_GF[ipt][ipphi] << " +/- " << R2_outside_err[ipt][ipphi] << endl;
	//}

	//clean up
	gsl_matrix_free (covariance_ptr);
	gsl_rng_free (rng_ptr);

	delete[] Correlfun3D_data.q_o;
	delete[] Correlfun3D_data.q_s;
	delete[] Correlfun3D_data.q_l;
	delete[] Correlfun3D_data.y;
	delete[] Correlfun3D_data.sigma;

	gsl_multifit_fdfsolver_free (solver_ptr);  // free up the solver

	return;
}


//Fourier transform of HBT radii once they're calculated
void CorrelationFunction::R2_Fourier_transform(int iKT, double plane_psi, int mode)
{
	const int interpMode = 1;
	//int mode: 0 - GF, 1 - QM
	for(int Morder = 0; Morder < n_order; ++Morder)
	{
		double cos_mK_phi[nKphi], sin_mK_phi[nKphi];

		for(int iKphi = 0; iKphi < nKphi; ++iKphi)
		{
			cos_mK_phi[iKphi] = cos(Morder*(K_phi[iKphi] - plane_psi));
			sin_mK_phi[iKphi] = sin(Morder*(K_phi[iKphi] - plane_psi));
		}

		double temp_sum_side_cos = 0.0, temp_sum_side_sin = 0.0;
		double temp_sum_out_cos = 0.0, temp_sum_out_sin = 0.0;
		double temp_sum_outside_cos = 0.0, temp_sum_outside_sin = 0.0;
		double temp_sum_long_cos = 0.0, temp_sum_long_sin = 0.0;
		double temp_sum_sidelong_cos = 0.0, temp_sum_sidelong_sin = 0.0;
		double temp_sum_outlong_cos = 0.0, temp_sum_outlong_sin = 0.0;

		for(int iKphi = 0; iKphi < nKphi; ++iKphi)
		{
			double local_R2s = interpolate2D(SP_pT, SP_pphi, R2_side_GF, K_T[iKT], K_phi[iKphi], n_pT_pts, n_pphi_pts, interpMode, false, true);
			double local_R2o = interpolate2D(SP_pT, SP_pphi, R2_out_GF, K_T[iKT], K_phi[iKphi], n_pT_pts, n_pphi_pts, interpMode, false, true);
			double local_R2os = interpolate2D(SP_pT, SP_pphi, R2_outside_GF, K_T[iKT], K_phi[iKphi], n_pT_pts, n_pphi_pts, interpMode, false, true);
			double local_R2l = interpolate2D(SP_pT, SP_pphi, R2_long_GF, K_T[iKT], K_phi[iKphi], n_pT_pts, n_pphi_pts, interpMode, false, true);
			double local_R2sl = interpolate2D(SP_pT, SP_pphi, R2_sidelong_GF, K_T[iKT], K_phi[iKphi], n_pT_pts, n_pphi_pts, interpMode, false, true);
			double local_Rol = interpolate2D(SP_pT, SP_pphi, R2_outlong_GF, K_T[iKT], K_phi[iKphi], n_pT_pts, n_pphi_pts, interpMode, false, true);
			temp_sum_side_cos += local_R2s*cos_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_side_sin += local_R2s*sin_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_out_cos += local_R2o*cos_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_out_sin += local_R2o*sin_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_outside_cos += local_R2os*cos_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_outside_sin += local_R2os*sin_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_long_cos += local_R2l*cos_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_long_sin += local_R2l*sin_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_sidelong_cos += local_R2sl*cos_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_sidelong_sin += local_R2sl*sin_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_outlong_cos += local_Rol*cos_mK_phi[iKphi]*K_phi_weight[iKphi];
			temp_sum_outlong_sin += local_Rol*sin_mK_phi[iKphi]*K_phi_weight[iKphi];
		}

		if (mode == 0)
		{
			R2_side_GF_C[iKT][Morder] = temp_sum_side_cos/(2.*M_PI);
			R2_side_GF_S[iKT][Morder] = temp_sum_side_sin/(2.*M_PI);
			R2_out_GF_C[iKT][Morder] = temp_sum_out_cos/(2.*M_PI);
			R2_out_GF_S[iKT][Morder] = temp_sum_out_sin/(2.*M_PI);
			R2_outside_GF_C[iKT][Morder] = temp_sum_outside_cos/(2.*M_PI);
			R2_outside_GF_S[iKT][Morder] = temp_sum_outside_sin/(2.*M_PI);
			R2_long_GF_C[iKT][Morder] = temp_sum_long_cos/(2.*M_PI);
			R2_long_GF_S[iKT][Morder] = temp_sum_long_sin/(2.*M_PI);
			R2_sidelong_GF_C[iKT][Morder] = temp_sum_sidelong_cos/(2.*M_PI);
			R2_sidelong_GF_S[iKT][Morder] = temp_sum_sidelong_sin/(2.*M_PI);
			R2_outlong_GF_C[iKT][Morder] = temp_sum_outlong_cos/(2.*M_PI);
			R2_outlong_GF_S[iKT][Morder] = temp_sum_outlong_sin/(2.*M_PI);
		}
		else if (mode == 1)
		{
			R2_side_QM_C[iKT][Morder] = temp_sum_side_cos/(2.*M_PI);
			R2_side_QM_S[iKT][Morder] = temp_sum_side_sin/(2.*M_PI);
			R2_out_QM_C[iKT][Morder] = temp_sum_out_cos/(2.*M_PI);
			R2_out_QM_S[iKT][Morder] = temp_sum_out_sin/(2.*M_PI);
			R2_outside_QM_C[iKT][Morder] = temp_sum_outside_cos/(2.*M_PI);
			R2_outside_QM_S[iKT][Morder] = temp_sum_outside_sin/(2.*M_PI);
			R2_long_QM_C[iKT][Morder] = temp_sum_long_cos/(2.*M_PI);
			R2_long_QM_S[iKT][Morder] = temp_sum_long_sin/(2.*M_PI);
			R2_sidelong_QM_C[iKT][Morder] = temp_sum_sidelong_cos/(2.*M_PI);
			R2_sidelong_QM_S[iKT][Morder] = temp_sum_sidelong_sin/(2.*M_PI);
			R2_outlong_QM_C[iKT][Morder] = temp_sum_outlong_cos/(2.*M_PI);
			R2_outlong_QM_S[iKT][Morder] = temp_sum_outlong_sin/(2.*M_PI);
		}
	}

	return;
}

void CorrelationFunction::Set_target_moments(int iqt, int iqz)
{
	bool include_thermal_pions = true;
	if (include_thermal_pions)
	{
		/*if ( !thermal_pions_only or
				find(osr.begin(), osr.end(), target_particle_id) != osr.end() )
			return;*/

		//Stopwatch sw;
		//sw.Start();
		cout << "Setting thermal target moments..." << endl;
		Set_thermal_target_moments(iqt, iqz);
		//sw.Stop();
		//cout << "done in " << sw.printTime() << " seconds." << endl;

		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int itrig = 0; itrig < ntrig; ++itrig)
			full_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, itrig)]
				= thermal_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, itrig)];
	}
	else
	{
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int itrig = 0; itrig < ntrig; ++itrig)
			full_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, itrig)]
				= 0.0;
	}

	bool include_resonance_pions = true;
	if (include_resonance_pions)
	{
		cout << "Setting full target moments..." << endl;
		//if not just doing thermal pions, include resonance contributions, too
		if ( !thermal_pions_only )
			Set_full_target_moments(iqt, iqz);
	}
	cout << "done." << endl;

	return;
}

//assume only want midrapidity (Y=0) pions for the timebeing
void CorrelationFunction::Set_thermal_target_moments(int iqt, int iqz)
{
	double loc_qz = qz_pts[iqz];
	double loc_qt = qt_pts[iqt];
	current_pY_shift = 0.5 * log(abs((loc_qt+loc_qz + 1.e-100)/(loc_qt-loc_qz + 1.e-100)));

	for (int ipY = 0; ipY < n_pY_pts; ++ipY)
	{
		ch_SP_pY[ipY] = cosh(SP_Del_pY[ipY] + current_pY_shift);
		sh_SP_pY[ipY] = sinh(SP_Del_pY[ipY] + current_pY_shift);
	}

	if ( USE_EXACT && !USE_CF )
		Cal_dN_dypTdpTdphi_with_weights_toy(target_particle_id, iqt, iqz, ipY0, thermal_target_Yeq0_moments);
	else
		Cal_dN_dypTdpTdphi_with_weights_Yeq0_alternate(iqt, iqz);

	return;
}

void CorrelationFunction::Set_full_target_moments(int iqt, int iqz)
{
	//this one includes resonance decay contributions,
	//needs thermal contributions added in
	int getHDFresonanceSpectra = Access_target_full_in_HDF_array(iqt, iqz, 1, current_dN_dypTdpTdphi_moments, true);

	bool use_interp = true;
	if (use_interp)
	{
		double loc_qz = qz_pts[iqz];
		double loc_qt = qt_pts[iqt];
		current_pY_shift = 0.5 * log(abs((loc_qt+loc_qz + 1.e-100)/(loc_qt-loc_qz + 1.e-100)));
		if (abs(current_pY_shift) > SP_Del_pY_max)	//if resonance decays basically contribute nothing at this (qt, qz)-pair
		{
			cout << "Set_full_target_moments(" << loc_qt << ", " << loc_qz << "): y_{sym} = " << current_pY_shift << " > SP_Del_pY_max = " << SP_Del_pY_max <<
					" ==>> No contributions to full target moments." << endl;
			return;	//don't do anything
		}

		gsl_cheb_series *cs_accel_expEdNd3p = gsl_cheb_alloc (n_pY_pts - 1);
		cs_accel_expEdNd3p->a = adjusted_SP_Del_pY_minimum;
		cs_accel_expEdNd3p->b = SP_Del_pY_max;

		double * chebyshev_a_cfs = new double[n_pY_pts];

		double trig_parity_factors[4] = {1.0, -1.0, -1.0, 1.0};

		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int itrig = 0; itrig < ntrig; ++itrig)
		{
			for (int ipY = 0; ipY < n_pY_pts; ++ipY)
			{
				chebyshev_a_cfs[ipY] = 0.0;
				for (int kpY = 0; kpY < n_pY_pts; ++kpY)
				{
					chebyshev_a_cfs[ipY] += exp(abs(SP_Del_pY[kpY])) * chebTcfs[ipY * n_pY_pts + kpY] * current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,kpY,iqx,iqy,itrig)];

					if (ipY==0 && ipT==0 && ipphi==0) cout << "CHECKinterp: " << iqz << "   " << itrig << "   " << SP_Del_pY[kpY]
						<< "   " << thermal_target_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,kpY,iqx,iqy,itrig)]
						<< "   " << current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,kpY,iqx,iqy,itrig)] << endl;
				}
			}

			cs_accel_expEdNd3p->c = chebyshev_a_cfs;
			double tmp_Del_pY = -current_pY_shift;	//interpolating to this point
			double current_trig_parity_factor = ( ( tmp_Del_pY < 0.0 )
												&& ( abs(loc_qt) < abs(loc_qz) )
												) ? trig_parity_factors[itrig] : 1.0;
			full_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, itrig)] += exp(-abs(tmp_Del_pY)) * current_trig_parity_factor
																						* gsl_cheb_eval (cs_accel_expEdNd3p, abs(tmp_Del_pY));
			if (ipT==0 && ipphi==0) cout << "CHECKinterp2: " << iqz << "   " << itrig << "   "
								<< exp(-abs(tmp_Del_pY)) * current_trig_parity_factor
										* gsl_cheb_eval (cs_accel_expEdNd3p, abs(tmp_Del_pY)) << endl;
		}

		delete [] chebyshev_a_cfs;
	}
	else
	{
		//alt: skip interp if only doing Y=0
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int itrig = 0; itrig < ntrig; ++itrig)
			full_target_Yeq0_moments[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, itrig)] += current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY0,iqx,iqy,itrig)];
	}

	return;
}


//*********************************************************************
// 3D case
//*********************************************************************
//  Simple function to print results of each iteration in nice format
int CorrelationFunction::print_fit_state_3D (size_t iteration, gsl_multifit_fdfsolver * solver_ptr)
{
	cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
	cout.precision (5);		// digits in doubles

	int width = 15;		// setw width for output
	cout << scientific
		<< "iteration " << iteration << ": "
		<< "  x = {" << setw (width) << gsl_vector_get (solver_ptr->x, 0)
		<< setw (width) << gsl_vector_get (solver_ptr->x, 1)
		<< setw (width) << gsl_vector_get (solver_ptr->x, 2)
		<< setw (width) << gsl_vector_get (solver_ptr->x, 3)
		<< "}, |f(x)| = " << scientific << gsl_blas_dnrm2 (solver_ptr->f) 
		<< endl << endl;

	return 0;
}
//  Simple function to print results of each iteration in nice format
int CorrelationFunction::print_fit_state_3D_withlambda (size_t iteration, gsl_multifit_fdfsolver * solver_ptr)
{
	cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
	cout.precision (5);		// digits in doubles

	int width = 15;		// setw width for output
	cout << scientific
		<< "iteration " << iteration << ": "
		<< "  x = {" << setw (width) << gsl_vector_get (solver_ptr->x, 0)
		<< setw (width) << gsl_vector_get (solver_ptr->x, 1)
		<< setw (width) << gsl_vector_get (solver_ptr->x, 2)
		<< setw (width) << gsl_vector_get (solver_ptr->x, 3)
		<< setw (width) << gsl_vector_get (solver_ptr->x, 4)
		<< "}, |f(x)| = " << scientific << gsl_blas_dnrm2 (solver_ptr->f) 
		<< endl << endl;

	return 0;
}

//*********************************************************************
//  Function returning the residuals for each point; that is, the 
//  difference of the fit function using the current parameters
//  and the data to be fit.
int Fittarget_correlfun3D_f (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr)
{
	size_t n = ((struct Correlationfunction3D_data *) params_ptr)->data_length;
	double * q_o = ((struct Correlationfunction3D_data *) params_ptr)->q_o;
	double * q_s = ((struct Correlationfunction3D_data *) params_ptr)->q_s;
	double * q_l = ((struct Correlationfunction3D_data *) params_ptr)->q_l;
	double * y = ((struct Correlationfunction3D_data *) params_ptr)->y;
	double * sigma = ((struct Correlationfunction3D_data *) params_ptr)->sigma;

	//fit parameters
	double R2_o = gsl_vector_get (xvec_ptr, 0);
	double R2_s = gsl_vector_get (xvec_ptr, 1);
	double R2_l = gsl_vector_get (xvec_ptr, 2);
	double R2_os = gsl_vector_get (xvec_ptr, 3);

	size_t i;

	for (i = 0; i < n; i++)
	{
		double Yi = 1.0 + exp(- q_l[i]*q_l[i]*R2_l - q_s[i]*q_s[i]*R2_s - q_o[i]*q_o[i]*R2_o - 2.*q_o[i]*q_s[i]*R2_os);
		gsl_vector_set (f_ptr, i, (Yi - y[i]) / sigma[i]);
//cout << "i = " << i << ": " << y[i] << "   " << Yi << "   " << (Yi - y[i]) / sigma[i]
//		<< "   " << q_o[i] << "   " << q_s[i] << "   " << q_l[i] << "   " << R2_o << "   " << R2_s << "   " << R2_l << "   " << R2_os << endl;
	}

	return GSL_SUCCESS;
}

int Fittarget_correlfun3D_f_withlambda (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr)
{
	size_t n = ((struct Correlationfunction3D_data *) params_ptr)->data_length;
	double * q_o = ((struct Correlationfunction3D_data *) params_ptr)->q_o;
	double * q_s = ((struct Correlationfunction3D_data *) params_ptr)->q_s;
	double * q_l = ((struct Correlationfunction3D_data *) params_ptr)->q_l;
	double * y = ((struct Correlationfunction3D_data *) params_ptr)->y;
	double * sigma = ((struct Correlationfunction3D_data *) params_ptr)->sigma;

	//fit parameters
	double lambda = gsl_vector_get (xvec_ptr, 0);
	double R2_o = gsl_vector_get (xvec_ptr, 1);
	double R2_s = gsl_vector_get (xvec_ptr, 2);
	double R2_l = gsl_vector_get (xvec_ptr, 3);
	double R2_os = gsl_vector_get (xvec_ptr, 4);

	size_t i;

	for (i = 0; i < n; i++)
	{
		//double Yi = lambda*exp(- q_l[i]*q_l[i]*R_l*R_l - q_s[i]*q_s[i]*R_s*R_s
		//             - q_o[i]*q_o[i]*R_o*R_o - q_o[i]*q_s[i]*R_os*R_os);
		double Yi = 1.0 + lambda*exp(- q_l[i]*q_l[i]*R2_l - q_s[i]*q_s[i]*R2_s - q_o[i]*q_o[i]*R2_o - 2.*q_o[i]*q_s[i]*R2_os);
		gsl_vector_set (f_ptr, i, (Yi - y[i]) / sigma[i]);
//cout << "i = " << i << ": " << y[i] << "   " << Yi << "   " << (Yi - y[i]) / sigma[i] << "   " << lambda
//		<< "   " << q_o[i] << "   " << q_s[i] << "   " << q_l[i] << "   " << R2_o << "   " << R2_s << "   " << R2_l << "   " << R2_os << endl;
	}

	return GSL_SUCCESS;
}

//*********************************************************************
//  Function returning the Jacobian of the residual function
int Fittarget_correlfun3D_df (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr)
{
	size_t n = ((struct Correlationfunction3D_data *) params_ptr)->data_length;
	double * q_o = ((struct Correlationfunction3D_data *) params_ptr)->q_o;
	double * q_s = ((struct Correlationfunction3D_data *) params_ptr)->q_s;
	double * q_l = ((struct Correlationfunction3D_data *) params_ptr)->q_l;
	double * sigma = ((struct Correlationfunction3D_data *) params_ptr)->sigma;

	//fit parameters
	double R2_o = gsl_vector_get (xvec_ptr, 0);
	double R2_s = gsl_vector_get (xvec_ptr, 1);
	double R2_l = gsl_vector_get (xvec_ptr, 2);
	double R2_os = gsl_vector_get (xvec_ptr, 3);

	size_t i;

	for (i = 0; i < n; i++)
	{
		double sig = sigma[i];

		// derivatives
		// double common_elemt = exp(- q_l[i]*q_l[i]*R_l*R_l - q_s[i]*q_s[i]*R_s*R_s - q_o[i]*q_o[i]*R_o*R_o - q_o[i]*q_s[i]*R_os*R_os);
		double common_elemt = exp(- q_l[i]*q_l[i]*R2_l - q_s[i]*q_s[i]*R2_s - q_o[i]*q_o[i]*R2_o - 2.*q_o[i]*q_s[i]*R2_os);
      
		gsl_matrix_set (Jacobian_ptr, i, 0, - q_o[i]*q_o[i]*common_elemt/sig);
		gsl_matrix_set (Jacobian_ptr, i, 1, - q_s[i]*q_s[i]*common_elemt/sig);
		gsl_matrix_set (Jacobian_ptr, i, 2, - q_l[i]*q_l[i]*common_elemt/sig);
		gsl_matrix_set (Jacobian_ptr, i, 3, - 2.*q_o[i]*q_s[i]*common_elemt/sig);
	}

	return GSL_SUCCESS;
}

int Fittarget_correlfun3D_df_withlambda (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr)
{
	size_t n = ((struct Correlationfunction3D_data *) params_ptr)->data_length;
	double * q_o = ((struct Correlationfunction3D_data *) params_ptr)->q_o;
	double * q_s = ((struct Correlationfunction3D_data *) params_ptr)->q_s;
	double * q_l = ((struct Correlationfunction3D_data *) params_ptr)->q_l;
	double * sigma = ((struct Correlationfunction3D_data *) params_ptr)->sigma;

	//fit parameters
	double lambda = gsl_vector_get (xvec_ptr, 0);
	double R2_o = gsl_vector_get (xvec_ptr, 1);
	double R2_s = gsl_vector_get (xvec_ptr, 2);
	double R2_l = gsl_vector_get (xvec_ptr, 3);
	double R2_os = gsl_vector_get (xvec_ptr, 4);

	size_t i;

	for (i = 0; i < n; i++)
	{
		double sig = sigma[i];

		//derivatives
		double common_elemt = exp(- q_l[i]*q_l[i]*R2_l - q_s[i]*q_s[i]*R2_s - q_o[i]*q_o[i]*R2_o - 2.*q_o[i]*q_s[i]*R2_os);
      
		gsl_matrix_set (Jacobian_ptr, i, 0, common_elemt/sig);
		gsl_matrix_set (Jacobian_ptr, i, 1, - lambda*q_o[i]*q_o[i]*common_elemt/sig);
		gsl_matrix_set (Jacobian_ptr, i, 2, - lambda*q_s[i]*q_s[i]*common_elemt/sig);
		gsl_matrix_set (Jacobian_ptr, i, 3, - lambda*q_l[i]*q_l[i]*common_elemt/sig);
		gsl_matrix_set (Jacobian_ptr, i, 4, - 2.*lambda*q_o[i]*q_s[i]*common_elemt/sig);
	}

	return GSL_SUCCESS;
}

//*********************************************************************
//  Function combining the residual function and its Jacobian
int Fittarget_correlfun3D_fdf (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr)
{
	Fittarget_correlfun3D_f(xvec_ptr, params_ptr, f_ptr);
	Fittarget_correlfun3D_df(xvec_ptr, params_ptr, Jacobian_ptr);

	return GSL_SUCCESS;
}

int Fittarget_correlfun3D_fdf_withlambda (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr)
{
	Fittarget_correlfun3D_f_withlambda(xvec_ptr, params_ptr, f_ptr);
	Fittarget_correlfun3D_df_withlambda(xvec_ptr, params_ptr, Jacobian_ptr);

	return GSL_SUCCESS;
}

//End of file
