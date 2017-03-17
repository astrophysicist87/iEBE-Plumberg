#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>

#include "fitCF.h"
#include "fitCF_lib.h"
#include "CPStopwatch.h"
#include "Stopwatch.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"
#include "chebyshev.h"
#include "stats.h"

using namespace std;

void FitCF::Get_GF_HBTradii()
{
	//if (!VARY_ALPHA)
		*global_out_stream_ptr << "--> Getting HBT radii by Gaussian fit method" << endl;
	//else
	//	*global_out_stream_ptr << "--> Getting HBT radii by Levy-stable fit method" << endl;

	if (FLESH_OUT_CF)
		Allocate_fleshed_out_CF();

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		*global_out_stream_ptr << "   --> Doing pT = " << SPinterp_pT[ipt] << ", pphi = " << SPinterp_pphi[ipphi] << "..." << endl;
		
		//determine whether to use fleshed out / projected CFvals
		double *** CF_for_fitting = CFvals[ipt][ipphi];
		if (FLESH_OUT_CF)
		{
			Flesh_out_CF(ipt, ipphi);
			CF_for_fitting = fleshed_out_CF;
		}

		//finally, do fits, depending on what kind you want to do
		if (USE_LAMBDA)
			Fit_Correlationfunction3D_withlambda( CF_for_fitting, ipt, ipphi, FLESH_OUT_CF );
		else
			Fit_Correlationfunction3D( CF_for_fitting, ipt, ipphi, FLESH_OUT_CF );
	}

	if (FLESH_OUT_CF)
		Delete_fleshed_out_CF();

	return;
}

//to save time, don't bother making grid large enough to interpolate to OSL
//this function leaves open the option of just interpolating over the qt-direction
void FitCF::Cal_correlationfunction()
{
	*global_out_stream_ptr << "Calculating the correlation function..." << endl;

	// Can't interpolate if there's only one point in qt-direction!
	if (qtnpts == 1)
		return;

	string local_name = all_particles[target_particle_id].name;
	replace_parentheses(local_name);

	ostringstream filename_stream_target_dN_dypTdpTdphi_evavg;
	filename_stream_target_dN_dypTdpTdphi_evavg << global_path
				<< "/total_" << local_name << "_eiqx_dN_dypTdpTdphi_evavg" << no_df_stem << ".dat";

	//now check if file exists
	bool avg_file_exists = fexists(filename_stream_target_dN_dypTdpTdphi_evavg.str().c_str());

	if (avg_file_exists)
	{
		Readin_total_target_eiqx_dN_dypTdpTdphi_evavg();
		Readin_resonance_fraction(1);
	}
	else
	{
		////////////////////////////////////////////////////////////
		//reading in events...
		for (int iEvent = 0; iEvent < nEvents; ++iEvent)
		{
			*global_out_stream_ptr << "Reading in event = " << chosen_events[iEvent] << endl;

			//construct needed filename
			string resultsDirectory = "";
			if (currentfolderindex == -1)
				resultsDirectory = "/results-" + patch::to_string(chosen_events[iEvent]);
			ostringstream filename_stream_target_dN_dypTdpTdphi;
			filename_stream_target_dN_dypTdpTdphi << global_path << resultsDirectory
						<< "/total_" << local_name << "_eiqx_dN_dypTdpTdphi" << no_df_stem << ".dat";
						//<< "/total_" << local_name << "_eiqx_dN_dypTdpTdphi_ev" << chosen_events[iEvent] << no_df_stem << ".dat";

			//now check if file exists
			bool file_exists = fexists(filename_stream_target_dN_dypTdpTdphi.str().c_str());
			if (file_exists)
				*global_out_stream_ptr << "    --> Note: file exists." << endl;
			else
				*global_out_stream_ptr << "    --> Note: file does not exist." << endl;

			//if not, unzip from zip file first, then remove when done
			if (not file_exists)
			{
				ostringstream local_cmd;
				//need to change directories if actually in a results directory...
				if (currentfolderindex != -1)
					local_cmd << "cd ..; ";
				//...unzip...
				//local_cmd << "unzip FTresults.zip results-" << chosen_events[iEvent]
				//			<< "/total_" << local_name << "_eiqx_dN_dypTdpTdphi_ev" << chosen_events[iEvent] << no_df_stem << ".dat; ";
				local_cmd << "unzip FTspectra_ev" << chosen_events[iEvent] << ".zip; ";
				//...then change back
				if (currentfolderindex != -1)
					local_cmd << "cd " << global_path;

				//now submit the full command
				*global_out_stream_ptr << "    --> Running this command: " << local_cmd.str().c_str() << endl;
				int cmd_result = system(local_cmd.str().c_str());
			}

			//now read the file in and average appropriately
			Readin_resonance_fraction(chosen_events[iEvent]);
			//Readin_total_target_eiqx_dN_dypTdpTdphi(chosen_events[iEvent]);
			Readin_total_target_eiqx_dN_dypTdpTdphi(filename_stream_target_dN_dypTdpTdphi.str());

			//delete file, as promised
			if (not file_exists)
			{
				ostringstream local_cmd;
				local_cmd << "/bin/rm " << filename_stream_target_dN_dypTdpTdphi.str();
				*global_out_stream_ptr << "    --> Running this command: " << local_cmd.str().c_str() << endl;
				int cmd_result = system(local_cmd.str().c_str());
			}

		}

		Average_total_target_eiqx_dN_dypTdpTdphi();	//if nEvents=1, does nothing
	}
	////////////////////////////////////////////////////////////////////////

	//now that we have the spectra (averaged or not), the rest of the calculation is the same

	// chooses the qo, qs, ql (or qx, qy, ql) points at which to evaluate correlation function,
	// and allocates the array to hold correlation function values
	Set_correlation_function_q_pts();
	Allocate_CFvals();

	double * q_interp = new double [4];

	// Then compute full correlation function
	*global_out_stream_ptr << "Computing the full correlator in XYZ coordinates..." << endl;
	CPStopwatch sw;
	sw.Start();
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	{
		Get_q_points(qx_pts[iqx], qy_pts[iqy], qz_pts[iqz], SPinterp_pT[ipt], SPinterp_pphi[ipphi], q_interp);

		//returns only projected value automatically if appropriate options are specified!
		double tmp1 = 0.0, tmp2 = 0.0, tmp2a = 0.0, tmp3 = 0.0;
		Compute_correlationfunction(&tmp1, &tmp2, &tmp2a, &tmp3, ipt, ipphi, iqx, iqy, iqz, q_interp[0], 0);
		CFvals[ipt][ipphi][iqx][iqy][iqz] = 1.0 + tmp1;		//C == Ct + Cct + Cr + 1
		thermalCFvals[ipt][ipphi][iqx][iqy][iqz] = tmp2;	//Ct
		crosstermCFvals[ipt][ipphi][iqx][iqy][iqz] = tmp2a;	//Cct
		resonancesCFvals[ipt][ipphi][iqx][iqy][iqz] = tmp3;	//Cr
	}
	sw.Stop();
	*global_out_stream_ptr << "Finished computing correlator in " << sw.printTime() << " seconds." << endl;

	//output the un-regulated correlation function to separate file for debugging purposes
	Output_correlationfunction();

	if (REGULATE_CF)
	{
		sw.Reset();
		sw.Start();
	
		*global_out_stream_ptr << "WARNING: Regulating computed correlator values using Hampel criterion for outlier detection..." << endl;
		double * pphi_CF_slice = new double [n_interp_pphi_pts];
		double * pphi_CF_slice_term1 = new double [n_interp_pphi_pts];
		double * pphi_CF_slice_term2 = new double [n_interp_pphi_pts];
		double * pphi_CF_slice_term3 = new double [n_interp_pphi_pts];
	
		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int iqz = 0; iqz < qznpts; ++iqz)
		{
			for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
			{
				pphi_CF_slice[ipphi] = CFvals[ipt][ipphi][iqx][iqy][iqz];
				pphi_CF_slice_term1[ipphi] = thermalCFvals[ipt][ipphi][iqx][iqy][iqz];
				pphi_CF_slice_term2[ipphi] = crosstermCFvals[ipt][ipphi][iqx][iqy][iqz];
				pphi_CF_slice_term3[ipphi] = resonancesCFvals[ipt][ipphi][iqx][iqy][iqz];
			}
			//Regulate_CF_Hampel(ipt, iqx, iqy, iqz, pphi_CF_slice, pphi_CF_slice_term1, pphi_CF_slice_term2, pphi_CF_slice_term3);
			Regulate_CF_Hampel_v2(ipt, iqx, iqy, iqz, pphi_CF_slice, pphi_CF_slice_term1, pphi_CF_slice_term2, pphi_CF_slice_term3);
		}
		sw.Stop();
		*global_out_stream_ptr << "WARNING: Finished regulating computed correlator values in " << sw.printTime() << " seconds." << endl;
	
		delete [] pphi_CF_slice;
		delete [] pphi_CF_slice_term1;
		delete [] pphi_CF_slice_term2;
		delete [] pphi_CF_slice_term3;
	}

	delete [] q_interp;

	return;
}

void FitCF::Get_total_target_eiqx_dN_dypTdpTdphi_on_pair_momentum_grid(double * eiqx_EdNd3p_in, double * eiqx_EdNd3p_out, double * KT, double * KPHI)
{
	const int arrayLength = n_interp_pT_pts * n_interp_pphi_pts * qtnpts * qxnpts * qynpts * qznpts * ntrig;

	double ** localGrid = new double * [n_interp_pT_pts];
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		localGrid[ipt] = new double [n_interp_pphi_pts];
	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	for (int itrig = 0; itrig < ntrig; ++itrig)
	{
		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
			localGrid[ipt][ipphi] = eiqx_EdNd3p_in[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, itrig)];

		for (int iKT = 0; iKT < n_localp_T; ++iKT)
		for (int iKPHI = 0; iKPHI < n_localp_phi; ++iKPHI)
			eiqx_EdNd3p_out[indexer2(iKT, iKPHI, iqt, iqx, iqy, iqz, itrig)]
				= interpolate2D(SPinterp_pT, SPinterp_pphi, localGrid, KT[iKT], KPHI[iKPHI], n_interp_pT_pts, n_interp_pphi_pts, 0, false);
	}

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		delete [] localGrid[ipt];
	delete [] localGrid;

	return;
}

void FitCF::Compute_correlationfunction(double * totalresult, double * thermalresult, double * CTresult, double * resonanceresult,
										int ipt, int ipphi, int iqx, int iqy, int iqz, double qt_interp, int interp_flag /*==0*/)
{
	int qidx = binarySearch(qt_pts, qtnpts, qt_interp);
	double q_min = qt_pts[0] / cos(M_PI / (2.*qtnpts)), q_max = qt_pts[qtnpts-1]/ cos(M_PI / (2.*qtnpts));

	bool q_point_is_outside_grid = ( qidx == -1 && ( qt_interp < q_min || qt_interp > q_max ) );

	if (!q_point_is_outside_grid)
	{
		double C_at_q[qtnpts], Ct_at_q[qtnpts], Cct_at_q[qtnpts], Cr_at_q[qtnpts];	//C - 1
		double tmpC = 0.0, tmpCt = 0.0, tmpCct = 0.0, tmpCr = 0.0;
	
		// set CF values along qt-slice for interpolation
		for (int iqtidx = 0; iqtidx < qtnpts; ++iqtidx)
		{
			//return C - 1!!!
			get_CF(&tmpC, &tmpCt, &tmpCct, &tmpCr, ipt, ipphi, iqtidx, iqx, iqy, iqz, FIT_WITH_PROJECTED_CFVALS && !thermal_pions_only);
			C_at_q[iqtidx] = tmpC;
			Ct_at_q[iqtidx] = tmpCt;
			Cct_at_q[iqtidx] = tmpCct;
			Cr_at_q[iqtidx] = tmpCr;
//cout << "Working with " << tmpC << "   " << tmpCt << "   " << tmpCct << "   " << tmpCr << "   " << FIT_WITH_PROJECTED_CFVALS << "   " << thermal_pions_only << endl;
		}

		//assumes qt-grid has already been computed at (adjusted) Chebyshev nodes!!!
		if (QT_POINTS_SPACING == 1 && interp_flag == 0)
		{
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

			//*global_out_stream_ptr << "CHECK: " << *totalresult << "   " << *thermalresult << "   " << *CTresult << "   " << *resonanceresult << endl;
		}
		else	//if not using Chebyshev nodes in qt-direction, just use straight-up linear(0) or cubic(1) interpolation
		{
			*totalresult = interpolate1D(qt_pts, C_at_q, qt_interp, qtnpts, 1, false);
			*thermalresult = interpolate1D(qt_pts, Ct_at_q, qt_interp, qtnpts, 1, false);
			*CTresult = interpolate1D(qt_pts, Cct_at_q, qt_interp, qtnpts, 1, false);
			*resonanceresult = interpolate1D(qt_pts, Cr_at_q, qt_interp, qtnpts, 1, false);
		}
	}
	else
	{
		*global_out_stream_ptr << "Warning: qt_interp point was outside of computed grid!" << endl
								<< "\t qt_interp = " << qt_interp << " out of {q_min, q_max} = {" << q_min << ", " << q_max << "}" << endl;
		*totalresult = 0.0;
	}
	return;
}


void FitCF::Average_total_target_eiqx_dN_dypTdpTdphi()
{
	*global_out_stream_ptr << "Ensemble-averaging single-particle and two-particle spectra..." << endl;

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		spectra[target_particle_id][ipt][ipphi] /= double(nEvents);
		thermal_spectra[target_particle_id][ipt][ipphi] /= double(nEvents);
	}

	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		//full_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0] /= double(nEvents);
		//full_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1] /= double(nEvents);
		//thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0] /= double(nEvents);
		//thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1] /= double(nEvents);
		full_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 0)] /= double(nEvents);
		full_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 1)] /= double(nEvents);
		thermal_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 0)] /= double(nEvents);
		thermal_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 1)] /= double(nEvents);
	}

	*global_out_stream_ptr << "Finished ensemble-averaging single-particle and two-particle spectra." << endl;
	return;
}

//******************************************************************
// Routines for refining correlation function grid via interpolation
//******************************************************************

void FitCF::Flesh_out_CF(int ipt, int ipphi, double sample_scale /*==1.0*/)
{
	//declare needed quantities here
	double qxmin = 0.9999*sample_scale*qx_pts[0], qxmax = 0.9999*sample_scale*qx_pts[qxnpts-1];
	double qymin = 0.9999*sample_scale*qy_pts[0], qymax = 0.9999*sample_scale*qy_pts[qynpts-1];
	double qzmin = 0.9999*sample_scale*qz_pts[0], qzmax = 0.9999*sample_scale*qz_pts[qznpts-1];

	double new_Del_qx = (qxmax - qxmin)/(double(new_nqpts-1)+1.e-100);
	double new_Del_qy = (qymax - qymin)/(double(new_nqpts-1)+1.e-100);
	double new_Del_qz = (qzmax - qzmin)/(double(new_nqpts-1)+1.e-100);

	double *** current_C_slice = CFvals[ipt][ipphi];

	//cout << "(qxmin, qxmax, new_Del_qx) = (" << qxmin << ", " << qxmax << ", " << new_Del_qx << ")" << endl;
	//cout << "(qymin, qymax, new_Del_qy) = (" << qymin << ", " << qymax << ", " << new_Del_qy << ")" << endl;
	//cout << "(qzmin, qzmax, new_Del_qz) = (" << qzmin << ", " << qzmax << ", " << new_Del_qz << ")" << endl;
	
	//if we chose a completely Chebyshev grid, exploit this...
	if (QX_POINTS_SPACING && QY_POINTS_SPACING && QZ_POINTS_SPACING)
	{
		double qx_rescale = 1./(cos(M_PI / (2.*qxnpts)));
		double qy_rescale = 1./(cos(M_PI / (2.*qynpts)));
		double qz_rescale = 1./(cos(M_PI / (2.*qznpts)));

		/////////////////////////////////////////////////////////////////////////////////////////////////
		// INITIALIZE CHEBYSHEV STUFF HERE
		/////////////////////////////////////////////////////////////////////////////////////////////////
		//for 3D Chebyshev
		//X-Y-Z
		const int dim_loc = 3;
		int npts_loc[3] = { qxnpts, qynpts, qznpts };
		int os[3] = { qxnpts - 1, qynpts - 1, qznpts - 1 };
		double lls[3] = { qx_rescale*qx_pts[0], qy_rescale*qy_pts[0], qz_rescale*qz_pts[0] };
		double uls[3] = { qx_rescale*qx_pts[qxnpts - 1], qy_rescale*qy_pts[qynpts - 1], qz_rescale*qz_pts[qznpts - 1] };

		//different combinations of directions for 2D Chebyshev
		//X-Y
		int npts_loc_XY[2] = { qxnpts, qynpts };
		int os_XY[2] = { qxnpts - 1, qynpts - 1 };
		double lls_XY[2] = { qx_rescale*qx_pts[0], qy_rescale*qy_pts[0] };
		double uls_XY[2] = { qx_rescale*qx_pts[qxnpts - 1], qy_rescale*qy_pts[qynpts - 1] };
		//X-Z
		int npts_loc_XZ[2] = { qxnpts, qznpts };
		int os_XZ[2] = { qxnpts - 1, qznpts - 1 };
		double lls_XZ[2] = { qx_rescale*qx_pts[0], qz_rescale*qz_pts[0] };
		double uls_XZ[2] = { qx_rescale*qx_pts[qxnpts - 1], qz_rescale*qz_pts[qznpts - 1] };
		//Y-Z
		int npts_loc_YZ[2] = { qynpts, qznpts };
		int os_YZ[2] = { qynpts - 1, qznpts - 1 };
		double lls_YZ[2] = { qy_rescale*qy_pts[0], qz_rescale*qz_pts[0] };
		double uls_YZ[2] = { qy_rescale*qy_pts[qynpts - 1], qz_rescale*qz_pts[qznpts - 1] };
		//different combinations of directions for 1D Chebyshev
		//X
		int npts_loc_X[1] = { qxnpts };
		int os_X[1] = { qxnpts - 1 };
		double lls_X[1] = { qx_rescale*qx_pts[0] };
		double uls_X[1] = { qx_rescale*qx_pts[qxnpts - 1] };
		//Y
		int npts_loc_Y[1] = { qynpts };
		int os_Y[1] = { qynpts - 1 };
		double lls_Y[1] = { qy_rescale*qy_pts[0] };
		double uls_Y[1] = { qy_rescale*qy_pts[qynpts - 1] };
		//Z
		int npts_loc_Z[1] = { qznpts };
		int os_Z[1] = { qznpts - 1 };
		double lls_Z[1] = { qz_rescale*qz_pts[0] };
		double uls_Z[1] = { qz_rescale*qz_pts[qznpts - 1] };

		double flat_C_at_q[qxnpts*qynpts*qznpts];
		double flat2D_C_at_q_XY[qxnpts*qynpts], flat2D_C_at_q_XZ[qxnpts*qznpts], flat2D_C_at_q_YZ[qynpts*qznpts];
		double flat1D_C_at_q_X[qxnpts], flat1D_C_at_q_Y[qynpts], flat1D_C_at_q_Z[qznpts];

		int qidx = 0;
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int iqz = 0; iqz < qznpts; ++iqz)
		{
			flat_C_at_q[qidx++] = current_C_slice[iqx][iqy][iqz] - 1.0;
			cout << "FC: " << current_C_slice[iqx][iqy][iqz] << endl;
		}

		Chebyshev cf(flat_C_at_q, npts_loc, os, lls, uls, 3);

		vector<Chebyshev*> cf_XY;
		for (int iqz = 0; iqz < qznpts; ++iqz)
		{
			int iqXY = 0;
			for (int iqx = 0; iqx < qxnpts; ++iqx)
			for (int iqy = 0; iqy < qynpts; ++iqy)
				flat2D_C_at_q_XY[iqXY++] = current_C_slice[iqx][iqy][iqz] - 1.0;
			cf_XY.push_back(new Chebyshev(flat2D_C_at_q_XY, npts_loc_XY, os_XY, lls_XY, uls_XY, 2) );
		}
		
		vector<Chebyshev*> cf_XZ;
		for (int iqy = 0; iqy < qynpts; ++iqy)
		{
			int iqXZ = 0;
			for (int iqx = 0; iqx < qxnpts; ++iqx)
			for (int iqz = 0; iqz < qznpts; ++iqz)
				flat2D_C_at_q_XZ[iqXZ++] = current_C_slice[iqx][iqy][iqz] - 1.0;
			cf_XZ.push_back(new Chebyshev(flat2D_C_at_q_XZ, npts_loc_XZ, os_XZ, lls_XZ, uls_XZ, 2) );
		}

		vector<Chebyshev*> cf_YZ;
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		{
			int iqYZ = 0;
			for (int iqy = 0; iqy < qynpts; ++iqy)
			for (int iqz = 0; iqz < qznpts; ++iqz)
				flat2D_C_at_q_YZ[iqYZ++] = current_C_slice[iqx][iqy][iqz] - 1.0;
			cf_YZ.push_back(new Chebyshev(flat2D_C_at_q_YZ, npts_loc_YZ, os_YZ, lls_YZ, uls_YZ, 2) );
		}

		vector< vector<Chebyshev*> > cf_X (qynpts);
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int iqz = 0; iqz < qznpts; ++iqz)
		{
			int iqX = 0;
			for (int iqx = 0; iqx < qxnpts; ++iqx)
				flat1D_C_at_q_X[iqX++] = current_C_slice[iqx][iqy][iqz] - 1.0;
			cf_X[iqy].push_back(new Chebyshev(flat1D_C_at_q_X, npts_loc_X, os_X, lls_X, uls_X, 1) );
		}
		
		vector< vector<Chebyshev*> > cf_Y (qxnpts);
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqz = 0; iqz < qznpts; ++iqz)
		{
			int iqY = 0;
			for (int iqy = 0; iqy < qynpts; ++iqy)
				flat1D_C_at_q_Y[iqY++] = current_C_slice[iqx][iqy][iqz] - 1.0;
			cf_Y[iqx].push_back(new Chebyshev(flat1D_C_at_q_Y, npts_loc_Y, os_Y, lls_Y, uls_Y, 1) );
		}

		vector< vector<Chebyshev*> > cf_Z (qxnpts);
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		{
			int iqZ = 0;
			for (int iqz = 0; iqz < qznpts; ++iqz)
				flat1D_C_at_q_Z[iqZ++] = current_C_slice[iqx][iqy][iqz] - 1.0;
			cf_Z[iqx].push_back(new Chebyshev(flat1D_C_at_q_Z, npts_loc_Z, os_Z, lls_Z, uls_Z, 1) );
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////
		// END ALL CHEBYSHEV INITIALIZATION
		/////////////////////////////////////////////////////////////////////////////////////////////////

		//first, enumerate various regions by region index
		int region_index[new_nqpts][new_nqpts][new_nqpts];
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

			if (abs(qx0) <= qxmax)
			{
				if (abs(qy0) <= qymax)
				{
					if (abs(qz0) <= qzmax)
						region_index[iqx][iqy][iqz] = 0;	//XYZ-Cheb, no Gauss
					else
						region_index[iqx][iqy][iqz] = 1;	//XY-Cheb, Z-Gauss
				}
				else
				{
					if (abs(qz0) <= qzmax)
						region_index[iqx][iqy][iqz] = 2;	//XZ-Cheb, Y-Gauss
					else
						region_index[iqx][iqy][iqz] = 3;	//X-Cheb, YZ-Gauss
				}
			}
			else
			{
				if (abs(qy0) <= qymax)
				{
					if (abs(qz0) <= qzmax)
						region_index[iqx][iqy][iqz] = 4;	//YZ-Cheb, X-Gauss
					else
						region_index[iqx][iqy][iqz] = 5;	//Y-Cheb, XZ-Gauss
				}
				else
				{
					if (abs(qz0) <= qzmax)
						region_index[iqx][iqy][iqz] = 6;	//Z-Cheb, XY-Gauss
					else
						region_index[iqx][iqy][iqz] = 7;	//No Cheb, XYZ-Gauss
				}
			}
		}

		//next, use switch statement on region_index to decide how to interpolate
		for (int iqx = 0; iqx < new_nqpts; ++iqx)
		for (int iqy = 0; iqy < new_nqpts; ++iqy)
		for (int iqz = 0; iqz < new_nqpts; ++iqz)
		{
			double qx0 = qx_fleshed_out_pts[iqx];
			double qy0 = qy_fleshed_out_pts[iqy];
			double qz0 = qz_fleshed_out_pts[iqz];
			int iqx0_loc = ( qx0 < qx_pts[0]) ? 0 : qxnpts - 2;
			int iqy0_loc = ( qy0 < qy_pts[0]) ? 0 : qynpts - 2;
			int iqz0_loc = ( qz0 < qz_pts[0]) ? 0 : qznpts - 2;

cout << "FLESH CHECK: " << qx0 << "   " << qy0 << "   " << qz0 << "   " << iqx0_loc << "   " << iqy0_loc << "   " << iqz0_loc << "   " << region_index[iqx][iqy][iqz] << "   ";

			switch(region_index[iqx][iqy][iqz])
			{
				case 0:	//X-Y-Z in range
					double point_XYZ[3] = { qx0, qy0, qz0 };
					fleshed_out_CF[iqx][iqy][iqz] = 1.0 + cf.eval(point_XYZ);
					break;
				case 1:	//X-Y in range, Z out of range
					double point_XY[2] = { qx0, qy0 };
					double fZ1 = (*cf_XY[iqz0_loc]).eval(point_XY);
					double fZ2 = (*cf_XY[iqz0_loc+1]).eval(point_XY);
					fleshed_out_CF[iqx][iqy][iqz] = 1.0 + Extrapolate_Gaussian_1D(qz0, qz_pts[iqz0_loc], qz_pts[iqz0_loc+1], fZ1, fZ2);
					break;
				case 2:	//X-Z in range, Y out of range
					double point_XZ[2] = { qx0, qz0 };
					double fY1 = (*cf_XZ[iqy0_loc]).eval(point_XZ);
					double fY2 = (*cf_XZ[iqy0_loc+1]).eval(point_XZ);
					fleshed_out_CF[iqx][iqy][iqz] = 1.0 + Extrapolate_Gaussian_1D(qy0, qy_pts[iqy0_loc], qy_pts[iqy0_loc+1], fY1, fY2);
					break;
				case 3:	//X in range, Y-Z out of range
					double point_X[1] = { qx0 };
					double Gpoint_YZ[2] = { qy0, qz0 };
					double Glower_YZ[2] = { qy_pts[iqy0_loc], qz_pts[iqz0_loc] };
					double Gupper_YZ[2] = { qy_pts[iqy0_loc+1], qz_pts[iqz0_loc+1] };
					double Gvals_YZ[2][2];
					for (int i1 = 0; i1 < 2; ++i1)
					for (int i2 = 0; i2 < 2; ++i2)
						Gvals_YZ[i1][i2] = (*cf_X[iqy0_loc+i1][iqz0_loc+i2]).eval(point_X);
					fleshed_out_CF[iqx][iqy][iqz] = 1.0 + Extrapolate_Gaussian_2D(Gpoint_YZ, Glower_YZ, Gupper_YZ, Gvals_YZ);
					break;
				case 4:	//Y-Z in range, X out of range
					double point_YZ[2] = { qy0, qz0 };
					double fX1 = (*cf_YZ[iqx0_loc]).eval(point_YZ);
					double fX2 = (*cf_YZ[iqx0_loc+1]).eval(point_YZ);
					fleshed_out_CF[iqx][iqy][iqz] = 1.0 + Extrapolate_Gaussian_1D(qx0, qx_pts[iqx0_loc], qx_pts[iqx0_loc+1], fX1, fX2);
					break;
				case 5:	//Y in range, X-Z out of range
					double point_Y[1] = { qy0 };
					double Gpoint_XZ[2] = { qx0, qz0 };
					double Glower_XZ[2] = { qx_pts[iqx0_loc], qz_pts[iqz0_loc] };
					double Gupper_XZ[2] = { qx_pts[iqx0_loc+1], qz_pts[iqz0_loc+1] };
					double Gvals_XZ[2][2];
					for (int i1 = 0; i1 < 2; ++i1)
					for (int i2 = 0; i2 < 2; ++i2)
						Gvals_XZ[i1][i2] = (*cf_Y[iqx0_loc+i1][iqz0_loc+i2]).eval(point_Y);
					fleshed_out_CF[iqx][iqy][iqz] = 1.0 + Extrapolate_Gaussian_2D(Gpoint_XZ, Glower_XZ, Gupper_XZ, Gvals_XZ);
					break;
				case 6:	//Z in range, X-Y out of range
					double point_Z[1] = { qz0 };
					double Gpoint_XY[2] = { qx0, qy0 };
					double Glower_XY[2] = { qx_pts[iqx0_loc], qy_pts[iqy0_loc] };
					double Gupper_XY[2] = { qx_pts[iqx0_loc+1], qy_pts[iqy0_loc+1] };
					double Gvals_XY[2][2];
					for (int i1 = 0; i1 < 2; ++i1)
					for (int i2 = 0; i2 < 2; ++i2)
						Gvals_XY[i1][i2] = (*cf_Z[iqx0_loc+i1][iqy0_loc+i2]).eval(point_Z);
					fleshed_out_CF[iqx][iqy][iqz] = 1.0 + Extrapolate_Gaussian_2D(Gpoint_XY, Glower_XY, Gupper_XY, Gvals_XY);
					break;
				case 7:	//X-Y-Z out of range
					double Gpoint_XYZ[3] = { qx0, qy0, qz0 };
					double Glower_XYZ[3] = { qx_pts[iqx0_loc], qy_pts[iqy0_loc], qz_pts[iqz0_loc] };
					double Gupper_XYZ[3] = { qx_pts[iqx0_loc+1], qy_pts[iqy0_loc+1], qz_pts[iqz0_loc+1] };
					double Gvals_XYZ[2][2][2];
					for (int i1 = 0; i1 < 2; ++i1)
					for (int i2 = 0; i2 < 2; ++i2)
					for (int i3 = 0; i3 < 2; ++i3)
						Gvals_XYZ[i1][i2][i3] = current_C_slice[iqx0_loc+i1][iqy0_loc+i2][iqz0_loc+i3] - 1.0;
					fleshed_out_CF[iqx][iqy][iqz] = 1.0 + Extrapolate_Gaussian_3D(Gpoint_XYZ, Glower_XYZ, Gupper_XYZ, Gvals_XYZ);
					break;
				default:
					cerr << "FitCF::Flesh_out_CF(" << ipt << ", " << ipphi << "): shouldn't have made it to this point!  Exiting..." << endl;
					exit(1);
					break;
			}
cout << fleshed_out_CF[iqx][iqy][iqz] << endl;
//exit(1);
		}

		for (int it = 0; it < cf_XY.size(); ++it)
			delete cf_XY[it];
		for (int it = 0; it < cf_XZ.size(); ++it)
			delete cf_XZ[it];
		for (int it = 0; it < cf_YZ.size(); ++it)
			delete cf_YZ[it];
		cf_XY.clear();
		cf_XZ.clear();
		cf_YZ.clear();

		for (int it = 0; it < cf_X.size(); ++it)
		{
			for (int it2 = 0; it2 < cf_X[it].size(); ++it2)
				delete cf_X[it][it2];
			cf_X[it].clear();
		}
		for (int it = 0; it < cf_Y.size(); ++it)
		{
			for (int it2 = 0; it2 < cf_Y[it].size(); ++it2)
				delete cf_Y[it][it2];
			cf_Y[it].clear();
		}
		for (int it = 0; it < cf_Z.size(); ++it)
		{
			for (int it2 = 0; it2 < cf_Z[it].size(); ++it2)
				delete cf_Z[it][it2];
			cf_Z[it].clear();
		}
		cf_X.clear();
		cf_Y.clear();
		cf_Z.clear();
	}
	else
	{
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
	}

	return;
}

double FitCF::Extrapolate_Gaussian_1D(double q0, double qi0, double qi1, double f1, double f2)
{
	//as a safeguard for the timebeing...
	if (f1 < 0.0 || f2 < 0.0)
		return lin_int(q0 - qi0, 1./(qi1-qi0), f1, f2);
	else	//this is what *should* execute
		return ( exp( lin_int( sgn(q0)*q0*q0 - sgn(qi0)*qi0*qi0, 1./( sgn(qi1)*qi1*qi1 - sgn(qi0)*qi0*qi0 ), log(f1), log(f2) ) ) );
}

double FitCF::Extrapolate_Gaussian_2D(double * q0, double * qi0, double * qi1, double (*vals) [2])
{
	double total = 0.0;
	double slopes[2][2];

	if (vals[0][0] <= 0.0 || vals[0][1] <= 0.0 || vals[1][0] <= 0.0 || vals[1][1] <= 0.0)
	{
		for (int i = 0; i < 2; ++i)
		{
			slopes[i][0] = (qi1[i] - q0[i]) / (qi1[i] - qi0[i]);
			slopes[i][1] = 1.0 - slopes[i][0];
		}
	
		for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
			total += vals[i][j] * slopes[0][i] * slopes[1][j];
	}
	else
	{
		for (int i = 0; i < 2; ++i)
		{
			double sqi1 = sgn(qi1[i])*qi1[i]*qi1[i];
			double sq0 = sgn(q0[i])*q0[i]*q0[i];
			double sqi0 = sgn(qi0[i])*qi0[i]*qi0[i];
			slopes[i][0] = (sqi1 - sq0) / (sqi1 - sqi0);
			slopes[i][1] = 1.0 - slopes[i][0];
		}
	
		for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
			total += log(vals[i][j]) * slopes[0][i] * slopes[1][j];

		total = exp(total);
	}
	
	return (total);
}

double FitCF::Extrapolate_Gaussian_3D(double * q0, double * qi0, double * qi1, double (*vals) [2][2])
{
	double total = 0.0;
	double slopes[3][2];

	if (vals[0][0][0] <= 0.0 || vals[0][0][1] <= 0.0 || vals[0][1][0] <= 0.0 || vals[0][1][1] <= 0.0
			|| vals[1][0][0] <= 0.0 || vals[1][0][1] <= 0.0 || vals[1][1][0] <= 0.0 || vals[1][1][1] <= 0.0)
	{
		for (int i = 0; i < 3; ++i)
		{
			slopes[i][0] = (qi1[i] - q0[i]) / (qi1[i] - qi0[i]);
			slopes[i][1] = 1.0 - slopes[i][0];
		}
	
		for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
		for (int k = 0; k < 2; ++k)
			total += vals[i][j][k] * slopes[0][i] * slopes[1][j] * slopes[2][k];
	}
	else
	{
		for (int i = 0; i < 3; ++i)
		{
			double sqi1 = sgn(qi1[i])*qi1[i]*qi1[i];
			double sq0 = sgn(q0[i])*q0[i]*q0[i];
			double sqi0 = sgn(qi0[i])*qi0[i]*qi0[i];
			slopes[i][0] = (sqi1 - sq0) / (sqi1 - sqi0);
			slopes[i][1] = 1.0 - slopes[i][0];
		}
	
		for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
		for (int k = 0; k < 2; ++k)
			total += log(vals[i][j][k]) * slopes[0][i] * slopes[1][j] * slopes[2][k];

		total = exp(total);
	}
	
	return (total);
}

double FitCF::interpolate_CF(double *** current_C_slice, double qx0, double qy0, double qz0, int ipt, int thermal_or_resonances)
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
		//fxiyizi = interpolate_qi(sqx02, sqx_02, sqx_12, fx0, fx1, false);
	}
	else
	{
		/*int cqx = (qxnpts - 1)/2;           //central qx index
		int cqy = (qynpts - 1)/2;
		int cqz = (qznpts - 1)/2;
		int iqxh0 = (qxnpts - cqx - 1)/2;   //halfway between cqx and nqpts-1
		int iqyh0 = (qynpts - cqy - 1)/2;
		int iqzh0 = (qznpts - cqz - 1)/2;

		//if any of these are false, use logarithmic interpolation in that direction
		//(assuming approximately exponential fall-off)
		bool use_linear_qx = (iqx0_loc < cqx + iqxh0) && (iqx0_loc >= cqx - iqxh0);
		bool use_linear_qy = (iqy0_loc < cqy + iqyh0) && (iqy0_loc >= cqy - iqyh0);
		bool use_linear_qz = (iqz0_loc < cqz + iqzh0) && (iqz0_loc >= cqz - iqzh0);

		//if (thermal_or_resonances == 1)	//if you're doing the crossterm, use linear, no matter what
		//{
			use_linear_qx = true;
			use_linear_qy = true;
			use_linear_qz = true;
		//}

		//compute minimum element of current_C_slice...
		double min = current_C_slice[0][0][0];
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int iqz = 0; iqz < qznpts; ++iqz)
			if (current_C_slice[iqx][iqy][iqz] < min) min = current_C_slice[iqx][iqy][iqz];

		//use minimum value to define offset...
		double offset = 1e-100;
		if (min <= offset)
			offset = abs(min) + 0.001;

		fx0 += offset;
		fx1 += offset;

		//interpolate over qz-points first
		///////////////////////
		//interpolate over first pair of qz-points
		fx0y0z0 += offset;
		fx0y0z1 += offset;
		double fx0y0zi = interpolate_qi(qz0, qz_0, qz_1, fx0y0z0, fx0y0z1, use_linear_qz);
	
		//interpolate over second pair of qz-points
		fx0y1z0 += offset;
		fx0y1z1 += offset;
		double fx0y1zi = interpolate_qi(qz0, qz_0, qz_1, fx0y1z0, fx0y1z1, use_linear_qz);    
	
		//interpolate over third pair of qz-points
		fx1y0z0 += offset;
		fx1y0z1 += offset;
		double fx1y0zi = interpolate_qi(qz0, qz_0, qz_1, fx1y0z0, fx1y0z1, use_linear_qz);
	
		//interpolate over fourth pair of qz-points
		fx1y1z0 += offset;
		fx1y1z1 += offset;
		double fx1y1zi = interpolate_qi(qz0, qz_0, qz_1, fx1y1z0, fx1y1z1, use_linear_qz);
		///////////////////////
	
		//interpolate over qy-points next
		double fx0yizi = interpolate_qi(qy0, qy_0, qy_1, fx0y0zi, fx0y1zi, use_linear_qy);
		double fx1yizi = interpolate_qi(qy0, qy_0, qy_1, fx1y0zi, fx1y1zi, use_linear_qy);
	
		//finally, interpolate over qx-points
		fxiyizi = interpolate_qi(qx0, qx_0, qx_1, fx0yizi, fx1yizi, use_linear_qx) - offset;*/

		fxiyizi = interpolate3D(qx_pts, qy_pts, qz_pts, current_C_slice, qx0, qy0, qz0, qxnpts, qynpts, qznpts, 1, true);

//cout << "offset = " << offset << "   " << qx0 << "   " << qy0 << "   " << qz0 << "   " << ipt << "   " << thermal_or_resonances 
//		<< "   " << fx0y0z0 << "   " << fx0y0z1 << "   " << fx0y1z0 << "   " << fx0y1z1 << "   " << fx1y0z0 << "   " << fx1y0z1 << "   " << fx1y1z0 << "   " << fx1y1z1
//		<< "   " << fx0y0zi << "   " << fx0y1zi << "   " << fx1y0zi << "   " << fx1y1zi << "   " << fx0yizi << "   " << fx1yizi << "   " << fxiyizi 
//		<< "   " << use_linear_qx << "   " << use_linear_qy << "   " << use_linear_qz << endl;
//cout << "Here: " << use_linear_qx << "   " << iqx0_loc << "   " << cqx << "   " << iqxh0 << "   " << thermal_or_resonances << "   ";
//		if (thermal_or_resonances == 1)
//			use_linear_qx = true;
//		fxiyizi = interpolate_qi(qx0, qx_0, qx_1, fx0, fx1, use_linear_qx) - offset;
	}

	return (fxiyizi);
}

double FitCF::interpolate_qi(double q0, double qi0, double qi1, double f1, double f2, bool use_linear)
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

void FitCF::Regulate_CF(int ipt, int iqt, int iqx, int iqy, int iqz, double * CF, double * projCF)
{
	double mu = target_pphiavgd_CFs[ipt][iqt][iqx][iqy][iqz];
	double sigma = sqrt(target_pphivar_CFs[ipt][iqt][iqx][iqy][iqz]);

	double zthreshhold = 3.0;
	double pTcutoff = 2.0;	//GeV

	double zact = abs(mu - *CF) / sigma;

	if (zact >= zthreshhold || abs(*CF-1.5) > 0.500001)
	{
		//if (SPinterp_pT[ipt] < pTcutoff)
			*global_out_stream_ptr << "WARNING: regulated CF point at pT = " << SPinterp_pT[ipt]
				<< ": (" << *CF << "," << *projCF << ") --> (";
		*CF = mu;	//if it's an outlier, replace with pphi-averaged value
		*projCF = mu;
		//if (SPinterp_pT[ipt] < pTcutoff)
			*global_out_stream_ptr << *CF << "," << *projCF << ")" << endl;
	}

	return;
}

void FitCF::Regulate_CF_Hampel(int ipt, int iqx, int iqy, int iqz,
												double * pphi_CF_slice, double * pphi_CF_slice_term1, double * pphi_CF_slice_term2, double * pphi_CF_slice_term3)
{
	bool * is_outlier = new bool [n_interp_pphi_pts];

	double pphi_median = 0.0;

	find_outliers_Hampel(pphi_CF_slice, n_interp_pphi_pts, is_outlier, &pphi_median, 7.5);	//last part is Hampel factor

	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	if (is_outlier[ipphi] || abs(pphi_CF_slice[ipphi]-1.5) > 0.500001)
	{
		//if (SPinterp_pT[ipt] < pTcutoff)
			*global_out_stream_ptr << "\t WARNING: regulated CF point at pT = " << SPinterp_pT[ipt] << ", pphi = " << SPinterp_pphi[ipphi]
				<< " (" << iqx << ", " << iqy << ", " << iqz << "): (" << scientific << setprecision(6) << setw(9)
				<< pphi_CF_slice[ipphi] << ", " << pphi_CF_slice_term1[ipphi] << ", " << pphi_CF_slice_term2[ipphi] << ", " << pphi_CF_slice_term3[ipphi] << ") --> (";
		pphi_CF_slice[ipphi] = pphi_median;	//if it's an outlier, replace with pphi-median value
		pphi_CF_slice_term1[ipphi] = get_median(pphi_CF_slice_term1, n_interp_pphi_pts);
		pphi_CF_slice_term2[ipphi] = get_median(pphi_CF_slice_term2, n_interp_pphi_pts);
		pphi_CF_slice_term3[ipphi] = get_median(pphi_CF_slice_term3, n_interp_pphi_pts);
		//if (SPinterp_pT[ipt] < pTcutoff)
			*global_out_stream_ptr << pphi_CF_slice[ipphi] << ", " << pphi_CF_slice_term1[ipphi] << ", " << pphi_CF_slice_term2[ipphi] << ", " << pphi_CF_slice_term3[ipphi] << ")" << endl;
	}

	delete [] is_outlier;

	return;
}

void FitCF::Regulate_CF_Hampel_v2(int ipt, int iqx, int iqy, int iqz,
												double * pphi_CF_slice, double * pphi_CF_slice_term1, double * pphi_CF_slice_term2, double * pphi_CF_slice_term3)
{
	bool * is_outlier = new bool [n_interp_pphi_pts];
	double * local_pphi_medians = new double [n_interp_pphi_pts];

	int window_width = 11;	//this seems to work fairly well

	find_outliers_window_Hampel(pphi_CF_slice, n_interp_pphi_pts, is_outlier, local_pphi_medians, 5.2, window_width);	//last part is Hampel factor

	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	if (is_outlier[ipphi] || abs(pphi_CF_slice[ipphi]-1.5) > 0.500001)
	{
		//if (SPinterp_pT[ipt] < pTcutoff)
			*global_out_stream_ptr << "\t WARNING: regulated CF point at pT = " << SPinterp_pT[ipt] << ", pphi = " << SPinterp_pphi[ipphi]
				<< " (" << iqx << ", " << iqy << ", " << iqz << "): (" << scientific << setprecision(6) << setw(9)
				<< pphi_CF_slice[ipphi] << ", " << pphi_CF_slice_term1[ipphi] << ", " << pphi_CF_slice_term2[ipphi] << ", " << pphi_CF_slice_term3[ipphi] << ") --> (";
		pphi_CF_slice[ipphi] = local_pphi_medians[ipphi];	//if it's an outlier, replace with pphi-median value
		pphi_CF_slice_term1[ipphi] = get_median(pphi_CF_slice_term1, n_interp_pphi_pts);
		pphi_CF_slice_term2[ipphi] = get_median(pphi_CF_slice_term2, n_interp_pphi_pts);
		pphi_CF_slice_term3[ipphi] = get_median(pphi_CF_slice_term3, n_interp_pphi_pts);
		//if (SPinterp_pT[ipt] < pTcutoff)
			*global_out_stream_ptr << pphi_CF_slice[ipphi] << ", " << pphi_CF_slice_term1[ipphi] << ", "
				<< pphi_CF_slice_term2[ipphi] << ", " << pphi_CF_slice_term3[ipphi] << ")" << endl;
	}


	delete [] is_outlier;

	return;
}

double FitCF::get_CF(int ipt, int ipphi, int iqt, int iqx, int iqy, int iqz, bool return_projected_value)
{
	double nonFTd_spectra = spectra[target_particle_id][ipt][ipphi];
	//double cos_transf_spectra = full_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0];
	//double sin_transf_spectra = full_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1];
	double cos_transf_spectra = full_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 0)];
	double sin_transf_spectra = full_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 1)];

	if (return_projected_value)
	{
		//with no resonances
		double nonFTd_tspectra = thermal_spectra[target_particle_id][ipt][ipphi];
		//double cos_transf_tspectra = thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0];
		//double sin_transf_tspectra = thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1];
		double cos_transf_tspectra = thermal_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 0)];
		double sin_transf_tspectra = thermal_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 1)];

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

void FitCF::get_CF(double * totalresult, double * thermalresult, double * crosstermresult, double * resonanceresult,
									int ipt, int ipphi, int iqt, int iqx, int iqy, int iqz, bool return_projected_value)
{
	//thermal
	double nonFTd_tspectra = thermal_spectra[target_particle_id][ipt][ipphi];
	//double cos_transf_tspectra = thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0];
	//double sin_transf_tspectra = thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1];
	double cos_transf_tspectra = thermal_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 0)];
	double sin_transf_tspectra = thermal_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 1)];
	//total
	double nonFTd_spectra = spectra[target_particle_id][ipt][ipphi];
	//double cos_transf_spectra = full_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0];
	//double sin_transf_spectra = full_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1];
	double cos_transf_spectra = full_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 0)];
	double sin_transf_spectra = full_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 1)];

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
	//den = nonFTd_tspectra*nonFTd_tspectra;
	*crosstermresult = num / den;

	num = cosNT_spectra*cosNT_spectra+sinNT_spectra*sinNT_spectra;								//resonances only
	//den = nonFTd_tspectra*nonFTd_tspectra;
	*resonanceresult = num / den;

	num = cos_transf_spectra*cos_transf_spectra + sin_transf_spectra*sin_transf_spectra;
	//den = nonFTd_spectra*nonFTd_spectra;
	*totalresult = num / den;

	//cout << "CHECK (" << ipt << ", " << ipphi << ", " << iqt << ", " << iqx << ", " << iqy << ", " << iqz << "): "
	//		<< fraction_of_resonances << "   " << nonFTd_tspectra << "   " << cos_transf_tspectra << "   " << sin_transf_tspectra << "   "
	//		<< nonFTd_spectra << "   " << cos_transf_spectra << "   " << sin_transf_spectra << "   "
	//		<< *thermalresult << "   " << *crosstermresult << "   " << *resonanceresult << "   " << *totalresult << endl;

	//this is just for debugging: get rid of it if you don't remember what it does!
	/*if (return_projected_value)
	{
		nonFTd_spectra = spectra[target_particle_id][ipt][ipphi];
		//cos_transf_spectra = full_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0];
		//sin_transf_spectra = full_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1];
		cos_transf_spectra = full_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 0)];
		sin_transf_spectra = full_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 1)];
		double CF0 = (cos_transf_tspectra*cos_transf_tspectra + sin_transf_tspectra*sin_transf_tspectra) / (nonFTd_tspectra*nonFTd_tspectra);
		double CF1 = (cos_transf_spectra*cos_transf_spectra + sin_transf_spectra*sin_transf_spectra) / (nonFTd_spectra*nonFTd_spectra);
		*totalresult = CF0 + (CF1 - CF0) / fraction_of_resonances;
	}*/


	return;
}

//**************************************************************
// Gaussian fit routines below
//**************************************************************

void FitCF::Fit_Correlationfunction3D(double *** Correl_3D, int ipt, int ipphi, bool fleshing_out_CF /*== true*/)
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
	const size_t n_para = 4;  // # of parameters

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
	double ckp = cos_SPinterp_pphi[ipphi], skp = sin_SPinterp_pphi[ipphi];
	for (int i = 0; i < q1npts; i++)
	for (int j = 0; j < q2npts; j++)
	for (int k = 0; k < q3npts; k++)
	{
		Correlfun3D_data.q_o[idx] = q1pts[i] * ckp + q2pts[j] * skp;
		Correlfun3D_data.q_s[idx] = -q1pts[i] * skp + q2pts[j] * ckp;
		Correlfun3D_data.q_l[idx] = q3pts[k];
		Correlfun3D_data.y[idx] = Correl_3D[i][j][k];
		//Correlfun3D_data.sigma[idx] = Correl_3D_err[i][j][k];
		Correlfun3D_data.sigma[idx] = 1.e-2;
if (i==(q1npts-1)/2 && j==(q2npts-1)/2 && k==(q3npts-1)/2)
	Correlfun3D_data.sigma[idx] = 1.e10;	//ignore central point
		idx++;
	}

	double para_init[n_para] = { 1.0, 1.0, 1.0, 1.0 };  // initial guesses of parameters

	gsl_vector_view xvec_ptr = gsl_vector_view_array (para_init, n_para);

	// set up the function to be fit 
	gsl_multifit_function_fdf target_func;
	target_func.f = &Fittarget_correlfun3D_f;        // the function of residuals
	target_func.df = &Fittarget_correlfun3D_df;      // the gradient of this function
	target_func.fdf = &Fittarget_correlfun3D_fdf;    // combined function and gradient
	target_func.n = data_length;              // number of points in the data set
	target_func.p = n_para;              // number of parameters in the fit function
	target_func.params = &Correlfun3D_data;  // structure with the data and error bars

	const gsl_multifit_fdfsolver_type *type_ptr = gsl_multifit_fdfsolver_lmsder;
	gsl_multifit_fdfsolver *solver_ptr = gsl_multifit_fdfsolver_alloc (type_ptr, data_length, n_para);
	gsl_multifit_fdfsolver_set (solver_ptr, &target_func, &xvec_ptr.vector);

	size_t iteration = 0;         // initialize iteration counter
	print_fit_state_3D (iteration, solver_ptr);
	int status;  		// return value from gsl function calls (e.g., error)
	do
	{
		iteration++;
      
		// perform a single iteration of the fitting routine
		status = gsl_multifit_fdfsolver_iterate (solver_ptr);

		// print out the status of the fit
		if (VERBOSE > 2) cout << "status = " << gsl_strerror (status) << endl;

		// customized routine to print out current parameters
		if (VERBOSE > 2) print_fit_state_3D (iteration, solver_ptr);

		if (status)    // check for a nonzero status code
		{
			break;  // this should only happen if an error code is returned 
		}

		// test for convergence with an absolute and relative error (see manual)
		status = gsl_multifit_test_delta (solver_ptr->dx, solver_ptr->x, fit_tolerance, fit_tolerance);
	}
	while (status == GSL_CONTINUE && iteration < fit_max_iterations);

	//cerr >> "iteration = " << iteration << endl;

	// calculate the covariance matrix of the best-fit parameters
	gsl_multifit_covar (solver_ptr->J, 0.0, covariance_ptr);

	// print out the covariance matrix using the gsl function (not elegant!)
	if (VERBOSE > 2) cout << endl << "Covariance matrix: " << endl;
	if (VERBOSE > 2) gsl_matrix_fprintf (stdout, covariance_ptr, "%g");

	cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
	cout.precision (5);		                // # of digits in doubles

	int width = 7;		// setw width for output
	/*cout << endl << "Best fit results:" << endl;
	cout << "R2o = " << setw (width) << get_fit_results (0, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (0, covariance_ptr) << endl;

	cout << "R2s      = " << setw (width) << get_fit_results (1, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (1, covariance_ptr) << endl;

	cout << "R2l      = " << setw (width) << get_fit_results (2, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (2, covariance_ptr) << endl;
  
	cout << "R2os      = " << setw (width) << get_fit_results (3, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (3, covariance_ptr) << endl;
    
	cout << "status = " << gsl_strerror (status) << endl;
	cout << "--------------------------------------------------------------------" << endl;*/

	double chi = gsl_blas_dnrm2(solver_ptr->f);
	double dof = data_length - n_para;
	double c = GSL_MAX_DBL(1, chi/sqrt(dof));

	lambda_Correl[ipt][ipphi] = 1.0;
	lambda_Correl_err[ipt][ipphi] = 0.0;
	R2_out_GF[ipt][ipphi] = fabs(get_fit_results(0, solver_ptr))*hbarC*hbarC;
	R2_side_GF[ipt][ipphi] = fabs(get_fit_results(1, solver_ptr))*hbarC*hbarC;
	R2_long_GF[ipt][ipphi] = fabs(get_fit_results(2, solver_ptr))*hbarC*hbarC;
	R2_outside_GF[ipt][ipphi] = get_fit_results(3, solver_ptr)*hbarC*hbarC;
	R2_out_err[ipt][ipphi] = c*get_fit_err(0, covariance_ptr)*hbarC*hbarC;
	R2_side_err[ipt][ipphi] = c*get_fit_err(1, covariance_ptr)*hbarC*hbarC;
	R2_long_err[ipt][ipphi] = c*get_fit_err(2, covariance_ptr)*hbarC*hbarC;
	R2_outside_err[ipt][ipphi] = c*get_fit_err(3, covariance_ptr)*hbarC*hbarC;

	if (VERBOSE > 1) 
	{
		cout << "final results: " << endl;
		cout << scientific << setw(10) << setprecision(5) 
			<< "chisq/dof = " << chi*chi/dof << endl;
		cout << scientific << setw(10) << setprecision(5) 
			<< " lambda[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << lambda_Correl[ipt][ipphi] << " +/- " << lambda_Correl_err[ipt][ipphi] << endl;
		cout << " R2_out[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << R2_out_GF[ipt][ipphi] << " +/- " << R2_out_err[ipt][ipphi] << endl;
		cout << " R2_side[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << R2_side_GF[ipt][ipphi] << " +/- " << R2_side_err[ipt][ipphi] << endl;
		cout << " R2_long[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << R2_long_GF[ipt][ipphi] << " +/- " << R2_long_err[ipt][ipphi] << endl;
		cout << " R2_outside[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << R2_outside_GF[ipt][ipphi] << " +/- " << R2_outside_err[ipt][ipphi] << endl;
	}

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

void FitCF::Fit_Correlationfunction3D_withlambda(double *** Correl_3D, int ipt, int ipphi, bool fleshing_out_CF /*== true*/)
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
	double ckp = cos_SPinterp_pphi[ipphi], skp = sin_SPinterp_pphi[ipphi];
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
if (i==(q1npts-1)/2 && j==(q2npts-1)/2 && k==(q3npts-1)/2)
	Correlfun3D_data.sigma[idx] = 1.e10;	//ignore central point
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
	print_fit_state_3D (iteration, solver_ptr);
	int status;  		// return value from gsl function calls (e.g., error)
	do
	{
		iteration++;
      
		// perform a single iteration of the fitting routine
		status = gsl_multifit_fdfsolver_iterate (solver_ptr);

		// print out the status of the fit
		if (VERBOSE > 2) cout << "status = " << gsl_strerror (status) << endl;

		// customized routine to print out current parameters
		if (VERBOSE > 2) print_fit_state_3D (iteration, solver_ptr);

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
	if (VERBOSE > 2) cout << endl << "Covariance matrix: " << endl;
	if (VERBOSE > 2) gsl_matrix_fprintf (stdout, covariance_ptr, "%g");

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

	if (VERBOSE > 1) 
	{
		cout << "final results: " << endl;
		cout << scientific << setw(10) << setprecision(5) 
			<< "chisq/dof = " << chi*chi/dof << endl;
		cout << scientific << setw(10) << setprecision(5) 
			<< " lambda = " << lambda_Correl[ipt][ipphi] << " +/- " << lambda_Correl_err[ipt][ipphi] << endl;
		cout << " R2_out[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << R2_out_GF[ipt][ipphi] << " +/- " << R2_out_err[ipt][ipphi] << endl;
		cout << " R2_side[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << R2_side_GF[ipt][ipphi] << " +/- " << R2_side_err[ipt][ipphi] << endl;
		cout << " R2_long[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << R2_long_GF[ipt][ipphi] << " +/- " << R2_long_err[ipt][ipphi] << endl;
		cout << " R2_outside[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << R2_outside_GF[ipt][ipphi] << " +/- " << R2_outside_err[ipt][ipphi] << endl;
	}

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

//*********************************************************************
// 3D case
//*********************************************************************
//  Simple function to print results of each iteration in nice format
int FitCF::print_fit_state_3D (size_t iteration, gsl_multifit_fdfsolver * solver_ptr)
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
int FitCF::print_fit_state_3D_withlambda (size_t iteration, gsl_multifit_fdfsolver * solver_ptr)
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
