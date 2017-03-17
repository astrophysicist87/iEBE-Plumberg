#ifndef MAIN_H
#define MAIN_H

#include <string>
#include <sstream>

#include "Arsenal.h"
#include "gauss_quadrature.h"

using namespace std;

#define USE_SCALING 	0	// quick bug fix - delete once events have been re-run
#define DO_SLICES_ONLY	1	// do OSL slices only vs. every bin of 3D correlation function

string workingDirectory = "";

const int nKT = 4;
const int n_bin_centers = 13;
const int n_qpts_per_bin = 15;
const int nqpts = n_qpts_per_bin*n_qpts_per_bin*n_qpts_per_bin;
const int n_pT_pts = 15;
const int n_pphi_pts = 48;
const int nqx = 11;
const int nqy = 11;
const int nqz = 11;

const double Kphi_min = 0.0, Kphi_max = 2.0*M_PI;

double SP_pT[n_pT_pts], SP_pT_wts[n_pT_pts];
double SP_pphi[n_pphi_pts], SP_pphi_wts[n_pphi_pts];
double qx_pts[nqx];
double qy_pts[nqy];
double qz_pts[nqz];

double spectra[n_pT_pts][n_pphi_pts];
double thermalCFvals[n_pT_pts][n_pphi_pts][nqx][nqy][nqz];
double crosstermCFvals[n_pT_pts][n_pphi_pts][nqx][nqy][nqz];
double resonancesCFvals[n_pT_pts][n_pphi_pts][nqx][nqy][nqz];
double correlation_function[n_pT_pts][n_pphi_pts][nqx][nqy][nqz];

double *** CF_KTKphi_grid, *** CF_KTKphi_grid_T;
double CF_KT_grid[n_pT_pts][nqpts], CF_KT_grid_T[nqpts][n_pT_pts];

inline double lin_int(double x_m_x1, double one_by_x2_m_x1, double f1, double f2)
{
	return ( f1 + (f2 - f1) * x_m_x1 * one_by_x2_m_x1 );
}

/*inline int indexer(const int ipt, const int ipphi, const int iq)
{
	return (
		( ipt * n_pphi_pts + ipphi ) * nqpts + iq
	);
}*/

void Read_in_correlationfunction()
{
	ifstream iCorrFunc;
	iCorrFunc.open( (workingDirectory + "/correlfunct3D_Pion_+.dat").c_str() );

	double scaling_factors[n_pT_pts][n_pphi_pts];
	ifstream iCorrFunc2;
	if (USE_SCALING)
		iCorrFunc2.open( (workingDirectory + "/results-1/Qeq0_correlfunct3D_Pion_+.dat").c_str() );

	double dummy;
	for (int ipt = 0; ipt < n_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	{
		if (USE_SCALING)
		{
			iCorrFunc2 >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> scaling_factors[ipt][ipphi];
			scaling_factors[ipt][ipphi] -= 1.0;
		}
		else
			scaling_factors[ipt][ipphi] = 1.0;
		for (int iqx = 0; iqx < nqx; ++iqx)
		for (int iqy = 0; iqy < nqy; ++iqy)
		for (int iqz = 0; iqz < nqz; ++iqz)
		{
			iCorrFunc
				>> SP_pT[ipt]
				>> SP_pphi[ipphi]
				>> qx_pts[iqx]
				>> qy_pts[iqy]
				>> qz_pts[iqz]
				>> spectra[ipt][ipphi]
				>> thermalCFvals[ipt][ipphi][iqx][iqy][iqz]
				>> crosstermCFvals[ipt][ipphi][iqx][iqy][iqz]
				>> resonancesCFvals[ipt][ipphi][iqx][iqy][iqz]
				>> correlation_function[ipt][ipphi][iqx][iqy][iqz];
			thermalCFvals[ipt][ipphi][iqx][iqy][iqz] *= scaling_factors[ipt][ipphi];
			crosstermCFvals[ipt][ipphi][iqx][iqy][iqz] *= scaling_factors[ipt][ipphi];
			resonancesCFvals[ipt][ipphi][iqx][iqy][iqz] *= scaling_factors[ipt][ipphi];
			correlation_function[ipt][ipphi][iqx][iqy][iqz] = 1.0 + scaling_factors[ipt][ipphi] *
																		(correlation_function[ipt][ipphi][iqx][iqy][iqz] - 1.0);
		}
	}

	iCorrFunc.close();
	if (USE_SCALING)
		iCorrFunc2.close();
				
	return;
}

double interpolate_qi(double q0, double qi0, double qi1, double f1, double f2, bool use_linear)
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

double interpolate_CF(double current_C_slice[][nqy][nqz], double qx0, double qy0, double qz0, int ipt, int ipphi, int thermal_or_resonances)
{
	//int thermal_or_resonances - 	0: interpolate thermal CF (assuming basically Gaussian)
	//								1: interpolate resonance contributions (linear near origin, exponential further out)

	int iqx0_loc = binarySearch(qx_pts, nqx, qx0, true, true);
	int iqy0_loc = binarySearch(qy_pts, nqy, qy0, true, true);
	int iqz0_loc = binarySearch(qz_pts, nqz, qz0, true, true);

	if (iqx0_loc == -1 || iqy0_loc == -1 || iqz0_loc == -1)
	{
		cerr << "Interpolation failed: exiting!" << endl;
		exit(1);
	}
	
	double *** alt_current_C_slice = new double ** [nqx];
	for (int iqx = 0; iqx < nqx; ++iqx)
	{
		alt_current_C_slice[iqx] = new double * [nqy];
		for (int iqy = 0; iqy < nqy; ++iqy)
		{
			alt_current_C_slice[iqx][iqy] = new double [nqz];
			for (int iqz = 0; iqz < nqz; ++iqz)
				alt_current_C_slice[iqx][iqy][iqz] = current_C_slice[iqx][iqy][iqz];
		}
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
		fxiyizi = interpolate3D(qx_pts, qy_pts, qz_pts, alt_current_C_slice, qx0, qy0, qz0, nqx, nqy, nqz, 1, true);
	}

	for (int iqx = 0; iqx < nqx; ++iqx)
	{
		for (int iqy = 0; iqy < nqy; ++iqy)
			delete [] alt_current_C_slice[iqx][iqy];
		delete [] alt_current_C_slice[iqx];
	}
	delete [] alt_current_C_slice;

	return (fxiyizi);
}

void Evaluate_CF(int ipt, int ipphi, double * q1pts, double * q2pts, double * q3pts, double * CFvals, int nqpts, int coords)
{
	// coords:
	// 		0 - (q1, q2, q3) = (qx, qy, qz)
	// 		1 - (q1, q2, q3) = (qo, qs, ql)
	
	double a = 1.0, b = 0.0, c = 1.0;
	if (coords == 1)
	{
		a = cos(SP_pphi[ipphi]);
		b = sin(SP_pphi[ipphi]);
		c = 1.0;
	}
	
	for (int iq = 0; iq < nqpts; ++iq)
	{
		double qx0 = a * q1pts[iq] - b * q2pts[iq];
		double qy0 = b * q1pts[iq] + a * q2pts[iq];
		double qz0 = q3pts[iq];

		double tmp_thermal = interpolate_CF(thermalCFvals[ipt][ipphi], qx0, qy0, qz0, ipt, ipphi, 0);
		double tmp_crossterm = interpolate_CF(crosstermCFvals[ipt][ipphi], qx0, qy0, qz0, ipt, ipphi, 1);
		double tmp_resonances = interpolate_CF(resonancesCFvals[ipt][ipphi], qx0, qy0, qz0, ipt, ipphi, 2);

		CFvals[iq] = 1.0 + tmp_thermal + tmp_crossterm + tmp_resonances;
	}

	return;
}

void Allocate_CF_grids(int nqpts)
{
	//KT-Kphi (2D) grids
	CF_KTKphi_grid = new double ** [n_pT_pts];
	for (int ipt = 0; ipt < n_pT_pts; ++ipt)
	{
		CF_KTKphi_grid[ipt] = new double * [n_pphi_pts];
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			CF_KTKphi_grid[ipt][ipphi] = new double [nqpts];
			for (int iq = 0; iq < nqpts; ++iq)
				CF_KTKphi_grid[ipt][ipphi][iq] = 0.0;
		}
	}
	CF_KTKphi_grid_T = new double ** [nqpts];
	for (int iq = 0; iq < nqpts; ++iq)
	{
		CF_KTKphi_grid_T[iq] = new double * [n_pT_pts];
		for (int ipt = 0; ipt < n_pT_pts; ++ipt)
		{
			CF_KTKphi_grid_T[iq][ipt] = new double [n_pphi_pts];
			for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
				CF_KTKphi_grid_T[iq][ipt][ipphi] = CF_KTKphi_grid[ipt][ipphi][iq];	//transposed grid is useful for interpolation
		}
	}
	return;
}

void Reset_CF_grids(int nqpts)
{
	//KT-Kphi (2D) grids
	for (int ipt = 0; ipt < n_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	for (int iq = 0; iq < nqpts; ++iq)
	{
		CF_KTKphi_grid[ipt][ipphi][iq] = 0.0;
		CF_KTKphi_grid_T[iq][ipt][ipphi] = 0.0;
	}

	return;
}

void Set_CF_grids(double * q1pts, double * q2pts, double * q3pts, int nqpts, int coords)
{
	//KT-Kphi (2D) grids
	for (int ipt = 0; ipt < n_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		Evaluate_CF(ipt, ipphi, q1pts, q2pts, q3pts, CF_KTKphi_grid[ipt][ipphi], nqpts, coords);

	for (int ipt = 0; ipt < n_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	for (int iq = 0; iq < nqpts; ++iq)
		CF_KTKphi_grid_T[iq][ipt][ipphi] = CF_KTKphi_grid[ipt][ipphi][iq];	//transposed grid is useful for interpolation

	return;
}

void Evaluate_CF(double KT, double Kphi, double * CFvals, int nqpts)
{
	for (int iq = 0; iq < nqpts; ++iq)
		CFvals[iq] = interpolate2D(SP_pT, SP_pphi, CF_KTKphi_grid_T[iq], KT, Kphi, n_pT_pts, n_pphi_pts, 0, false);

	return;
}

void Set_Kphiavgd_CF_grids()
{
	for (int ipt = 0; ipt < n_pT_pts; ++ipt)
	for (int iq = 0; iq < nqpts; ++iq)
	{
		CF_KT_grid[ipt][iq] = 0.0;
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
			CF_KT_grid[ipt][iq] += SP_pphi_wts[ipphi] * CF_KTKphi_grid[ipt][ipphi][iq] / (2.0*M_PI);
		CF_KT_grid_T[iq][ipt] = CF_KT_grid[ipt][iq];								//transposed grid is useful for interpolation
	}

	return;
}

void Evaluate_Kphiavgd_CF(double KT, double * CFvals, int nqpts)
{
	for (int iq = 0; iq < nqpts; ++iq)
		CFvals[iq] = interpolate1D(SP_pT, CF_KT_grid_T[iq], KT, n_pT_pts, 0, false);

	return;
}

#endif
