#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include<math.h>
#include<sys/time.h>
#include<algorithm>

#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#include "main.h"

using namespace std;

int main(int argc, char *argv[])
{
    //workingDirectory = "../RESULTS_Edec300/results";
	workingDirectory = argv[1];

	////////////////////////////////////////////
	// Actual calculations start here...
	////////////////////////////////////////////

	//do calculations
	const int coordinates = 1;	//1: input is (qo, qs, ql)
	double bin_centers[n_bin_centers];	//same in each of the three directions
	const double delta_q = 0.01;	//bin width of 10 MeV
	linspace(bin_centers, -0.05, 0.05, n_bin_centers);

	Read_in_correlationfunction();

	gauss_quadrature(n_pT_pts, 5, 0.0, 0.0, 0.0, 13.0, SP_pT, SP_pT_wts);
	gauss_quadrature(n_pphi_pts, 1, 0.0, 0.0, Kphi_min, Kphi_max, SP_pphi, SP_pphi_wts);

	gauss_quadrature(n_KTpts_per_bin, 1, 0.0, 0.0, -1.0, 1.0, base_xpts, base_xwts);
	double KTlls[nKTbins], KTuls[nKTbins];
	linspace(KTlls, 0.0, 0.8, nKTbins);
	linspace(KTuls, 0.2, 1.0, nKTbins);
	KTlls[0] += SP_pT[0]+1.e-6;
	
	double * qo_slice_pts = new double [n_bin_centers*n_qpts_per_bin];
	double * qs_slice_pts = new double [n_bin_centers*n_qpts_per_bin];
	double * ql_slice_pts = new double [n_bin_centers*n_qpts_per_bin];
	double * qo_slice_wts = new double [n_bin_centers*n_qpts_per_bin];
	double * qs_slice_wts = new double [n_bin_centers*n_qpts_per_bin];
	double * ql_slice_wts = new double [n_bin_centers*n_qpts_per_bin];
	double * dummy_pts = new double [n_qpts_per_bin];
	double * dummy_wts = new double [n_qpts_per_bin];

	for (int ibo = 0; ibo < n_bin_centers; ++ibo)
	{
		gauss_quadrature(n_qpts_per_bin, 1, 0.0, 0.0, bin_centers[ibo] - 0.5*delta_q, bin_centers[ibo] + 0.5*delta_q, dummy_pts, dummy_wts);
		for (int iqo = 0; iqo < n_qpts_per_bin; ++iqo)
		{
			qo_slice_pts[ibo*n_qpts_per_bin + iqo] = dummy_pts[iqo];
			qo_slice_wts[ibo*n_qpts_per_bin + iqo] = dummy_wts[iqo];
		}
	}

	for (int ibs = 0; ibs < n_bin_centers; ++ibs)
	{
		gauss_quadrature(n_qpts_per_bin, 1, 0.0, 0.0, bin_centers[ibs] - 0.5*delta_q, bin_centers[ibs] + 0.5*delta_q, dummy_pts, dummy_wts);
		for (int iqs = 0; iqs < n_qpts_per_bin; ++iqs)
		{
			qs_slice_pts[ibs*n_qpts_per_bin + iqs] = dummy_pts[iqs];
			qs_slice_wts[ibs*n_qpts_per_bin + iqs] = dummy_wts[iqs];
		}
	}

	for (int ibl = 0; ibl < n_bin_centers; ++ibl)
	{
		gauss_quadrature(n_qpts_per_bin, 1, 0.0, 0.0, bin_centers[ibl] - 0.5*delta_q, bin_centers[ibl] + 0.5*delta_q, dummy_pts, dummy_wts);
		for (int iql = 0; iql < n_qpts_per_bin; ++iql)
		{
			ql_slice_pts[ibl*n_qpts_per_bin + iql] = dummy_pts[iql];
			ql_slice_wts[ibl*n_qpts_per_bin + iql] = dummy_wts[iql];
		}
	}

	Allocate_CF_grids(nqpts);
	Set_Kphiavgd_spectra();

	double * CFresults = new double [nqpts];

	//////////////
	// Do qo-slice
	//////////////
	for (int ibo = 0; ibo < n_bin_centers; ++ibo)
	{
		int ibs = (n_bin_centers - 1) / 2;
		int ibl = (n_bin_centers - 1) / 2;

		double * qopts = new double [nqpts];
		double * qspts = new double [nqpts];
		double * qlpts = new double [nqpts];
		double * qowts = new double [nqpts];
		double * qswts = new double [nqpts];
		double * qlwts = new double [nqpts];

		Reset_CF_grids(nqpts);
		
		int iq = 0;
		for (int iqo = 0; iqo < n_qpts_per_bin; ++iqo)
		for (int iqs = 0; iqs < n_qpts_per_bin; ++iqs)
		for (int iql = 0; iql < n_qpts_per_bin; ++iql)
		{
			qopts[iq] = qo_slice_pts[ibo*n_qpts_per_bin + iqo];
			qowts[iq] = qo_slice_wts[ibo*n_qpts_per_bin + iqo];
			qspts[iq] = qs_slice_pts[ibs*n_qpts_per_bin + iqs];
			qswts[iq] = qs_slice_wts[ibs*n_qpts_per_bin + iqs];
			qlpts[iq] = ql_slice_pts[ibl*n_qpts_per_bin + iql];
			qlwts[iq] = ql_slice_wts[ibl*n_qpts_per_bin + iql];
			++iq;
		}

		//cout << "  --> Setting grids" << endl;
		//Try calculating any q(o,s,l), and any KT, Kphi (or just any KT, for the Kphi-averaged version)
		//set up full grids in KT-Kphi
		Set_CF_grids(qopts, qspts, qlpts, nqpts, coordinates);

		//cout << "  --> Doing K_phi-averaging" << endl;
		//set up Kphi-averaged grids as well
		Set_Kphiavgd_CF_grids();

		for (int iKT = 0; iKT < nKTbins; ++iKT)
		{
			//reset array to hold interpolated results
			for (int iq = 0; iq < nqpts; ++iq)
				CFresults[iq] = 0.0;

			//cout << "  --> Evaluating at K_T = " << KTpts[iKT] << endl;
			//evaluate at appropriate points
			Get_KTbinaveraged_Kphiavgd_CF(KTlls[iKT], KTuls[iKT], CFresults, nqpts);

			//compute average by taking weighted sum over 
			double CF_averaged_over_qbin = 0.0;
			for (int iq = 0; iq < nqpts; ++iq)
				CF_averaged_over_qbin += qowts[iq]*qswts[iq]*qlwts[iq]*(CFresults[iq] + 1.0);

			CF_averaged_over_qbin /= delta_q*delta_q*delta_q;

			cout << KTlls[iKT] << "   " << KTuls[iKT] << "   " << bin_centers[ibo] << "   " << bin_centers[ibs] << "   " << bin_centers[ibl] << "   " << CF_averaged_over_qbin << endl;
		}

		delete [] qopts;
		delete [] qowts;
		delete [] qspts;
		delete [] qswts;
		delete [] qlpts;
		delete [] qlwts;
	}

	//////////////
	// Do qs-slice
	//////////////
	for (int ibs = 0; ibs < n_bin_centers; ++ibs)
	{
		int ibo = (n_bin_centers - 1) / 2;
		int ibl = (n_bin_centers - 1) / 2;

		double * qopts = new double [nqpts];
		double * qspts = new double [nqpts];
		double * qlpts = new double [nqpts];
		double * qowts = new double [nqpts];
		double * qswts = new double [nqpts];
		double * qlwts = new double [nqpts];

		Reset_CF_grids(nqpts);
		
		int iq = 0;
		for (int iqo = 0; iqo < n_qpts_per_bin; ++iqo)
		for (int iqs = 0; iqs < n_qpts_per_bin; ++iqs)
		for (int iql = 0; iql < n_qpts_per_bin; ++iql)
		{
			qopts[iq] = qo_slice_pts[ibo*n_qpts_per_bin + iqo];
			qowts[iq] = qo_slice_wts[ibo*n_qpts_per_bin + iqo];
			qspts[iq] = qs_slice_pts[ibs*n_qpts_per_bin + iqs];
			qswts[iq] = qs_slice_wts[ibs*n_qpts_per_bin + iqs];
			qlpts[iq] = ql_slice_pts[ibl*n_qpts_per_bin + iql];
			qlwts[iq] = ql_slice_wts[ibl*n_qpts_per_bin + iql];
			++iq;
		}

		//cout << "  --> Setting grids" << endl;
		//Try calculating any q(o,s,l), and any KT, Kphi (or just any KT, for the Kphi-averaged version)
		//set up full grids in KT-Kphi
		Set_CF_grids(qopts, qspts, qlpts, nqpts, coordinates);

		//cout << "  --> Doing K_phi-averaging" << endl;
		//set up Kphi-averaged grids as well
		Set_Kphiavgd_CF_grids();

		for (int iKT = 0; iKT < nKTbins; ++iKT)
		{
			//reset array to hold interpolated results
			for (int iq = 0; iq < nqpts; ++iq)
				CFresults[iq] = 0.0;

			//cout << "  --> Evaluating at K_T = " << KTpts[iKT] << endl;
			//evaluate at appropriate points
			Get_KTbinaveraged_Kphiavgd_CF(KTlls[iKT], KTuls[iKT], CFresults, nqpts);

			//compute average by taking weighted sum over 
			double CF_averaged_over_qbin = 0.0;
			for (int iq = 0; iq < nqpts; ++iq)
				CF_averaged_over_qbin += qowts[iq]*qswts[iq]*qlwts[iq]*(CFresults[iq] + 1.0);

			CF_averaged_over_qbin /= delta_q*delta_q*delta_q;

			cout << KTlls[iKT] << "   " << KTuls[iKT] << "   " << bin_centers[ibo] << "   " << bin_centers[ibs] << "   " << bin_centers[ibl] << "   " << CF_averaged_over_qbin << endl;
		}

		delete [] qopts;
		delete [] qowts;
		delete [] qspts;
		delete [] qswts;
		delete [] qlpts;
		delete [] qlwts;
	}

	//////////////
	// Do ql-slice
	//////////////
	for (int ibl = 0; ibl < n_bin_centers; ++ibl)
	{
		int ibo = (n_bin_centers - 1) / 2;
		int ibs = (n_bin_centers - 1) / 2;

		double * qopts = new double [nqpts];
		double * qspts = new double [nqpts];
		double * qlpts = new double [nqpts];
		double * qowts = new double [nqpts];
		double * qswts = new double [nqpts];
		double * qlwts = new double [nqpts];

		Reset_CF_grids(nqpts);
		
		int iq = 0;
		for (int iqo = 0; iqo < n_qpts_per_bin; ++iqo)
		for (int iqs = 0; iqs < n_qpts_per_bin; ++iqs)
		for (int iql = 0; iql < n_qpts_per_bin; ++iql)
		{
			qopts[iq] = qo_slice_pts[ibo*n_qpts_per_bin + iqo];
			qowts[iq] = qo_slice_wts[ibo*n_qpts_per_bin + iqo];
			qspts[iq] = qs_slice_pts[ibs*n_qpts_per_bin + iqs];
			qswts[iq] = qs_slice_wts[ibs*n_qpts_per_bin + iqs];
			qlpts[iq] = ql_slice_pts[ibl*n_qpts_per_bin + iql];
			qlwts[iq] = ql_slice_wts[ibl*n_qpts_per_bin + iql];
			++iq;
		}

		//cout << "  --> Setting grids" << endl;
		//Try calculating any q(o,s,l), and any KT, Kphi (or just any KT, for the Kphi-averaged version)
		//set up full grids in KT-Kphi
		Set_CF_grids(qopts, qspts, qlpts, nqpts, coordinates);

		//cout << "  --> Doing K_phi-averaging" << endl;
		//set up Kphi-averaged grids as well
		Set_Kphiavgd_CF_grids();

		for (int iKT = 0; iKT < nKTbins; ++iKT)
		{
			//reset array to hold interpolated results
			for (int iq = 0; iq < nqpts; ++iq)
				CFresults[iq] = 0.0;

			//cout << "  --> Evaluating at K_T = " << KTpts[iKT] << endl;
			//evaluate at appropriate points
			Get_KTbinaveraged_Kphiavgd_CF(KTlls[iKT], KTuls[iKT], CFresults, nqpts);

			//compute average by taking weighted sum over 
			double CF_averaged_over_qbin = 0.0;
			for (int iq = 0; iq < nqpts; ++iq)
				CF_averaged_over_qbin += qowts[iq]*qswts[iq]*qlwts[iq]*(CFresults[iq] + 1.0);

			CF_averaged_over_qbin /= delta_q*delta_q*delta_q;

			cout << KTlls[iKT] << "   " << KTuls[iKT] << "   " << bin_centers[ibo] << "   " << bin_centers[ibs] << "   " << bin_centers[ibl] << "   " << CF_averaged_over_qbin << endl;
		}

		delete [] qopts;
		delete [] qowts;
		delete [] qspts;
		delete [] qswts;
		delete [] qlpts;
		delete [] qlwts;
	}

	return 0;
}

//End of file
