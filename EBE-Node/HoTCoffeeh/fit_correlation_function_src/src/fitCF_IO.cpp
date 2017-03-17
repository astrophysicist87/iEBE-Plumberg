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
#include "Arsenal.h"
#include "gauss_quadrature.h"
#include "stats.h"

using namespace std;

void FitCF::replace_parentheses(std::string & tempstring)
{
	int len = tempstring.length();
	for(unsigned int i = 0; i < len; i++)
	{
		char c = tempstring[i];
		if (c == '(' || c == ')')
			tempstring[i] = '_';
	}
	
	if (tempstring[len - 1] == '_')
		tempstring.erase( len - 1 );
	
	return;
}

//allows possibility of dumping thermal_spectra, spectra, log_spectra, etc...
void FitCF::Dump_spectra_array(string output_filename, double *** array_to_dump)
{
	ostringstream filename_stream;
	filename_stream << global_path << "/" << output_filename;
	ofstream out(filename_stream.str().c_str());

	for (int ir = 0; ir < Nparticle; ++ir)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		for (int ipT = 0; ipT < n_interp_pT_pts; ++ipT)
			out << scientific << setprecision(8) << setw(12) << array_to_dump[ir][ipT][ipphi] << "   ";
		out << endl;
	}

	out.close();
}

//allows possibility of reading in thermal_spectra, spectra, log_spectra, etc...
void FitCF::Load_spectra_array(string input_filename, double *** array_to_read)
{
	ostringstream filename_stream;
	filename_stream << global_path << "/" << input_filename;
	ifstream in(filename_stream.str().c_str());

	for (int ir = 0; ir < Nparticle; ++ir)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	for (int ipT = 0; ipT < n_interp_pT_pts; ++ipT)
		in >> array_to_read[ir][ipT][ipphi];

	in.close();
}

void FitCF::Output_results(int mode)
{
	string modeString = "";
	if (mode == 0)
		modeString = "GF_";
	else if (mode == 1)
		modeString = "SVWR_";

	string folderindexstring = "";
	if (currentfolderindex < 0)
		folderindexstring = "avg";
	else
		folderindexstring = patch::to_string(currentfolderindex);

	ostringstream filename_stream_HBT_g0;
	filename_stream_HBT_g0 << global_path << "/HBTradii_" << modeString << "ev" << folderindexstring << no_df_stem << "_grid0.dat";
	ofstream outputHBT_g0;
	outputHBT_g0.open(filename_stream_HBT_g0.str().c_str());
	ostringstream filename_stream_HBT;
	filename_stream_HBT << global_path << "/HBTradii_" << modeString << "ev" << folderindexstring << no_df_stem << ".dat";
	ofstream outputHBT;
	outputHBT.open(filename_stream_HBT.str().c_str());
	ostringstream filename_stream_HBTcfs;
	filename_stream_HBTcfs << global_path << "/HBTradii_" << modeString << "cfs_ev" << folderindexstring << no_df_stem << ".dat";
	ofstream outputHBTcfs;
	outputHBTcfs.open(filename_stream_HBTcfs.str().c_str());

	//at this point, take Chebyshev-spaced pT-pphi grid of R2ij and use to compute R2ij at the KT-Kphi points we want to study
	double flat_R2s[n_interp_pT_pts*n_interp_pphi_pts];
	double flat_R2o[n_interp_pT_pts*n_interp_pphi_pts];
	double flat_R2l[n_interp_pT_pts*n_interp_pphi_pts];
	double flat_R2os[n_interp_pT_pts*n_interp_pphi_pts];
	double flat_R2sl[n_interp_pT_pts*n_interp_pphi_pts];
	double flat_R2ol[n_interp_pT_pts*n_interp_pphi_pts];

	int npts_loc[2] = { n_interp_pT_pts, n_interp_pphi_pts };
	int os[2] = { n_interp_pT_pts-1, n_interp_pphi_pts-1 };
	double lls[2] = { interp_pT_min, interp_pphi_min };
	double uls[2] = { interp_pT_max, interp_pphi_max };
	int modes_loc[2] = { 0, 0 };

	int iptipphi = 0;
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		flat_R2s[iptipphi] = R2_side_GF[ipt][ipphi];
		flat_R2o[iptipphi] = R2_out_GF[ipt][ipphi];
		flat_R2l[iptipphi] = R2_long_GF[ipt][ipphi];
		flat_R2os[iptipphi] = R2_outside_GF[ipt][ipphi];
		flat_R2sl[iptipphi] = R2_sidelong_GF[ipt][ipphi];
		flat_R2ol[iptipphi] = R2_outlong_GF[ipt][ipphi];
		iptipphi++;
	}

	/*approx_R2s = new Chebyshev (flat_R2s, npts_loc, os, lls, uls, 2, modes_loc);
	approx_R2o = new Chebyshev (flat_R2o, npts_loc, os, lls, uls, 2, modes_loc);
	approx_R2l = new Chebyshev (flat_R2l, npts_loc, os, lls, uls, 2, modes_loc);
	approx_R2os = new Chebyshev (flat_R2os, npts_loc, os, lls, uls, 2, modes_loc);
	approx_R2sl = new Chebyshev (flat_R2sl, npts_loc, os, lls, uls, 2, modes_loc);
	approx_R2ol = new Chebyshev (flat_R2ol, npts_loc, os, lls, uls, 2, modes_loc);*/

	//output R2ij on original pT-pphi grid
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		outputHBT_g0 << SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi]
			<< "   " << R2_side_GF[ipt][ipphi] << "   " << R2_out_GF[ipt][ipphi]
			<< "   " << R2_outside_GF[ipt][ipphi] << "   " << R2_long_GF[ipt][ipphi]
			<< "   " << R2_sidelong_GF[ipt][ipphi] << "   " << R2_outlong_GF[ipt][ipphi] << endl;
	}

	const int interpMode = 1;
	//output R2ij them on the desired KT-Kphi grid, and Fourier transform
	for (int iKT = 0; iKT < n_localp_T; ++iKT)
	{
		//output actual extracted R2ij
		for (int iKphi = 0; iKphi < n_localp_phi; ++iKphi)
		{
			double point[2] = { K_T[iKT], K_phi[iKphi] };
			/*outputHBT << K_T[iKT] << "   " << K_phi[iKphi]
				<< "   " << (*approx_R2s).eval(point) << "   " << (*approx_R2o).eval(point)
				<< "   " << (*approx_R2os).eval(point) << "   " << (*approx_R2l).eval(point)
				<< "   " << (*approx_R2sl).eval(point) << "   " << (*approx_R2ol).eval(point) << endl;*/
			outputHBT << K_T[iKT] << "   " << K_phi[iKphi]
				<< "   " << interpolate2D(SPinterp_pT, SPinterp_pphi, R2_side_GF, K_T[iKT], K_phi[iKphi], n_interp_pT_pts, n_interp_pphi_pts, interpMode, false, true)
				<< "   " << interpolate2D(SPinterp_pT, SPinterp_pphi, R2_out_GF, K_T[iKT], K_phi[iKphi], n_interp_pT_pts, n_interp_pphi_pts, interpMode, false, true)
				<< "   " << interpolate2D(SPinterp_pT, SPinterp_pphi, R2_outside_GF, K_T[iKT], K_phi[iKphi], n_interp_pT_pts, n_interp_pphi_pts, interpMode, false, true)
				<< "   " << interpolate2D(SPinterp_pT, SPinterp_pphi, R2_long_GF, K_T[iKT], K_phi[iKphi], n_interp_pT_pts, n_interp_pphi_pts, interpMode, false, true)
				<< "   " << interpolate2D(SPinterp_pT, SPinterp_pphi, R2_sidelong_GF, K_T[iKT], K_phi[iKphi], n_interp_pT_pts, n_interp_pphi_pts, interpMode, false, true)
				<< "   " << interpolate2D(SPinterp_pT, SPinterp_pphi, R2_outlong_GF, K_T[iKT], K_phi[iKphi], n_interp_pT_pts, n_interp_pphi_pts, interpMode, false, true) << endl;
		}

		//do Fourier transforming here for now...
		double plane_psi = 0.0;
		R2_Fourier_transform(iKT, plane_psi, mode);

		//output Fourier coefficients
		if (mode == 0)
		{
			for (int Morder = 0; Morder < n_order; Morder++)
			{
				outputHBTcfs << currentfolderindex << "  " << K_T[iKT] << "  " << Morder
					<< "  " << R2_side_GF_C[iKT][Morder] << "   " << R2_side_GF_S[iKT][Morder] << "  " << R2_out_GF_C[iKT][Morder] << "  " << R2_out_GF_S[iKT][Morder]
					<< "  " << R2_outside_GF_C[iKT][Morder] << "   " << R2_outside_GF_S[iKT][Morder] << "  " << R2_long_GF_C[iKT][Morder] << "  " << R2_long_GF_S[iKT][Morder]
					<< "  " << R2_sidelong_GF_C[iKT][Morder] << "   " << R2_sidelong_GF_S[iKT][Morder] << "  " << R2_outlong_GF_C[iKT][Morder] << "  " << R2_outlong_GF_S[iKT][Morder] << endl;
			}
		}
	}

	/*delete approx_R2s;
	delete approx_R2o;
	delete approx_R2l;
	delete approx_R2os;
	delete approx_R2sl;
	delete approx_R2ol;*/

	outputHBT_g0.close();
	outputHBT.close();
	outputHBTcfs.close();

	return;
}

void FitCF::Output_correlationfunction()
{
	ostringstream oCorrFunc_stream;
	string temp_particle_name = particle_name;
	replace_parentheses(temp_particle_name);

	string CF_reg_string = "";
	if (REGULATE_CF)
		CF_reg_string = "regulated_";

	string CF_proj_string = "";
	if (!FIT_WITH_PROJECTED_CFVALS)
		CF_reg_string = "unprojected_";

	oCorrFunc_stream << global_path << "/correlfunct3D_" << CF_reg_string << CF_proj_string << temp_particle_name << ".dat";
	ofstream oCorrFunc;
	oCorrFunc.open(oCorrFunc_stream.str().c_str());

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	{
		double ckp = cos_SPinterp_pphi[ipphi], skp = sin_SPinterp_pphi[ipphi];
		oCorrFunc << scientific << setprecision(8) << setw(12)
			<< SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi] << "   " << qx_pts[iqx] << "   "
			<< qy_pts[iqy] << "   " << qz_pts[iqz] << "   "
			//<< qx_pts[iqx] * ckp + qy_pts[iqy] * skp << "   "
			//<< -qx_pts[iqx] * skp + qy_pts[iqy] * ckp << "   "
			//<< qz_pts[iqz] << "   "
			<< spectra[target_particle_id][ipt][ipphi] << "   "
			<< thermalCFvals[ipt][ipphi][iqx][iqy][iqz] << "   "
			<< crosstermCFvals[ipt][ipphi][iqx][iqy][iqz] << "   "
			<< resonancesCFvals[ipt][ipphi][iqx][iqy][iqz] << "   "
			<< CFvals[ipt][ipphi][iqx][iqy][iqz] << endl;
	}

	oCorrFunc.close();
				
	return;
}

void FitCF::Output_fleshed_out_correlationfunction(int ipt, int ipphi)
{
	ostringstream oCorrFunc_stream;
	string temp_particle_name = particle_name;
	replace_parentheses(temp_particle_name);
	oCorrFunc_stream << global_path << "/correlfunct3D" << "_" << temp_particle_name << "_fleshed_out.dat";
	ofstream oCorrFunc;
	if (ipt==0 && ipphi==0)
		oCorrFunc.open(oCorrFunc_stream.str().c_str());
	else
		oCorrFunc.open(oCorrFunc_stream.str().c_str(), ios::app);

	for (int iqx = 0; iqx < new_nqpts; ++iqx)
	for (int iqy = 0; iqy < new_nqpts; ++iqy)
	for (int iqz = 0; iqz < new_nqpts; ++iqz)
	{
		double ckp = cos_SPinterp_pphi[ipphi], skp = sin_SPinterp_pphi[ipphi];
		oCorrFunc << scientific << setprecision(7) << setw(15)
			<< SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi] << "   " << qx_fleshed_out_pts[iqx] << "   "
			<< qy_fleshed_out_pts[iqy] << "   " << qz_fleshed_out_pts[iqz] << "   "
			<< qx_fleshed_out_pts[iqx] * ckp + qy_fleshed_out_pts[iqy] * skp << "   "
			<< -qx_fleshed_out_pts[iqx] * skp + qy_fleshed_out_pts[iqy] * ckp << "   "
			<< qz_fleshed_out_pts[iqz] << "   "
			<< fleshed_out_thermal[iqx][iqy][iqz] << "   " << fleshed_out_crossterm[iqx][iqy][iqz] << "   " << fleshed_out_resonances[iqx][iqy][iqz] << "   " << fleshed_out_CF[iqx][iqy][iqz] << endl;
	}

	oCorrFunc.close();
				
	return;
}

void FitCF::Readin_results(int mode)
{
	string modeString = "";
	if (mode == 0)
		modeString = "GF_";

	double dummy;
	ostringstream filename_stream_HBT;
	filename_stream_HBT << global_path << "/HBTradii_" << modeString << "ev" << currentfolderindex << no_df_stem << ".dat";
	ifstream inputHBT(filename_stream_HBT.str().c_str());

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		inputHBT >> dummy;	//pt value
		inputHBT >> dummy;	//pphi value
		inputHBT >> R2_side_GF[ipt][ipphi];
		inputHBT >> R2_out_GF[ipt][ipphi];
		inputHBT >> R2_outside_GF[ipt][ipphi];
		inputHBT >> R2_long_GF[ipt][ipphi];
		inputHBT >> R2_sidelong_GF[ipt][ipphi];
		inputHBT >> R2_outlong_GF[ipt][ipphi];
	}

	inputHBT.close();

	return;
}



void FitCF::Readin_total_target_eiqx_dN_dypTdpTdphi(int folderindex)
{
	string resultsDirectory = "";
	if (currentfolderindex == -1)	//if you're not in any results directory (i.e., if you're averaging), add this
		resultsDirectory = "/results-" + patch::to_string(folderindex);

	string local_name = all_particles[target_particle_id].name;
	replace_parentheses(local_name);
	ostringstream filename_stream_target_dN_dypTdpTdphi;
	filename_stream_target_dN_dypTdpTdphi << global_path << resultsDirectory << "/total_" << local_name << "_eiqx_dN_dypTdpTdphi_ev" << folderindex << no_df_stem << ".dat";
	ifstream input_target_dN_dypTdpTdphi(filename_stream_target_dN_dypTdpTdphi.str().c_str());

	cout << filename_stream_target_dN_dypTdpTdphi.str().c_str() << endl;

	double dummy = 0.0;

	double temp_thermal_spectra[n_interp_pT_pts * n_interp_pphi_pts];
	double temp_spectra[n_interp_pT_pts * n_interp_pphi_pts];

	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		////////////////////////////////////////////////////////////
		//actually start reading stuff here
		//read full spectra
		input_target_dN_dypTdpTdphi >> dummy;
		temp_spectra[ipt * n_interp_pphi_pts + ipphi] = dummy;
		//read full FTd spectra (cos)
		input_target_dN_dypTdpTdphi >> dummy;
		full_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 0)] += dummy;
		//read full FTd spectra (sin)
		input_target_dN_dypTdpTdphi >> dummy;
		full_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 1)] += dummy;
		//read thermal spectra
		input_target_dN_dypTdpTdphi >> dummy;
		temp_thermal_spectra[ipt * n_interp_pphi_pts + ipphi] = dummy;
		//read thermal FTd spectra (cos)
		input_target_dN_dypTdpTdphi >> dummy;
		thermal_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 0)] += dummy;
		//read thermal FTd spectra (sin)
		input_target_dN_dypTdpTdphi >> dummy;
		thermal_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 1)] += dummy;
		//actually finish reading stuff here
		////////////////////////////////////////////////////////////
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;

		/*
		//thermal
		double nonFTd_tspectra = thermal_spectra[target_particle_id][ipt][ipphi];
		double cos_transf_tspectra = thermal_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 0)];
		double sin_transf_tspectra = thermal_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 1)];
		//total
		double nonFTd_spectra = spectra[target_particle_id][ipt][ipphi];
		double cos_transf_spectra = full_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 0)];
		double sin_transf_spectra = full_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 1)];

		if (FIT_WITH_PROJECTED_CFVALS && !thermal_pions_only)
		{
			nonFTd_spectra = nonFTd_tspectra + (nonFTd_spectra - nonFTd_tspectra) / fraction_of_resonances;
			cos_transf_spectra = cos_transf_tspectra + (cos_transf_spectra - cos_transf_tspectra) / fraction_of_resonances;
			sin_transf_spectra = sin_transf_tspectra + (sin_transf_spectra - sin_transf_tspectra) / fraction_of_resonances;
		}

		cout << "CHECK (" << ipt << ", " << ipphi << ", " << iqt << ", " << iqx << ", " << iqy << ", " << iqz << "): "
				<< fraction_of_resonances << "   " << nonFTd_tspectra << "   " << cos_transf_tspectra << "   " << sin_transf_tspectra << "   "
				<< nonFTd_spectra << "   " << cos_transf_spectra << "   " << sin_transf_spectra << endl;
		*/
	}

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		thermal_spectra[target_particle_id][ipt][ipphi] += temp_thermal_spectra[ipt * n_interp_pphi_pts + ipphi];
		spectra[target_particle_id][ipt][ipphi] += temp_spectra[ipt * n_interp_pphi_pts + ipphi];
	}

	input_target_dN_dypTdpTdphi.close();

	return;
}

void FitCF::Readin_total_target_eiqx_dN_dypTdpTdphi(string filename)
{
	ifstream input_target_dN_dypTdpTdphi(filename.c_str());

	//cout << filename_stream_target_dN_dypTdpTdphi.str().c_str() << endl;

	double dummy = 0.0;

	double temp_thermal_spectra[n_interp_pT_pts * n_interp_pphi_pts];
	double temp_spectra[n_interp_pT_pts * n_interp_pphi_pts];

	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		////////////////////////////////////////////////////////////
		//actually start reading stuff here
		//read full spectra
		input_target_dN_dypTdpTdphi >> dummy;
		temp_spectra[ipt * n_interp_pphi_pts + ipphi] = dummy;
		//read full FTd spectra (cos)
		input_target_dN_dypTdpTdphi >> dummy;
		full_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 0)] += dummy;
		//read full FTd spectra (sin)
		input_target_dN_dypTdpTdphi >> dummy;
		full_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 1)] += dummy;
		//read thermal spectra
		input_target_dN_dypTdpTdphi >> dummy;
		temp_thermal_spectra[ipt * n_interp_pphi_pts + ipphi] = dummy;
		//read thermal FTd spectra (cos)
		input_target_dN_dypTdpTdphi >> dummy;
		thermal_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 0)] += dummy;
		//read thermal FTd spectra (sin)
		input_target_dN_dypTdpTdphi >> dummy;
		thermal_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 1)] += dummy;
		//actually finish reading stuff here
		////////////////////////////////////////////////////////////
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
	}

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		thermal_spectra[target_particle_id][ipt][ipphi] += temp_thermal_spectra[ipt * n_interp_pphi_pts + ipphi];
		spectra[target_particle_id][ipt][ipphi] += temp_spectra[ipt * n_interp_pphi_pts + ipphi];
	}

	input_target_dN_dypTdpTdphi.close();

	return;
}


void FitCF::Readin_total_target_eiqx_dN_dypTdpTdphi_evavg()
{
	string local_name = all_particles[target_particle_id].name;
	replace_parentheses(local_name);
	ostringstream filename_stream_target_dN_dypTdpTdphi;
	filename_stream_target_dN_dypTdpTdphi << global_path << "/total_" << local_name << "_eiqx_dN_dypTdpTdphi_evavg" << no_df_stem << ".dat";
	ifstream input_target_dN_dypTdpTdphi(filename_stream_target_dN_dypTdpTdphi.str().c_str());

	cout << filename_stream_target_dN_dypTdpTdphi.str().c_str() << endl;

	double dummy = 0.0;

	double temp_thermal_spectra[n_interp_pT_pts * n_interp_pphi_pts];
	double temp_spectra[n_interp_pT_pts * n_interp_pphi_pts];

	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		for (int iDummy = 0; iDummy < 6; ++iDummy)
			input_target_dN_dypTdpTdphi >> dummy;
		////////////////////////////////////////////////////////////
		//actually start reading stuff here
		input_target_dN_dypTdpTdphi >> temp_spectra[ipt * n_interp_pphi_pts + ipphi];
		input_target_dN_dypTdpTdphi >> full_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 0)];
		input_target_dN_dypTdpTdphi >> full_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 1)];
		input_target_dN_dypTdpTdphi >> temp_thermal_spectra[ipt * n_interp_pphi_pts + ipphi];
		input_target_dN_dypTdpTdphi >> thermal_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 0)];
		input_target_dN_dypTdpTdphi >> thermal_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 1)];
		//actually finish reading stuff here
		////////////////////////////////////////////////////////////
		for (int iDummy = 0; iDummy < 5; ++iDummy)
			input_target_dN_dypTdpTdphi >> dummy;
	}

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		thermal_spectra[target_particle_id][ipt][ipphi] = temp_thermal_spectra[ipt * n_interp_pphi_pts + ipphi];
		spectra[target_particle_id][ipt][ipphi] = temp_spectra[ipt * n_interp_pphi_pts + ipphi];
	}

	input_target_dN_dypTdpTdphi.close();

	return;
}

void FitCF::Readin_resonance_fraction(int folderindex)
{
	string resultsDirectory = "";
	if (currentfolderindex == -1)
		resultsDirectory = "/results-" + patch::to_string(folderindex);

	ostringstream filename_stream_rf;
	filename_stream_rf << global_path << resultsDirectory << "/resonance_fraction.dat";
	ifstream input_rf(filename_stream_rf.str().c_str());

	input_rf >> fraction_of_resonances;

	//cout << "CHECK: fraction_of_resonances = " << fraction_of_resonances << endl;

	input_rf.close();

	return;
}

void FitCF::Output_total_target_eiqx_dN_dypTdpTdphi()
{
	string local_name = all_particles[target_particle_id].name;
	replace_parentheses(local_name);
	ostringstream filename_stream_target_dN_dypTdpTdphi;
	filename_stream_target_dN_dypTdpTdphi << global_path << "/total_" << local_name << "_eiqx_dN_dypTdpTdphi_evavg" << no_df_stem << ".dat";
	ofstream output_target_dN_dypTdpTdphi(filename_stream_target_dN_dypTdpTdphi.str().c_str());

	// addresses NaN issue in sin component when all q^{\mu} == 0
	if (qtnpts%2==1 && qxnpts%2==1 && qynpts%2==1 && qznpts%2==1)
	{	//if all q-ranges are odd and centered on q=0 ==> q=0 is included!
		int iqt0 = (qtnpts-1)/2;
		int iqx0 = (qxnpts-1)/2;
		int iqy0 = (qynpts-1)/2;
		int iqz0 = (qznpts-1)/2;
		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			full_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt0, iqx0, iqy0, iqz0, 1)] = 0.0;
			thermal_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt0, iqx0, iqy0, iqz0, 1)] = 0.0;
		}
	}

	//use this to help look for outliers when CF calculation gets really choppy and noisy
	//Set_target_pphiavgd_CFs();

	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		//first, get CF and projected CF
		double CF = get_CF(ipt, ipphi, iqt, iqx, iqy, iqz, false);				//false means don't return projected value
		//double projected_CF = get_CF(ipt, ipphi, iqt, iqx, iqy, iqz, true && !thermal_pions_only);	//true means do return projected value

		//now, regulate results
		//Regulate_CF(ipt, iqt, iqx, iqy, iqz, &CF, &projected_CF);

		//!!!!!!!!!!!!should get projected_CF AFTER regulating CF...!!!!!!!!!!!!
		double projected_CF = get_CF(ipt, ipphi, iqt, iqx, iqy, iqz, true && !thermal_pions_only);	//true means do return projected value

		double nonFTd_spectra = spectra[target_particle_id][ipt][ipphi];
		double cos_transf_spectra = full_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 0)];
		double sin_transf_spectra = full_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 1)];

		output_target_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12)
			<< qt_pts[iqt] << "   " << qx_pts[iqx] << "   " << qy_pts[iqy] << "   " << qz_pts[iqz] << "   "
			<< SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi] << "   "
			<< nonFTd_spectra << "   "																								//non-thermal + thermal
			<< cos_transf_spectra << "   "																							//non-thermal + thermal (cos)
			<< sin_transf_spectra << "   "																							//non-thermal + thermal (sin)
			<< thermal_spectra[target_particle_id][ipt][ipphi] << "   "																//thermal only
			<< thermal_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 0)] << "   "							//thermal only (cos)
			<< thermal_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 1)] << "   "							//thermal only (sin)
			<< nonFTd_spectra - thermal_spectra[target_particle_id][ipt][ipphi] << "   "											//non-thermal only
			<< cos_transf_spectra - thermal_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 0)] << "   "		//non-thermal only (cos)
			<< sin_transf_spectra - thermal_target_dN_dypTdpTdphi_moments[indexer(ipt, ipphi, iqt, iqx, iqy, iqz, 1)] << "   "		//non-thermal only (sin)
			<< CF << "   " << projected_CF << endl;
	}

	output_target_dN_dypTdpTdphi.close();

	return;
}

void FitCF::Output_lambdas()
{
	ostringstream filename_stream_lambdas;
	filename_stream_lambdas << global_path << "/lambdas.dat";
	ofstream output_lambdas(filename_stream_lambdas.str().c_str());

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		output_lambdas << lambda_Correl[ipt][ipphi] << endl;

	return;
}

//End of file
