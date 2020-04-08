#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>

#include "../include/cfwr.h"
#include "../include/cfwr_lib.h"
#include "../include/Arsenal.h"
#include "../include/gauss_quadrature.h"

using namespace std;

void replace_parentheses(std::string & tempstring)
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
void CorrelationFunction::Dump_FOcells(int local_pid)
{
	ostringstream filename_stream;
	filename_stream << path << "/" << "__res_pid_" << local_pid << "_FOcells.tmp";
	ofstream out(filename_stream.str().c_str(), ios::binary);

	for (int iFOipT = 0; iFOipT < FO_length * n_pT_pts; ++iFOipT)
	{
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
			out << FOcells_to_include[iFOipT][ipphi] << "   ";
		out << endl;
	}

	out.close();
	return;
}

//allows possibility of dumping thermal_spectra, spectra, log_spectra, etc...
void CorrelationFunction::Load_FOcells(int local_pid)
{
	ostringstream filename_stream;
	filename_stream << path << "/" << "__res_pid_" << local_pid << "_FOcells.tmp";
	ifstream in(filename_stream.str().c_str(), ios::binary);

	for (int iFOipT = 0; iFOipT < FO_length * n_pT_pts; ++iFOipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	{
		in >> FOcells_to_include[iFOipT][ipphi];
		/*if (FOcells_to_include[iFOipT][ipphi] < 0)
		{
			cerr << "Shouldn't have stored a negative value: " << iFOipT << "   " << ipphi << "   " << FOcells_to_include[iFOipT][ipphi] << endl;
			exit(1);
		}*/
	}

	in.close();
	return;
}

//allows possibility of dumping thermal_spectra, spectra, log_spectra, etc...
void CorrelationFunction::Dump_spectra_array(string output_filename, double *** array_to_dump)
{
	ostringstream filename_stream;
	filename_stream << path << "/" << output_filename;
	ofstream out(filename_stream.str().c_str());

	for (int ir = 0; ir < Nparticle; ++ir)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	{
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		{
			//cout << output_filename << ": " << ir << "   " << ipT << "   " << ipphi << endl;
			out << scientific << setprecision(8) << setw(12) << array_to_dump[ir][ipT][ipphi] << "   ";
		}
		out << endl;
	}

	out.close();
}

//allows possibility of reading in thermal_spectra, spectra, log_spectra, etc...
void CorrelationFunction::Load_spectra_array(string input_filename, double *** array_to_read)
{
	ostringstream filename_stream;
	filename_stream << path << "/" << input_filename;
	ifstream in(filename_stream.str().c_str());

	for (int ir = 0; ir < Nparticle; ++ir)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		in >> array_to_read[ir][ipT][ipphi];

	in.close();
}

void CorrelationFunction::Output_results(int mode)
{
	string modeString = "";
	if (mode == 0)
		modeString = "GF";
	else if (mode == 1)
		modeString = "QM";

	ostringstream filename_stream_HBT_g0;
	filename_stream_HBT_g0 << path << "/HBTradii_" << modeString << no_df_stem << "_grid0.dat";
	ofstream outputHBT_g0;
	outputHBT_g0.open(filename_stream_HBT_g0.str().c_str());
	ostringstream filename_stream_HBT;
	filename_stream_HBT << path << "/HBTradii_" << modeString << no_df_stem << ".dat";
	ofstream outputHBT;
	outputHBT.open(filename_stream_HBT.str().c_str());
	ostringstream filename_stream_HBTcfs;
	filename_stream_HBTcfs << path << "/HBTradii_" << modeString << "_cfs" << no_df_stem << ".dat";
	ofstream outputHBTcfs;
	outputHBTcfs.open(filename_stream_HBTcfs.str().c_str());

	//at this point, take Chebyshev-spaced pT-pphi grid of R2ij and use to compute R2ij at the KT-Kphi points we want to study
	double flat_R2s[n_pT_pts*n_pphi_pts];
	double flat_R2o[n_pT_pts*n_pphi_pts];
	double flat_R2l[n_pT_pts*n_pphi_pts];
	double flat_R2os[n_pT_pts*n_pphi_pts];
	double flat_R2sl[n_pT_pts*n_pphi_pts];
	double flat_R2ol[n_pT_pts*n_pphi_pts];

	int npts_loc[2] = { n_pT_pts, n_pphi_pts };
	int os[2] = { n_pT_pts-1, n_pphi_pts-1 };
	double lls[2] = { KT_min, Kphi_min };
	double uls[2] = { KT_max, Kphi_max };
	int modes_loc[2] = { 0, 0 };

	int iptipphi = 0;
	if (mode == 0)
	{
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			flat_R2s[iptipphi] = R2_side_GF[ipT][ipphi];
			flat_R2o[iptipphi] = R2_out_GF[ipT][ipphi];
			flat_R2l[iptipphi] = R2_long_GF[ipT][ipphi];
			flat_R2os[iptipphi] = R2_outside_GF[ipT][ipphi];
			flat_R2sl[iptipphi] = R2_sidelong_GF[ipT][ipphi];
			flat_R2ol[iptipphi] = R2_outlong_GF[ipT][ipphi];
			iptipphi++;
		}
	}
	else if (mode == 1)
	{
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			flat_R2s[iptipphi] = R2_side_QM[ipT][ipphi];
			flat_R2o[iptipphi] = R2_out_QM[ipT][ipphi];
			flat_R2l[iptipphi] = R2_long_QM[ipT][ipphi];
			flat_R2os[iptipphi] = R2_outside_QM[ipT][ipphi];
			flat_R2sl[iptipphi] = R2_sidelong_QM[ipT][ipphi];
			flat_R2ol[iptipphi] = R2_outlong_QM[ipT][ipphi];
			iptipphi++;
		}
	}

	//output R2ij on original pT-pphi grid
	if (mode == 0)
	{
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			outputHBT_g0 << SP_pT[ipT] << "   " << SP_pphi[ipphi]
				<< "   " << R2_side_GF[ipT][ipphi] << "   " << R2_out_GF[ipT][ipphi]
				<< "   " << R2_outside_GF[ipT][ipphi] << "   " << R2_long_GF[ipT][ipphi]
				<< "   " << R2_sidelong_GF[ipT][ipphi] << "   " << R2_outlong_GF[ipT][ipphi] << endl;
		}
	}
	else if (mode == 1)
	{
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			outputHBT_g0 << SP_pT[ipT] << "   " << SP_pphi[ipphi]
				<< "   " << R2_side_QM[ipT][ipphi] << "   " << R2_out_QM[ipT][ipphi]
				<< "   " << R2_outside_QM[ipT][ipphi] << "   " << R2_long_QM[ipT][ipphi]
				<< "   " << R2_sidelong_QM[ipT][ipphi] << "   " << R2_outlong_QM[ipT][ipphi] << endl;
		}
	}

	const int interpMode = 1;
	//output R2ij them on the desired KT-Kphi grid, and Fourier transform
	for (int iKT = 0; iKT < nKT; ++iKT)
	{
		//output actual extracted R2ij
		for (int iKphi = 0; iKphi < nKphi; ++iKphi)
		{
			double point[2] = { K_T[iKT], K_phi[iKphi] };
			outputHBT << K_T[iKT] << "   " << K_phi[iKphi]
				<< "   " << interpolate2D(SP_pT, SP_pphi, R2_side_GF, K_T[iKT], K_phi[iKphi], n_pT_pts, n_pphi_pts, interpMode, false, true)
				<< "   " << interpolate2D(SP_pT, SP_pphi, R2_out_GF, K_T[iKT], K_phi[iKphi], n_pT_pts, n_pphi_pts, interpMode, false, true)
				<< "   " << interpolate2D(SP_pT, SP_pphi, R2_outside_GF, K_T[iKT], K_phi[iKphi], n_pT_pts, n_pphi_pts, interpMode, false, true)
				<< "   " << interpolate2D(SP_pT, SP_pphi, R2_long_GF, K_T[iKT], K_phi[iKphi], n_pT_pts, n_pphi_pts, interpMode, false, true)
				<< "   " << interpolate2D(SP_pT, SP_pphi, R2_sidelong_GF, K_T[iKT], K_phi[iKphi], n_pT_pts, n_pphi_pts, interpMode, false, true)
				<< "   " << interpolate2D(SP_pT, SP_pphi, R2_outlong_GF, K_T[iKT], K_phi[iKphi], n_pT_pts, n_pphi_pts, interpMode, false, true) << endl;
		}

		//do Fourier transforming here for now...
		double plane_psi = global_plane_psi;
		R2_Fourier_transform(iKT, plane_psi, mode);

		//output Fourier coefficients
		if (mode == 0)
		{
			for (int Morder = 0; Morder < n_order; Morder++)
			{
				outputHBTcfs << K_T[iKT] << "  " << Morder
					<< "  " << R2_side_GF_C[iKT][Morder]
					<< "   " << R2_side_GF_S[iKT][Morder]
					<< "  " << R2_out_GF_C[iKT][Morder]
					<< "  " << R2_out_GF_S[iKT][Morder]
					<< "  " << R2_outside_GF_C[iKT][Morder]
					<< "   " << R2_outside_GF_S[iKT][Morder]
					<< "  " << R2_long_GF_C[iKT][Morder]
					<< "  " << R2_long_GF_S[iKT][Morder]
					<< "  " << R2_sidelong_GF_C[iKT][Morder]
					<< "   " << R2_sidelong_GF_S[iKT][Morder]
					<< "  " << R2_outlong_GF_C[iKT][Morder]
					<< "  " << R2_outlong_GF_S[iKT][Morder] << endl;
			}
		}
		else if (mode == 1)
		{
			for (int Morder = 0; Morder < n_order; Morder++)
			{
				outputHBTcfs << SP_pT[iKT] << "  " << Morder
					<< "  " << R2_side_QM_C[iKT][Morder]
					<< "   " << R2_side_QM_S[iKT][Morder]
					<< "  " << R2_out_QM_C[iKT][Morder]
					<< "  " << R2_out_QM_S[iKT][Morder]
					<< "  " << R2_outside_QM_C[iKT][Morder]
					<< "   " << R2_outside_QM_S[iKT][Morder]
					<< "  " << R2_long_QM_C[iKT][Morder]
					<< "  " << R2_long_QM_S[iKT][Morder]
					<< "  " << R2_sidelong_QM_C[iKT][Morder]
					<< "   " << R2_sidelong_QM_S[iKT][Morder]
					<< "  " << R2_outlong_QM_C[iKT][Morder]
					<< "  " << R2_outlong_QM_S[iKT][Morder] << endl;
			}
		}
	}

	outputHBT_g0.close();
	outputHBT.close();
	outputHBTcfs.close();

	return;
}

void CorrelationFunction::Output_lambdas()
{
	ostringstream filename_stream_lambdas;
	filename_stream_lambdas << path << "/lambdas.dat";
	ofstream outputlambdas;
	outputlambdas.open(filename_stream_lambdas.str().c_str());
	
	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		outputlambdas << SP_pT[ipT] << "   " << SP_pphi[ipphi] << "   " << lambda_Correl[ipT][ipphi] << endl;

	outputlambdas.close();
	return;
}

void CorrelationFunction::Output_correlationfunction(bool project_CF /*==true*/)
{
	ostringstream oCorrFunc_stream;
	string temp_particle_name = particle_name;
	replace_parentheses(temp_particle_name);

	string CF_proj_string = "";
	if (!project_CF)
		CF_proj_string = "unprojected_";

	oCorrFunc_stream << path << "/correlfunct3D_" << CF_proj_string << temp_particle_name << ".dat";
	ofstream oCorrFunc;
	oCorrFunc.open(oCorrFunc_stream.str().c_str());

	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	{
		double ckp = cos_SP_pphi[ipphi], skp = sin_SP_pphi[ipphi];
		oCorrFunc << scientific << setprecision(8) << setw(12)
			<< SP_pT[ipT] << "   " << SP_pphi[ipphi] << "   "
			<< setprecision(3) << setw(5) << qx_pts[iqx] << "   "
			<< qy_pts[iqy] << "   " << qz_pts[iqz] << "   "
			//<< qx_pts[iqx] * ckp + qy_pts[iqy] * skp << "   "
			//<< -qx_pts[iqx] * skp + qy_pts[iqy] * ckp << "   "
			//<< qz_pts[iqz] << "   "
			<< setprecision(6) << setw(10)
			<< spectra[target_particle_id][ipT][ipphi] << "   "
			<< thermalCFvals[ipT][ipphi][iqx][iqy][iqz] << "   "
			<< crosstermCFvals[ipT][ipphi][iqx][iqy][iqz] << "   "
			<< resonancesCFvals[ipT][ipphi][iqx][iqy][iqz] << "   "
			<< CFvals[ipT][ipphi][iqx][iqy][iqz] << endl;
	}

	oCorrFunc.close();
				
	return;
}

void CorrelationFunction::Output_fleshed_out_correlationfunction(int ipT, int ipphi, bool project_CF /*==true*/)
{
	ostringstream oCorrFunc_stream;
	string temp_particle_name = particle_name;
	replace_parentheses(temp_particle_name);

	string CF_proj_string = "";
	if (!project_CF)
		CF_proj_string = "unprojected_";

	oCorrFunc_stream << path << "/correlfunct3D_" << CF_proj_string << temp_particle_name << "_fleshed_out.dat";
	ofstream oCorrFunc;
	if (ipT==0 && ipphi==0)
		oCorrFunc.open(oCorrFunc_stream.str().c_str());
	else
		oCorrFunc.open(oCorrFunc_stream.str().c_str(), ios::app);

	if (SLICE_OF_FLESH_ONLY)
	{
		for (int iqx = 0; iqx < new_nqxpts; ++iqx)
		{
			double ckp = cos_SP_pphi[ipphi], skp = sin_SP_pphi[ipphi];
			int iqy = (new_nqpts-1)/2;
			int iqz = (new_nqpts-1)/2;
			oCorrFunc << scientific << setprecision(7) << setw(15)
				<< SP_pT[ipT] << "   "
				<< SP_pphi[ipphi] << "   "
				<< qx_fleshed_out_pts[iqx] << "   "
				<< qy_fleshed_out_pts[iqy] << "   "
				<< qz_fleshed_out_pts[iqz] << "   "
				<< qx_fleshed_out_pts[iqx] * ckp + qy_fleshed_out_pts[iqy] * skp << "   "
				<< -qx_fleshed_out_pts[iqx] * skp + qy_fleshed_out_pts[iqy] * ckp << "   "
				<< qz_fleshed_out_pts[iqz] << "   "
				<< fleshed_out_thermal[iqx][iqy][iqz] << "   "
				<< fleshed_out_crossterm[iqx][iqy][iqz] << "   "
				<< fleshed_out_resonances[iqx][iqy][iqz] << "   "
				<< fleshed_out_CF[iqx][iqy][iqz] << endl;
		}
		for (int iqy = 0; iqy < new_nqypts; ++iqy)
		{
			double ckp = cos_SP_pphi[ipphi], skp = sin_SP_pphi[ipphi];
			int iqx = (new_nqpts-1)/2;
			int iqz = (new_nqpts-1)/2;
			oCorrFunc << scientific << setprecision(7) << setw(15)
				<< SP_pT[ipT] << "   "
				<< SP_pphi[ipphi] << "   "
				<< qx_fleshed_out_pts[iqx] << "   "
				<< qy_fleshed_out_pts[iqy] << "   "
				<< qz_fleshed_out_pts[iqz] << "   "
				<< qx_fleshed_out_pts[iqx] * ckp + qy_fleshed_out_pts[iqy] * skp << "   "
				<< -qx_fleshed_out_pts[iqx] * skp + qy_fleshed_out_pts[iqy] * ckp << "   "
				<< qz_fleshed_out_pts[iqz] << "   "
				<< fleshed_out_thermal[iqx][iqy][iqz] << "   "
				<< fleshed_out_crossterm[iqx][iqy][iqz] << "   "
				<< fleshed_out_resonances[iqx][iqy][iqz] << "   "
				<< fleshed_out_CF[iqx][iqy][iqz] << endl;
		}
		for (int iqz = 0; iqz < new_nqzpts; ++iqz)
		{
			double ckp = cos_SP_pphi[ipphi], skp = sin_SP_pphi[ipphi];
			int iqx = (new_nqpts-1)/2;
			int iqy = (new_nqpts-1)/2;
			oCorrFunc << scientific << setprecision(7) << setw(15)
				<< SP_pT[ipT] << "   "
				<< SP_pphi[ipphi] << "   "
				<< qx_fleshed_out_pts[iqx] << "   "
				<< qy_fleshed_out_pts[iqy] << "   "
				<< qz_fleshed_out_pts[iqz] << "   "
				<< qx_fleshed_out_pts[iqx] * ckp + qy_fleshed_out_pts[iqy] * skp << "   "
				<< -qx_fleshed_out_pts[iqx] * skp + qy_fleshed_out_pts[iqy] * ckp << "   "
				<< qz_fleshed_out_pts[iqz] << "   "
				<< fleshed_out_thermal[iqx][iqy][iqz] << "   "
				<< fleshed_out_crossterm[iqx][iqy][iqz] << "   "
				<< fleshed_out_resonances[iqx][iqy][iqz] << "   "
				<< fleshed_out_CF[iqx][iqy][iqz] << endl;
		}
	}
	else
	{
		for (int iqx = 0; iqx < new_nqxpts; ++iqx)
		for (int iqy = 0; iqy < new_nqypts; ++iqy)
		for (int iqz = 0; iqz < new_nqzpts; ++iqz)
		{
			double ckp = cos_SP_pphi[ipphi], skp = sin_SP_pphi[ipphi];
			oCorrFunc << scientific << setprecision(7) << setw(15)
				<< SP_pT[ipT] << "   " << SP_pphi[ipphi] << "   " << qx_fleshed_out_pts[iqx] << "   "
				<< qy_fleshed_out_pts[iqy] << "   " << qz_fleshed_out_pts[iqz] << "   "
				<< qx_fleshed_out_pts[iqx] * ckp + qy_fleshed_out_pts[iqy] * skp << "   "
				<< -qx_fleshed_out_pts[iqx] * skp + qy_fleshed_out_pts[iqy] * ckp << "   "
				<< qz_fleshed_out_pts[iqz] << "   "
				<< fleshed_out_thermal[iqx][iqy][iqz] << "   "
				<< fleshed_out_crossterm[iqx][iqy][iqz] << "   "
				<< fleshed_out_resonances[iqx][iqy][iqz] << "   "
				<< fleshed_out_CF[iqx][iqy][iqz] << endl;
		}
	}

	oCorrFunc.close();
				
	return;
}

void CorrelationFunction::Readin_results(int mode)
{
	string modeString = "";
	if (mode == 0)
		modeString = "GF";
	else if (mode == 1)
		modeString = "QM";

	double dummy;
	ostringstream filename_stream_HBT;
	filename_stream_HBT << path << "/HBTradii_" << modeString << no_df_stem << ".dat";
	ifstream inputHBT(filename_stream_HBT.str().c_str());

	if (mode == 0)
	{
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			inputHBT >> dummy;	//pt value
			inputHBT >> dummy;	//pphi value
			inputHBT >> R2_side_GF[ipT][ipphi];
			inputHBT >> R2_out_GF[ipT][ipphi];
			inputHBT >> R2_outside_GF[ipT][ipphi];
			inputHBT >> R2_long_GF[ipT][ipphi];
			inputHBT >> R2_sidelong_GF[ipT][ipphi];
			inputHBT >> R2_outlong_GF[ipT][ipphi];
		}
	}
	else if (mode == 1)
	{
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			inputHBT >> dummy;	//pt value
			inputHBT >> dummy;	//pphi value
			inputHBT >> R2_side_QM[ipT][ipphi];
			inputHBT >> R2_out_QM[ipT][ipphi];
			inputHBT >> R2_outside_QM[ipT][ipphi];
			inputHBT >> R2_long_QM[ipT][ipphi];
			inputHBT >> R2_sidelong_QM[ipT][ipphi];
			inputHBT >> R2_outlong_QM[ipT][ipphi];
		}
	}

	inputHBT.close();

	return;
}

void CorrelationFunction::Output_total_target_dN_dypTdpTdphi()
{
	string local_name = all_particles[target_particle_id].name;
	replace_parentheses(local_name);
	ostringstream filename_stream_target_dN_dypTdpTdphi;
	filename_stream_target_dN_dypTdpTdphi << path << "/total_" << local_name << "_dN_dypTdpTdphi.dat";
	ofstream output_target_dN_dypTdpTdphi(filename_stream_target_dN_dypTdpTdphi.str().c_str());

	for(int ipphi = 0; ipphi < n_pphi_pts; ipphi++)
	{
		for(int ipT = 0; ipT < n_pT_pts; ipT++)
		{
			double fs = spectra[target_particle_id][ipT][ipphi];
			double ts = thermal_spectra[target_particle_id][ipT][ipphi];
			double result = (fs - ts) / fraction_of_resonances + ts;
			output_target_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12) << spectra[target_particle_id][ipT][ipphi] << "   ";
		}
		output_target_dN_dypTdpTdphi << endl;
	}

	output_target_dN_dypTdpTdphi.close();

	return;
}

void CorrelationFunction::Output_total_target_eiqx_dN_dypTdpTdphi(double current_fraction /*==-1.0*/)
{
	string local_name = all_particles[target_particle_id].name;
	string current_fraction_string = "";
	replace_parentheses(local_name);
	ostringstream filename_stream_target_dN_dypTdpTdphi;
	filename_stream_target_dN_dypTdpTdphi << path << "/total_" << local_name << current_fraction_string << "_eiqx_dN_dypTdpTdphi" << no_df_stem << ".dat";
	ofstream output_target_dN_dypTdpTdphi(filename_stream_target_dN_dypTdpTdphi.str().c_str());

	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	{
		//first, get CF and projected CF
		double CF = get_CF(ipT, ipphi, iqt, iqx, iqy, iqz, false);				//false means don't return projected value

		//!!!!!!!!!!!!should get projected_CF AFTER regulating CF...!!!!!!!!!!!!
		double projected_CF = get_CF(ipT, ipphi, iqt, iqx, iqy, iqz, true && !thermal_pions_only);	//true means do return projected value

		double nonFTd_spectra = spectra[target_particle_id][ipT][ipphi];
		//
		double cos_transf_spectra = full_target_Yeq0_moments[indexer(ipT,ipphi,iqt,iqx,iqy,iqz,0)]
										+ full_target_Yeq0_moments[indexer(ipT,ipphi,iqt,iqx,iqy,iqz,3)];
		//
		double sin_transf_spectra = full_target_Yeq0_moments[indexer(ipT,ipphi,iqt,iqx,iqy,iqz,1)]
										+ full_target_Yeq0_moments[indexer(ipT,ipphi,iqt,iqx,iqy,iqz,2)];
		//
		double cos_transf_tspectra = thermal_target_Yeq0_moments[indexer(ipT,ipphi,iqt,iqx,iqy,iqz,0)]
										+ thermal_target_Yeq0_moments[indexer(ipT,ipphi,iqt,iqx,iqy,iqz,3)];
		//
		double sin_transf_tspectra = thermal_target_Yeq0_moments[indexer(ipT,ipphi,iqt,iqx,iqy,iqz,1)]
										+ thermal_target_Yeq0_moments[indexer(ipT,ipphi,iqt,iqx,iqy,iqz,2)];
		//
		/*cout << "IOcheck: " << qt_pts[iqt] << "   " << qx_pts[iqx] << "   " << qy_pts[iqy] << "   " << qz_pts[iqz] << "   "
                        << SP_pT[ipT] << "   " << SP_pphi[ipphi] << endl
						<< "\t\t" << thermal_target_Yeq0_moments[indexer(ipT,ipphi,iqt,iqx,iqy,iqz,0)] << "   "
                        << thermal_target_Yeq0_moments[indexer(ipT,ipphi,iqt,iqx,iqy,iqz,1)] << "   "
                        << thermal_target_Yeq0_moments[indexer(ipT,ipphi,iqt,iqx,iqy,iqz,2)] << "   "
                        << thermal_target_Yeq0_moments[indexer(ipT,ipphi,iqt,iqx,iqy,iqz,3)] << endl
						<< "\t\t" << full_target_Yeq0_moments[indexer(ipT,ipphi,iqt,iqx,iqy,iqz,0)] << "   "
                        << full_target_Yeq0_moments[indexer(ipT,ipphi,iqt,iqx,iqy,iqz,1)] << "   "
                        << full_target_Yeq0_moments[indexer(ipT,ipphi,iqt,iqx,iqy,iqz,2)] << "   "
                        << full_target_Yeq0_moments[indexer(ipT,ipphi,iqt,iqx,iqy,iqz,3)] << endl;*/

		output_target_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12)
			<< qt_pts[iqt] << "   " << qx_pts[iqx] << "   " << qy_pts[iqy] << "   " << qz_pts[iqz] << "   "
			<< SP_pT[ipT] << "   " << SP_pphi[ipphi] << "   "
			<< nonFTd_spectra << "   "														//non-thermal + thermal
			<< cos_transf_spectra << "   "													//non-thermal + thermal (cos)
			<< sin_transf_spectra << "   "													//non-thermal + thermal (sin)
			<< thermal_spectra[target_particle_id][ipT][ipphi] << "   "						//thermal only
			<< cos_transf_tspectra << "   "													//thermal only (cos)
			<< sin_transf_tspectra << "   "													//thermal only (sin)
			<< nonFTd_spectra - thermal_spectra[target_particle_id][ipT][ipphi] << "   "	//non-thermal only
			<< cos_transf_spectra - cos_transf_tspectra << "   "							//non-thermal only (cos)
			<< sin_transf_spectra - sin_transf_tspectra << "   "							//non-thermal only (sin)
			<< CF << "   " << projected_CF << endl;
	}

	output_target_dN_dypTdpTdphi.close();

	return;
}

void CorrelationFunction::Output_thermal_target_eiqx_dN_dypTdpTdphi(int iqt, int iqz)
{
	string local_name = all_particles[target_particle_id].name;
	replace_parentheses(local_name);
	ostringstream filename_stream_target_dN_dypTdpTdphi;
	filename_stream_target_dN_dypTdpTdphi << path << "/thermal_" << local_name << "_iqt_" << iqt << "_iqz_" << iqz << "_eiqx_dN_dypTdpTdphi.dat";
	ofstream output_target_dN_dypTdpTdphi(filename_stream_target_dN_dypTdpTdphi.str().c_str());

	int HDFOpenSuccess = Administrate_target_thermal_HDF_array(1);	// 1 - open

	int accessHDFresonanceSpectra = Access_target_thermal_in_HDF_array(iqt, iqz, 1, current_dN_dypTdpTdphi_moments);		//get

	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int ipY = 0; ipY < n_pY_pts; ++ipY)
	{
		double loc_qt = qt_pts[iqt];
		double loc_qz = qz_pts[iqz];
		current_pY_shift = 0.5 * log(abs((loc_qt+loc_qz + 1.e-100)/(loc_qt-loc_qz + 1.e-100)));

		output_target_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12)
			<< qt_pts[iqt] << "   " << qx_pts[iqx] << "   " << qy_pts[iqy] << "   " << qz_pts[iqz] << "   "
			<< SP_pT[ipT] << "   " << SP_pphi[ipphi] << "   " << current_pY_shift + SP_Del_pY[ipY] << "   " << SP_Del_pY[ipY] << "   "
			<< current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,0)] << "   "
			<< current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,1)] << "   "
			<< current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,2)] << "   "
			<< current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,3)] << endl;
	}

	int HDFCloseSuccess = Administrate_target_thermal_HDF_array(2);	// 2 - close

	output_target_dN_dypTdpTdphi.close();

	return;
}

void CorrelationFunction::Output_total_eiqx_dN_dypTdpTdphi(int local_pid, int iqt, int iqz)
{
	string local_name = all_particles[local_pid].name;
	replace_parentheses(local_name);
	ostringstream filename_stream_dN_dypTdpTdphi;
	filename_stream_dN_dypTdpTdphi << path << "/total_" << local_name << "_iqt_" << iqt << "_iqz_" << iqz << "_eiqx_dN_dypTdpTdphi.dat";
	ofstream output_dN_dypTdpTdphi(filename_stream_dN_dypTdpTdphi.str().c_str());

	int HDFOpenSuccess = Administrate_resonance_HDF_array(1);	// 1 - open

	int accessHDFresonanceSpectra = Access_resonance_in_HDF_array(local_pid, iqt, iqz, 1, current_dN_dypTdpTdphi_moments);		//get

	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int ipY = 0; ipY < n_pY_pts; ++ipY)
	{
		double loc_qt = qt_pts[iqt];
		double loc_qz = qz_pts[iqz];
		current_pY_shift = 0.5 * log(abs((loc_qt+loc_qz + 1.e-100)/(loc_qt-loc_qz + 1.e-100)));

		double nonFTd_spectra = spectra[local_pid][ipT][ipphi];
		double cos_transf_spectra = current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,0)]
									+ current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,3)];
		double sin_transf_spectra = current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,1)]
									+ current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,2)];

		output_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12)
			<< qt_pts[iqt] << "   " << qx_pts[iqx] << "   " << qy_pts[iqy] << "   " << qz_pts[iqz] << "   "
			<< SP_pT[ipT] << "   " << SP_pphi[ipphi] << "   " << current_pY_shift + SP_Del_pY[ipY] << "   " << SP_Del_pY[ipY] << "   "
			<< current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,0)] << "   "
			<< current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,1)] << "   "
			<< current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,2)] << "   "
			<< current_dN_dypTdpTdphi_moments[fixQTQZ_indexer(ipT,ipphi,ipY,iqx,iqy,3)] << endl;
	}

	int HDFCloseSuccess = Administrate_resonance_HDF_array(2);	// 2 - close

	output_dN_dypTdpTdphi.close();

	return;
}

void CorrelationFunction::Readin_total_target_eiqx_dN_dypTdpTdphi()
{
	if (1)
	{
		cerr << "This function needs to be re-written!  Exiting..." << endl;
		exit(8);
	}

	string local_name = all_particles[target_particle_id].name;
	replace_parentheses(local_name);
	ostringstream filename_stream_target_dN_dypTdpTdphi;
	filename_stream_target_dN_dypTdpTdphi << path << "/total_" << local_name << "_eiqx_dN_dypTdpTdphi_" << no_df_stem << ".dat";
	ifstream input_target_dN_dypTdpTdphi(filename_stream_target_dN_dypTdpTdphi.str().c_str());

	double dummy = 0.0;

	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	{
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> spectra[target_particle_id][ipT][ipphi];
		input_target_dN_dypTdpTdphi >> current_dN_dypTdpTdphi_moments[indexer(ipT,ipphi,iqt,iqx,iqy,iqz,0)];
		input_target_dN_dypTdpTdphi >> current_dN_dypTdpTdphi_moments[indexer(ipT,ipphi,iqt,iqx,iqy,iqz,1)];
		input_target_dN_dypTdpTdphi >> thermal_spectra[target_particle_id][ipT][ipphi];
		input_target_dN_dypTdpTdphi >> thermal_target_dN_dypTdpTdphi_moments[indexer(ipT,ipphi,iqt,iqx,iqy,iqz,0)];
		input_target_dN_dypTdpTdphi >> thermal_target_dN_dypTdpTdphi_moments[indexer(ipT,ipphi,iqt,iqx,iqy,iqz,1)];
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
	}

	input_target_dN_dypTdpTdphi.close();

	return;
}

void CorrelationFunction::Read_in_correlationfunction()
{
	ostringstream iCorrFunc_stream;
	string temp_particle_name = particle_name;
	replace_parentheses(temp_particle_name);

	string CF_proj_string = "";
	if (!FIT_WITH_PROJECTED_CFVALS)
		CF_proj_string = "unprojected_";

	iCorrFunc_stream << path << "/correlfunct3D_" << CF_proj_string << temp_particle_name << ".dat";
	ifstream iCorrFunc;
	iCorrFunc.open(iCorrFunc_stream.str().c_str());

	double dummy;
	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	{
		iCorrFunc >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
					>> thermalCFvals[ipT][ipphi][iqx][iqy][iqz] >> crosstermCFvals[ipT][ipphi][iqx][iqy][iqz] >> resonancesCFvals[ipT][ipphi][iqx][iqy][iqz] >> CFvals[ipT][ipphi][iqx][iqy][iqz];
		//*out << "Read in CFvals[" << ipT << "][" << ipphi << "][" << iqx << "][" << iqy << "][" << iqz << "] = " << CFvals[ipT][ipphi][iqx][iqy][iqz] << endl;
	}

	iCorrFunc.close();
				
	return;
}

void CorrelationFunction::Output_chosen_resonances()
{
	ostringstream filename_stream_crf;
	filename_stream_crf << path << "/chosen_resonances.dat";
	ofstream output_crf(filename_stream_crf.str().c_str());

	output_crf << particle_monval << endl;
	for (int icr = 0; icr < (int)chosen_resonances.size(); icr++)
		output_crf << all_particles[chosen_resonances[icr]].monval << endl;

	output_crf.close();

	return;
}

void CorrelationFunction::Output_resonance_fraction()
{
	ostringstream filename_stream_rf;
	filename_stream_rf << path << "/resonance_fraction.dat";
	ofstream output_rf(filename_stream_rf.str().c_str());

	output_rf << fraction_of_resonances << endl;

	output_rf.close();

	return;
}

/*
void CorrelationFunction::Output_parameters_record()
{
	ostringstream filename_stream_HCparams;
	filename_stream_HCparams << path << "/HotCoffeeh_parameters.dat";
	ofstream output_HCparams(filename_stream_HCparams.str().c_str());

	string tentab = "\t\t\t\t\t\t\t\t\t\t";

	output_HCparams
		<< "#############################################" << endl
		<< "## Default parameters for HoTCoffeeh codes ##" << endl
		<< "#############################################" << endl << endl

		<< "###########################" << endl
		<< "## SVWR and CFWR parameters" << endl << endl;
		
	output_HCparams
		<< "# Model options" << endl;

	output_HCparams
		<< "grouping_particles = " << GROUPING_PARTICLES << "\t\t\t\t\t# set to 1 to perform calculations for similar" << endl
		<< tentab << "# particles together" << endl << endl;

	output_HCparams
		<< "use_plane_psi_order = " << USE_PLANE_PSI_ORDER << "\t\t\t\t\t# specifies whether to do HBT relative to flow-plane angle," << endl
		<< tentab << "# and at what order." << endl
		<< tentab << "#  0 - use plane_psi = 0.0" << endl
		<< tentab << "# !0 - use flow-plane angle at given order" << endl << endl;

	output_HCparams
		<< "chosenParticlesMode = " << chosenParticlesMode << "\t\t\t\t\t# 0 - set with threshold" << endl
		<< tentab << "# 1 - read in chosen particles from file" << endl << endl;

	output_HCparams
		<< "ignore_long_lived_resonances = " << IGNORE_LONG_LIVED_RESONANCES << "\t\t# do all decay channels (or don't)" << endl << endl;

	output_HCparams
		<< "include_delta_f	= " << INCLUDE_DELTA_F << "\t\t\t\t\t\t# include viscous CF corrections" << endl << endl;

	output_HCparams
		<< "flag_negative_S = " << flagneg << "\t\t\t\t\t\t# neglect all FO cells with negative emission function" << endl
		<< tentab << "# choose flag_negative_S == 0 to agree with iS.e" << endl
		<< tentab << "# choose flag_negative_S == 1 for the real world" << endl << endl;

	output_HCparams
		<< "include_bulk_pi = " << INCLUDE_BULK_PI << "\t\t\t\t\t\t# 1 means that there's an extra column in decdat2.dat," << endl
		<< tentab << "# i.e., bulk pi is NOT(!!!) used in current calculation!!!" << endl << endl;

	output_HCparams
		<< "# Model parameters " << endl
		<< "n_order = " << n_order << "\t\t\t\t\t\t\t\t# sets max order for azimuthal averaging (0...nOrder-1)" << endl << endl;

	output_HCparams
		<< "tolerance = " << tol << "\t\t\t\t\t\t# tolerance for flagging negative FO cells" << endl << endl;

	output_HCparams
		<< "max_lifetime = " << max_lifetime << "\t\t\t\t\t# resonance is ignored if its lifetime is longer than this value" << endl
		<< tentab << "# AND ignore_long_lived_resonances is set to 1" << endl << endl;

	output_HCparams
		<< "particle_diff_tolerance = " << PARTICLE_DIFF_TOLERANCE << "\t\t\t# particles with mass and chemical potential" << endl
		<< tentab << "# (for each FO-cell) difference less than this value" << endl
		<< tentab << "# will be considered to be identical (b/c Cooper-Frye)" << endl << endl;

	output_HCparams
		<< "# Pair momenta" << endl
		<< "nKT = " << nKT << "\t\t\t\t\t\t\t\t# size of pair momentum KT grid" << endl
		<< "nKphi = " << nKphi << "\t\t\t\t\t\t\t\t# size of pair momentum Kphi grid" << endl
		<< "KTmin = " << KTmin << "\t\t\t\t\t\t\t# minimum of pair momentum KT grid" << endl
		<< "KTmax = " << KTmax << "\t\t\t\t\t\t\t# maximum of pair momentum KT grid" << endl << endl;

	output_HCparams
		<< "##################" << endl
		<< "## SVWR parameters" << endl << endl

		<< "# Sub-model parameters" << endl
		<< "SV_resonanceThreshold = " << SV_resonanceThreshold << "\t\t\t# (approximate) fraction of all resonance contributions to do" << endl
		<< tentab << "# code sums over all resonances in (approximately) decreasing" << endl
		<< tentab << "# order until this threshhold is reached," << endl
		<< tentab << "# only skipping resonances with lifetimes > max_lifetime" << endl << endl;

	output_HCparams
		<< "# Single-particle momenta" << endl
		<< "SV_npT = " << SV_npT << "\t\t\t\t\t\t\t\t# size of pT grid for interpolation of (SV-weighted)" << endl
		<< tentab << "# single-particle spectra" << endl << endl;

	output_HCparams
		<< "SV_npphi = " << SV_npphi << "\t\t\t\t\t\t\t# size of pphi grid for interpolation of (SV-weighted)" << endl
		<< tentab << "# single-particle spectra" << endl << endl;


	output_HCparams
		<< "##################" << endl
		<< "## CFWR parameters" << endl << endl

		<< "# Sub-model options" << endl;

	output_HCparams
		<< "calculate_CF_mode = " << CALCULATE_CF_MODE << "\t\t\t\t\t# 0 (default) - compute resonance decays, correlation function, and fit" << endl
		<< tentab << "# 1 - read in resonance decays, compute correlation function, and fit" << endl
		<< tentab << "# 2 - read in correlation function and fit" << endl << endl;

	output_HCparams
		<< "fit_with_projected_cfvals = " << FIT_WITH_PROJECTED_CFVALS << "\t\t\t# enable linear extrapolation of sum over chosen" << endl
		<< tentab << "# resonances to estimate sum over remaining resonances" << endl << endl;

	output_HCparams
		<< "flesh_out_cf = " << FLESH_OUT_CF << "\t\t\t\t\t\t# refine grid via interpolation before fitting to enable" << endl
		<< tentab << "# more reliable fit results" << endl << endl;

	output_HCparams
		<< "# Sub-model parameters" << endl
		<< "CF_resonanceThreshold = " << CF_resonanceThreshold << "\t\t\t# (approximate) fraction of all resonance contributions to do" << endl
		<< tentab << "# code sums over all resonances in (approximately) decreasing" << endl
		<< tentab << "# order until this threshhold is reached," << endl
		<< tentab << "# only skipping resonances with lifetimes > max_lifetime" << endl << endl;

	output_HCparams
		<< "# Single-particle momenta" << endl
		<< "CF_npT = " << n_pT_pts << "\t\t\t\t\t\t\t\t# size of pT grid for interpolation of" << endl
		<< tentab << "# (Fourier-transformed) single-particle spectra" << endl << endl;

	output_HCparams
		<< "CF_npphi = " << n_pphi_pts << "\t\t\t\t\t\t\t# size of pphi grid for interpolation of" << endl
		<< tentab << "# (Fourier-transformed) single-particle spectra" << endl << endl;

	output_HCparams
		<< "CF_npY = " << n_pY_pts << "\t\t\t\t\t\t\t\t# size of pY grid for interpolation of" << endl
		<< tentab << "# (Fourier-transformed) single-particle spectra" << endl << endl;

	output_HCparams
		<< "# Relative momenta" << endl << endl;

	output_HCparams
		<< "qtnpts = " << qtnpts << endl
		<< "qxnpts = " << qxnpts << endl
		<< "qynpts = " << qynpts << endl
		<< "qznpts = " << qznpts << endl
		<< "delta_qx = " << delta_qx << endl
		<< "delta_qy = " << delta_qy << endl
		<< "delta_qz = " << delta_qz << endl;

	output_HCparams
		<< "## End of file" << endl;

	output_HCparams.close();

	return;

}
*/


//End of file
