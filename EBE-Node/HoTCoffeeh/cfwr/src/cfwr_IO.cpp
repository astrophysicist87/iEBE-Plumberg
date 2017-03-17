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
#include "Arsenal.h"
#include "gauss_quadrature.h"

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
void CorrelationFunction::Dump_spectra_array(string output_filename, double *** array_to_dump)
{
	ostringstream filename_stream;
	filename_stream << path << "/" << output_filename;
	ofstream out(filename_stream.str().c_str());

	for (int ir = 0; ir < Nparticle; ++ir)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	{
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
			out << scientific << setprecision(8) << setw(12) << array_to_dump[ir][ipT][ipphi] << "   ";
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
		for (int ipt = 0; ipt < n_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			flat_R2s[iptipphi] = R2_side_GF[ipt][ipphi];
			flat_R2o[iptipphi] = R2_out_GF[ipt][ipphi];
			flat_R2l[iptipphi] = R2_long_GF[ipt][ipphi];
			flat_R2os[iptipphi] = R2_outside_GF[ipt][ipphi];
			flat_R2sl[iptipphi] = R2_sidelong_GF[ipt][ipphi];
			flat_R2ol[iptipphi] = R2_outlong_GF[ipt][ipphi];
			iptipphi++;
		}
	}
	else if (mode == 1)
	{
		for (int ipt = 0; ipt < n_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			flat_R2s[iptipphi] = R2_side_QM[ipt][ipphi];
			flat_R2o[iptipphi] = R2_out_QM[ipt][ipphi];
			flat_R2l[iptipphi] = R2_long_QM[ipt][ipphi];
			flat_R2os[iptipphi] = R2_outside_QM[ipt][ipphi];
			flat_R2sl[iptipphi] = R2_sidelong_QM[ipt][ipphi];
			flat_R2ol[iptipphi] = R2_outlong_QM[ipt][ipphi];
			iptipphi++;
		}
	}

	//output R2ij on original pT-pphi grid
	if (mode == 0)
	{
		for (int ipt = 0; ipt < n_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			outputHBT_g0 << SP_pT[ipt] << "   " << SP_pphi[ipphi]
				<< "   " << R2_side_GF[ipt][ipphi] << "   " << R2_out_GF[ipt][ipphi]
				<< "   " << R2_outside_GF[ipt][ipphi] << "   " << R2_long_GF[ipt][ipphi]
				<< "   " << R2_sidelong_GF[ipt][ipphi] << "   " << R2_outlong_GF[ipt][ipphi] << endl;
		}
	}
	else if (mode == 1)
	{
		for (int ipt = 0; ipt < n_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			outputHBT_g0 << SP_pT[ipt] << "   " << SP_pphi[ipphi]
				<< "   " << R2_side_QM[ipt][ipphi] << "   " << R2_out_QM[ipt][ipphi]
				<< "   " << R2_outside_QM[ipt][ipphi] << "   " << R2_long_QM[ipt][ipphi]
				<< "   " << R2_sidelong_QM[ipt][ipphi] << "   " << R2_outlong_QM[ipt][ipphi] << endl;
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
					<< "  " << R2_side_GF_C[iKT][Morder] << "   " << R2_side_GF_S[iKT][Morder] << "  " << R2_out_GF_C[iKT][Morder] << "  " << R2_out_GF_S[iKT][Morder]
					<< "  " << R2_outside_GF_C[iKT][Morder] << "   " << R2_outside_GF_S[iKT][Morder] << "  " << R2_long_GF_C[iKT][Morder] << "  " << R2_long_GF_S[iKT][Morder]
					<< "  " << R2_sidelong_GF_C[iKT][Morder] << "   " << R2_sidelong_GF_S[iKT][Morder] << "  " << R2_outlong_GF_C[iKT][Morder] << "  " << R2_outlong_GF_S[iKT][Morder] << endl;
			}
		}
		else if (mode == 1)
		{
			for (int Morder = 0; Morder < n_order; Morder++)
			{
				outputHBTcfs << SP_pT[iKT] << "  " << Morder
					<< "  " << R2_side_QM_C[iKT][Morder] << "   " << R2_side_QM_S[iKT][Morder] << "  " << R2_out_QM_C[iKT][Morder] << "  " << R2_out_QM_S[iKT][Morder]
					<< "  " << R2_outside_QM_C[iKT][Morder] << "   " << R2_outside_QM_S[iKT][Morder] << "  " << R2_long_QM_C[iKT][Morder] << "  " << R2_long_QM_S[iKT][Morder]
					<< "  " << R2_sidelong_QM_C[iKT][Morder] << "   " << R2_sidelong_QM_S[iKT][Morder] << "  " << R2_outlong_QM_C[iKT][Morder] << "  " << R2_outlong_QM_S[iKT][Morder] << endl;
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
	
	for (int ipt = 0; ipt < n_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		outputlambdas << SP_pT[ipt] << "   " << SP_pphi[ipphi] << "   " << lambda_Correl[ipt][ipphi] << endl;

	outputlambdas.close();
	return;
}

void CorrelationFunction::Output_correlationfunction()
{
	ostringstream oCorrFunc_stream;
	string temp_particle_name = particle_name;
	replace_parentheses(temp_particle_name);

	string CF_proj_string = "";
	if (!FIT_WITH_PROJECTED_CFVALS)
		CF_proj_string = "unprojected_";

	oCorrFunc_stream << path << "/correlfunct3D_" << CF_proj_string << temp_particle_name << ".dat";
	ofstream oCorrFunc;
	oCorrFunc.open(oCorrFunc_stream.str().c_str());

	for (int ipt = 0; ipt < n_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	{
		double ckp = cos_SP_pphi[ipphi], skp = sin_SP_pphi[ipphi];
		oCorrFunc << scientific << setprecision(8) << setw(12)
			<< SP_pT[ipt] << "   " << SP_pphi[ipphi] << "   " << qx_pts[iqx] << "   "
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

void CorrelationFunction::Output_fleshed_out_correlationfunction(int ipt, int ipphi)
{
	ostringstream oCorrFunc_stream;
	string temp_particle_name = particle_name;
	replace_parentheses(temp_particle_name);
	oCorrFunc_stream << path << "/correlfunct3D" << "_" << temp_particle_name << "_fleshed_out.dat";
	ofstream oCorrFunc;
	if (ipt==0 && ipphi==0)
		oCorrFunc.open(oCorrFunc_stream.str().c_str());
	else
		oCorrFunc.open(oCorrFunc_stream.str().c_str(), ios::app);

	for (int iqx = 0; iqx < new_nqpts; ++iqx)
	for (int iqy = 0; iqy < new_nqpts; ++iqy)
	for (int iqz = 0; iqz < new_nqpts; ++iqz)
	{
		double ckp = cos_SP_pphi[ipphi], skp = sin_SP_pphi[ipphi];
		oCorrFunc << scientific << setprecision(7) << setw(15)
			<< SP_pT[ipt] << "   " << SP_pphi[ipphi] << "   " << qx_fleshed_out_pts[iqx] << "   "
			<< qy_fleshed_out_pts[iqy] << "   " << qz_fleshed_out_pts[iqz] << "   "
			<< qx_fleshed_out_pts[iqx] * ckp + qy_fleshed_out_pts[iqy] * skp << "   "
			<< -qx_fleshed_out_pts[iqx] * skp + qy_fleshed_out_pts[iqy] * ckp << "   "
			<< qz_fleshed_out_pts[iqz] << "   "
			<< fleshed_out_thermal[iqx][iqy][iqz] << "   " << fleshed_out_crossterm[iqx][iqy][iqz] << "   " << fleshed_out_resonances[iqx][iqy][iqz] << "   " << fleshed_out_CF[iqx][iqy][iqz] << endl;
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
		for (int ipt = 0; ipt < n_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
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
	}
	else if (mode == 1)
	{
		for (int ipt = 0; ipt < n_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			inputHBT >> dummy;	//pt value
			inputHBT >> dummy;	//pphi value
			inputHBT >> R2_side_QM[ipt][ipphi];
			inputHBT >> R2_out_QM[ipt][ipphi];
			inputHBT >> R2_outside_QM[ipt][ipphi];
			inputHBT >> R2_long_QM[ipt][ipphi];
			inputHBT >> R2_sidelong_QM[ipt][ipphi];
			inputHBT >> R2_outlong_QM[ipt][ipphi];
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
	filename_stream_target_dN_dypTdpTdphi << path << "/total_" << local_name << "_dN_dypTdpTdphi_" << no_df_stem << ".dat";
	ofstream output_target_dN_dypTdpTdphi(filename_stream_target_dN_dypTdpTdphi.str().c_str());

	for(int ipphi = 0; ipphi < n_pphi_pts; ipphi++)
	{
		for(int ipt = 0; ipt < n_pT_pts; ipt++)
		{
			double fs = spectra[target_particle_id][ipt][ipphi];
			double ts = thermal_spectra[target_particle_id][ipt][ipphi];
			double result = (fs - ts) / fraction_of_resonances + ts;
			output_target_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12) << spectra[target_particle_id][ipt][ipphi] << "   ";
		}
		output_target_dN_dypTdpTdphi << endl;
	}

	output_target_dN_dypTdpTdphi.close();

	return;
}

void CorrelationFunction::Output_total_target_eiqx_dN_dypTdpTdphi(double current_fraction /*==-1.0*/)
{
	string local_name = all_particles[target_particle_id].name;
	string current_fraction_string = (current_fraction >= 0.0) ? "_" + patch::to_string(current_fraction) : "";
	replace_parentheses(local_name);
	ostringstream filename_stream_target_dN_dypTdpTdphi;
	filename_stream_target_dN_dypTdpTdphi << path << "/total_" << local_name << current_fraction_string << "_eiqx_dN_dypTdpTdphi" << no_df_stem << ".dat";
	ofstream output_target_dN_dypTdpTdphi(filename_stream_target_dN_dypTdpTdphi.str().c_str());

	//int HDFloadTargetSuccess = Get_resonance_from_HDF_array(target_particle_id, current_dN_dypTdpTdphi_moments);
	Set_full_target_moments();

	// addresses NaN issue in sin component when all q^{\mu} == 0
	if (qtnpts%2==1 && qxnpts%2==1 && qynpts%2==1 && qznpts%2==1)
	{	//if all q-ranges are odd and centered on q=0 ==> q=0 is included!
		int iqt0 = (qtnpts-1)/2;
		int iqx0 = (qxnpts-1)/2;
		int iqy0 = (qynpts-1)/2;
		int iqz0 = (qznpts-1)/2;
		for (int ipt = 0; ipt < n_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		{
			current_dN_dypTdpTdphi_moments[indexer(ipt,ipphi,iqt0,iqx0,iqy0,iqz0,1)] = 0.0;
			thermal_target_dN_dypTdpTdphi_moments[indexer(ipt,ipphi,iqt0,iqx0,iqy0,iqz0,1)] = 0.0;
		}
	}

	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	for (int ipt = 0; ipt < n_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	{
		//first, get CF and projected CF
		double CF = get_CF(ipt, ipphi, iqt, iqx, iqy, iqz, false);				//false means don't return projected value

		//!!!!!!!!!!!!should get projected_CF AFTER regulating CF...!!!!!!!!!!!!
		double projected_CF = get_CF(ipt, ipphi, iqt, iqx, iqy, iqz, true && !thermal_pions_only);	//true means do return projected value

		double nonFTd_spectra = spectra[target_particle_id][ipt][ipphi];
		double cos_transf_spectra = current_dN_dypTdpTdphi_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,0)];
		double sin_transf_spectra = current_dN_dypTdpTdphi_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,1)];

		output_target_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12)
			<< qt_pts[iqt] << "   " << qx_pts[iqx] << "   " << qy_pts[iqy] << "   " << qz_pts[iqz] << "   "
			<< SP_pT[ipt] << "   " << SP_pphi[ipphi] << "   "
			<< nonFTd_spectra << "   "																								//non-thermal + thermal
			<< cos_transf_spectra << "   "																							//non-thermal + thermal (cos)
			<< sin_transf_spectra << "   "																							//non-thermal + thermal (sin)
			<< thermal_spectra[target_particle_id][ipt][ipphi] << "   "																//thermal only
			<< thermal_target_dN_dypTdpTdphi_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,0)] << "   "									//thermal only (cos)
			<< thermal_target_dN_dypTdpTdphi_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,1)] << "   "									//thermal only (sin)
			<< nonFTd_spectra - thermal_spectra[target_particle_id][ipt][ipphi] << "   "											//non-thermal only
			<< cos_transf_spectra - thermal_target_dN_dypTdpTdphi_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,0)] << "   "				//non-thermal only (cos)
			<< sin_transf_spectra - thermal_target_dN_dypTdpTdphi_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,1)] << "   "				//non-thermal only (sin)
			<< CF << "   " << projected_CF << endl;
	}

	output_target_dN_dypTdpTdphi.close();

	return;
}

void CorrelationFunction::Output_total_eiqx_dN_dypTdpTdphi(int local_pid)
{
	string local_name = all_particles[local_pid].name;
	replace_parentheses(local_name);
	ostringstream filename_stream_dN_dypTdpTdphi;
	filename_stream_dN_dypTdpTdphi << path << "/total_" << local_name << "_eiqx_dN_dypTdpTdphi_" << no_df_stem << ".dat";
	ofstream output_dN_dypTdpTdphi(filename_stream_dN_dypTdpTdphi.str().c_str());

	int HDFOpenSuccess = Open_resonance_HDF_array("resonance_spectra.h5");
	int HDFloadTargetSuccess = Get_resonance_from_HDF_array(local_pid, current_dN_dypTdpTdphi_moments);
	int HDFCloseSuccess = Close_resonance_HDF_array();

	// addresses NaN issue in sin component when all q^{\mu} == 0
	if (qtnpts%2==1 && qxnpts%2==1 && qynpts%2==1 && qznpts%2==1)
	{	//if all q-ranges are odd and centered on q=0 ==> q=0 is included!
		int iqt0 = (qtnpts-1)/2;
		int iqx0 = (qxnpts-1)/2;
		int iqy0 = (qynpts-1)/2;
		int iqz0 = (qznpts-1)/2;
		for (int ipt = 0; ipt < n_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
			current_dN_dypTdpTdphi_moments[indexer(ipt,ipphi,iqt0,iqx0,iqy0,iqz0,1)] = 0.0;
	}

	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	for (int ipt = 0; ipt < n_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	{
		double nonFTd_spectra = spectra[local_pid][ipt][ipphi];
		double cos_transf_spectra = current_dN_dypTdpTdphi_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,0)];
		double sin_transf_spectra = current_dN_dypTdpTdphi_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,1)];

		output_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12)
			<< qt_pts[iqt] << "   " << qx_pts[iqx] << "   " << qy_pts[iqy] << "   " << qz_pts[iqz] << "   "
			<< SP_pT[ipt] << "   " << SP_pphi[ipphi] << "   "
			<< nonFTd_spectra << "   "																								//non-thermal + thermal
			<< cos_transf_spectra << "   "																							//non-thermal + thermal (cos)
			<< sin_transf_spectra << endl;
	}

	output_dN_dypTdpTdphi.close();

	return;
}

void CorrelationFunction::Readin_total_target_eiqx_dN_dypTdpTdphi()
{
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
	for (int ipt = 0; ipt < n_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	{
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> spectra[target_particle_id][ipt][ipphi];
		input_target_dN_dypTdpTdphi >> current_dN_dypTdpTdphi_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,0)];
		input_target_dN_dypTdpTdphi >> current_dN_dypTdpTdphi_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,1)];
		input_target_dN_dypTdpTdphi >> thermal_spectra[target_particle_id][ipt][ipphi];
		input_target_dN_dypTdpTdphi >> thermal_target_dN_dypTdpTdphi_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,0)];
		input_target_dN_dypTdpTdphi >> thermal_target_dN_dypTdpTdphi_moments[indexer(ipt,ipphi,iqt,iqx,iqy,iqz,1)];
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
	for (int ipt = 0; ipt < n_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	{
		iCorrFunc >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
					>> thermalCFvals[ipt][ipphi][iqx][iqy][iqz] >> crosstermCFvals[ipt][ipphi][iqx][iqy][iqz] >> resonancesCFvals[ipt][ipphi][iqx][iqy][iqz] >> CFvals[ipt][ipphi][iqx][iqy][iqz];
		//*global_out_stream_ptr << "Read in CFvals[" << ipt << "][" << ipphi << "][" << iqx << "][" << iqy << "][" << iqz << "] = " << CFvals[ipt][ipphi][iqx][iqy][iqz] << endl;
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

//End of file
