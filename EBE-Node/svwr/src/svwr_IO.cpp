#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>

#include<gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>            // gsl random number generators
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting

#include "svwr.h"
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

void SourceVariances::Output_results(int folderindex)
{
	ostringstream filename_stream_HBT;
	filename_stream_HBT << global_path << "/HBTradii_SVWR_ev" << folderindex << no_df_stem << ".dat";
	ofstream outputHBT;
	outputHBT.open(filename_stream_HBT.str().c_str());
	ostringstream filename_stream_HBTcfs;
	filename_stream_HBTcfs << global_path << "/HBTradii_SVWR_cfs_ev" << folderindex << no_df_stem << ".dat";
	ofstream outputHBTcoeffs(filename_stream_HBTcfs.str().c_str());
	ostringstream filename_stream_S;
	filename_stream_S << global_path << "/Sourcefunction_variances_WR" << no_df_stem << ".dat";
	ofstream output_Svars(filename_stream_S.str().c_str());

	for(int iKT = 0; iKT < n_localp_T; iKT++)
	{
		for(int Morder=0; Morder<n_order; Morder++)
		{
			outputHBTcoeffs << folderindex << "  " << K_T[iKT] << "  " << Morder
				<< "  " << R2_side_C[iKT][Morder] << "   " << R2_side_S[iKT][Morder] << "  " << R2_out_C[iKT][Morder] << "  " << R2_out_S[iKT][Morder]
				<< "  " << R2_outside_C[iKT][Morder] << "   " << R2_outside_S[iKT][Morder] << "  " << R2_long_C[iKT][Morder] << "  " << R2_long_S[iKT][Morder]
				<< "  " << R2_sidelong_C[iKT][Morder] << "   " << R2_sidelong_S[iKT][Morder] << "  " << R2_outlong_C[iKT][Morder] << "  " << R2_outlong_S[iKT][Morder] << endl;
		}
		for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
		{
			outputHBT << folderindex << "  " << K_T[iKT] << "  " << K_phi[iKphi]
				<< "  " << R2_side[iKT][iKphi] << "  " << R2_out[iKT][iKphi]
				<< "  " << R2_outside[iKT][iKphi] << "  " << R2_long[iKT][iKphi]
				<< "  " << R2_sidelong[iKT][iKphi] << "  " << R2_outlong[iKT][iKphi] << endl;
		}
	}
	for(int ipt = 0; ipt < n_interp_pT_pts; ipt++)
	for(int ipphi = 0; ipphi < n_interp_pphi_pts; ipphi++)
		output_Svars << setprecision(8) << setw(15) 
			<< SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi] << "   " << S_func[ipt][ipphi] << "   "
			<< xs_S[ipt][ipphi] << "   " << xo_S[ipt][ipphi] << "   " << xl_S[ipt][ipphi] << "   "
			<< t_S[ipt][ipphi]  << "   " << xs_t_S[ipt][ipphi] << "   "
			<< xo_t_S[ipt][ipphi] << "   " << xl_t_S[ipt][ipphi] << "   "
			<< xo_xs_S[ipt][ipphi] << "   " << xl_xs_S[ipt][ipphi] << "   "
			<< xo_xl_S[ipt][ipphi] << "   " << xs2_S[ipt][ipphi] << "   " << xo2_S[ipt][ipphi] << "   "
			<< xl2_S[ipt][ipphi] << "   " << t2_S[ipt][ipphi] << endl;

	outputHBT.close();
	output_Svars.close();

	return;
}

void SourceVariances::Readin_results(int folderindex)
{
	double dummy;
	ostringstream filename_stream_HBT;
	filename_stream_HBT << global_path << "/HBTradii_SVWR_ev" << folderindex << no_df_stem << ".dat";
	ifstream inputHBT(filename_stream_HBT.str().c_str());
	ostringstream filename_stream_S;
	filename_stream_S << global_path << "/Sourcefunction_variances_WR" << no_df_stem << ".dat";
	ifstream input_Svars(filename_stream_S.str().c_str());

	for(int iKT = 0; iKT < n_localp_T; iKT++)
	{
		for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
		{
			inputHBT >> dummy;
			inputHBT >> dummy;
			inputHBT >> dummy;
			inputHBT >> R2_side[iKT][iKphi];
			inputHBT >> R2_out[iKT][iKphi];
			inputHBT >> R2_outside[iKT][iKphi];
			inputHBT >> R2_long[iKT][iKphi];
			inputHBT >> R2_sidelong[iKT][iKphi];
			inputHBT >> R2_outlong[iKT][iKphi];
		}
	}
	for(int ipt = 0; ipt < n_interp_pT_pts; ipt++)
	for(int ipphi = 0; ipphi < n_interp_pphi_pts; ipphi++)
	{
         	input_Svars >> dummy;
        	input_Svars >> dummy;
        	input_Svars >> S_func[ipt][ipphi];
        	input_Svars >> xs_S[ipt][ipphi];
        	input_Svars >> xo_S[ipt][ipphi];
        	input_Svars >> xl_S[ipt][ipphi];
        	input_Svars >> t_S[ipt][ipphi];
        	input_Svars >> xs_t_S[ipt][ipphi];
        	input_Svars >> xo_t_S[ipt][ipphi];
        	input_Svars >> xl_t_S[ipt][ipphi];
        	input_Svars >> xo_xs_S[ipt][ipphi];
        	input_Svars >> xl_xs_S[ipt][ipphi];
        	input_Svars >> xo_xl_S[ipt][ipphi];
        	input_Svars >> xs2_S[ipt][ipphi];
        	input_Svars >> xo2_S[ipt][ipphi];
        	input_Svars >> xl2_S[ipt][ipphi];
        	input_Svars >> t2_S[ipt][ipphi];
	}

	inputHBT.close();
	input_Svars.close();

	return;
}

/*void SourceVariances::Readin_ev_plane_psi(int folderindex)
{
	ostringstream filename_stream_planepsi;
	//filename_stream_planepsi << path << folderindex << "/plane_psi_ev" << folderindex << ".dat";
	filename_stream_planepsi << global_path << "/plane_psi_ev" << folderindex << no_df_stem << ".dat";
	ifstream inputplanepsi(filename_stream_planepsi.str().c_str());

	inputplanepsi >> global_plane_psi;

	inputplanepsi.close();

	return;
}

void SourceVariances::Output_ev_plane_psi(int folderindex)
{
	ostringstream filename_stream_planepsi;
	//filename_stream_planepsi << path << folderindex << "/plane_psi_ev" << folderindex << ".dat";
	filename_stream_planepsi << global_path << "/plane_psi_ev" << folderindex << no_df_stem << ".dat";
	ofstream outputplanepsi(filename_stream_planepsi.str().c_str());

	outputplanepsi << global_plane_psi << endl;

	outputplanepsi.close();

	return;
}

void SourceVariances::Output_ev_plane_psis(int folderindex)
{
	ostringstream filename_stream_planepsis;
	//filename_stream_planepsis << path << folderindex << "/plane_psis_ev" << folderindex << ".dat";
	filename_stream_planepsis << global_path << "/plane_psis_ev" << folderindex << no_df_stem << ".dat";
	ofstream outputplanepsis(filename_stream_planepsis.str().c_str());

	for (int i = 0; i < n_order; i++)
		outputplanepsis << i << "   " << plane_angle[i] << endl;

	outputplanepsis.close();

	return;
}*/

void SourceVariances::Output_dN_dypTdpTdphi(int folderindex)
{
	ostringstream filename_stream_dN_dypTdpTdphi;
	filename_stream_dN_dypTdpTdphi << global_path << "/dN_dypTdpTdphi_ev" << folderindex << no_df_stem << ".dat";
	ofstream output_dN_dypTdpTdphi(filename_stream_dN_dypTdpTdphi.str().c_str());

	for(int iphi=0; iphi<n_SP_pphi; iphi++)
	for(int ipt=0; ipt<n_SP_pT; ipt++)
		output_dN_dypTdpTdphi << SP_pT[ipt] << "   " << SP_pphi[iphi] << "   " << dN_dypTdpTdphi[ipt][iphi] << endl;

	output_dN_dypTdpTdphi.close();

	return;
}

void SourceVariances::Output_dN_dypTdpT(int folderindex)
{
	ostringstream filename_stream_dN_dypTdpT;
	filename_stream_dN_dypTdpT << global_path << "/dN_dypTdpT_ev" << folderindex << no_df_stem << ".dat";
	ofstream output_dN_dypTdpT(filename_stream_dN_dypTdpT.str().c_str());

	for(int ipt=0; ipt<n_SP_pT; ipt++)
		output_dN_dypTdpT << SP_pT[ipt] << "   " << dN_dypTdpT[ipt] << endl;

	output_dN_dypTdpT.close();

	return;
}

void SourceVariances::Output_all_dN_dypTdpTdphi(int folderindex)
{
	for(int wfi = 0; wfi < n_weighting_functions; wfi++)
	{
		ostringstream filename_stream_all_dN_dypTdpTdphi;
		filename_stream_all_dN_dypTdpTdphi << global_path << "/all_res_dN_dypTdpTdphi_mom_"
						<< setfill('0') << setw(2) << wfi << "_ev" << folderindex << no_df_stem << ".dat";
		ofstream output_all_dN_dypTdpTdphi(filename_stream_all_dN_dypTdpTdphi.str().c_str());
		for(int ii = 0; ii < Nparticle; ii++)
		for(int iphi = 0; iphi < n_interp_pphi_pts; iphi++)
		{
			for(int ipt = 0; ipt < n_interp_pT_pts; ipt++)
				output_all_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12) << dN_dypTdpTdphi_moments[ii][wfi][ipt][iphi] << "   ";
			output_all_dN_dypTdpTdphi << endl;
		}
		output_all_dN_dypTdpTdphi.close();
	}

	return;
}

void SourceVariances::Output_total_target_dN_dypTdpTdphi(int folderindex)
{
	string local_name = all_particles[target_particle_id].name;
	replace_parentheses(local_name);

	for(int wfi = 0; wfi < n_weighting_functions; wfi++)
	{
		ostringstream filename_stream_target_dN_dypTdpTdphi;
		filename_stream_target_dN_dypTdpTdphi << global_path << "/total_" << local_name << "_dN_dypTdpTdphi_mom_"
								<< setfill('0') << setw(2) << wfi << "_ev" << folderindex << no_df_stem << ".dat";
		ofstream output_target_dN_dypTdpTdphi(filename_stream_target_dN_dypTdpTdphi.str().c_str());
	
		for(int iphi = 0; iphi < n_interp_pphi_pts; iphi++)
		{
			for(int ipt = 0; ipt < n_interp_pT_pts; ipt++)
				output_target_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12) << dN_dypTdpTdphi_moments[target_particle_id][wfi][ipt][iphi] << "   ";
			output_target_dN_dypTdpTdphi << endl;
		}
	
		output_target_dN_dypTdpTdphi.close();
	}

	return;
}

void SourceVariances::Output_chosen_resonances()
{
	ostringstream filename_stream_crf;
	filename_stream_crf << global_path << "/chosen_resonances.dat";
	ofstream output_crf(filename_stream_crf.str().c_str());

	output_crf << particle_monval << endl;
	for (int icr = 0; icr < (int)chosen_resonances.size(); icr++)
		output_crf << all_particles[chosen_resonances[icr]].monval << endl;

	output_crf.close();

	return;
}

void SourceVariances::Read_in_all_dN_dypTdpTdphi(int folderindex)
{
	for(int wfi = 0; wfi < n_weighting_functions; wfi++)
	{
		ostringstream filename_stream_all_dN_dypTdpTdphi;
		filename_stream_all_dN_dypTdpTdphi << global_path << "/all_res_dN_dypTdpTdphi_mom_"
								<< setfill('0') << setw(2) << wfi << "_ev" << folderindex << no_df_stem << ".dat";
		ifstream input_all_dN_dypTdpTdphi(filename_stream_all_dN_dypTdpTdphi.str().c_str());
	
		int local_filelength = get_filelength(filename_stream_all_dN_dypTdpTdphi.str().c_str());
		int local_filewidth = get_filewidth(filename_stream_all_dN_dypTdpTdphi.str().c_str());
		if (VERBOSE > 0) *global_out_stream_ptr << "Read_in_all_dN_dypTdpTdphi(): nrows = "
							<< local_filelength << " and ncols = " << local_filewidth << endl;
		if ((Nparticle * n_interp_pphi_pts != local_filelength) || (n_interp_pT_pts != local_filewidth))
		{
			cerr << "Read_in_all_dN_dypTdpTdphi(): Mismatch in dimensions in file "
				<< "all_res_dN_dypTdpTdphi_mom_" << setfill('0') << setw(2) << wfi
				<< "_ev" << folderindex << no_df_stem << ".dat!" << endl;
			exit(1);
		}
	
		for(int ii = 0; ii < Nparticle; ii++)
		for(int iphi = 0; iphi < n_interp_pphi_pts; iphi++)
		for(int ipt = 0; ipt < n_interp_pT_pts; ipt++)
		{
			input_all_dN_dypTdpTdphi >> dN_dypTdpTdphi_moments[ii][wfi][ipt][iphi];
			if (abs(dN_dypTdpTdphi_moments[ii][wfi][ipt][iphi]) > 1.e-100)
			{
				ln_dN_dypTdpTdphi_moments[ii][wfi][ipt][iphi] = log(abs(dN_dypTdpTdphi_moments[ii][wfi][ipt][iphi]));
				sign_of_dN_dypTdpTdphi_moments[ii][wfi][ipt][iphi] = sgn(dN_dypTdpTdphi_moments[ii][wfi][ipt][iphi]);
			}
			//cout << "Read in (pT, pphi, EdNdp3 ST moms) = " << SPinterp_pT[ipt] << "   " << SPinterp_pphi[iphi] << "   "
			//	<< scientific << setprecision(8) << setw(12) << dN_dypTdpTdphi_moments[ii][wfi][ipt][iphi] << endl;
		}
	
		input_all_dN_dypTdpTdphi.close();
		if (VERBOSE > 0) *global_out_stream_ptr << "Successfully read in " << filename_stream_all_dN_dypTdpTdphi.str().c_str() << endl;
	}


	ostringstream filename_stream_pTpts;
	filename_stream_pTpts << global_path << "/pT_gauss_table.dat";
	ifstream input_pTpts(filename_stream_pTpts.str().c_str());
	ostringstream filename_stream_pphipts;
	filename_stream_pphipts << global_path << "/phi_gauss_table.dat";
	ifstream input_pphipts(filename_stream_pphipts.str().c_str());

	double * dummy_pT_wts = new double [n_interp_pT_pts];
	double * dummy_pphi_wts = new double [n_interp_pphi_pts];

	for(int ipt = 0; ipt < n_interp_pT_pts; ipt++)
		input_pTpts >> SPinterp_pT[ipt] >> dummy_pT_wts[ipt];

	for(int ipphi = 0; ipphi < n_interp_pphi_pts; ipphi++)
		input_pphipts >> SPinterp_pphi[ipphi] >> dummy_pphi_wts[ipphi];

	if (VERBOSE > 0) *global_out_stream_ptr << "Read_in_all_dN_dypTdpTdphi(): read in pT and pphi points!" << endl;

	input_pTpts.close();
	input_pphipts.close();

	return;
}

//End of file
