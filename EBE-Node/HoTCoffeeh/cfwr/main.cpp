#include <omp.h>
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

#include "src/Stopwatch.h"
#include "src/parameters.h"
#include "src/readindata.h"
#include "src/cfwr.h"
#include "src/generate_processing_record.h"
#include "src/lib.h"
#include "src/ParameterReader.h"
#include "main.h"

using namespace std;

int main(int argc, char *argv[])
{
   cout << endl
        << "            Correlation functions with resonances            " << endl
        << endl
        << "  Ver 1.0   ----- Christopher Plumberg, 07/2016   " << endl;
   cout << endl << "**********************************************************" << endl;
   display_logo(2); // Hail to the king~
   cout << endl << "**********************************************************" << endl << endl;
   
   // Read-in parameters
   ParameterReader * paraRdr = new ParameterReader;
   paraRdr->readFromFile("parameters.dat");
   paraRdr->readFromArguments(argc, argv);
   paraRdr->echo();

////////////////////////////////////////////////////////////////////
// below are Chris' set-up for computing the correlation function
////////////////////////////////////////////////////////////////////

	Stopwatch sw;
	Stopwatch sw_total;
	sw_total.Start();
	sw.Start();

	bool flag_nans = false;

	//flags NaNs when they crop up
	if (flag_nans)
		feenableexcept(FE_INVALID | FE_OVERFLOW);

	bool generatedcorrfuncs = false;
    string workingDirectory = "./results";

	initialize_PRfile(paraRdr, workingDirectory);

	ostringstream filename_stream;
	filename_stream << workingDirectory << "/Processing_record.txt";
	ofstream output(filename_stream.str().c_str(), ios::app);

	output << "/**********Processing output**********/" << endl;
	output << "entering folder: " << workingDirectory << endl;
	if (flag_nans)
		output << "***NOTE: flagging NaNs!***" << endl;

	//load freeze out and particle information
	int FO_length = 0;
	int particle_idx = 1;  //for pion+

	ostringstream decdatfile;
	output << "Loading the decoupling data...." << endl;
	decdatfile << workingDirectory << "/decdat2.dat";
	output << decdatfile.str() << endl;
	FO_length=get_filelength(decdatfile.str().c_str());
	output << "Total number of freeze out fluid cell: " <<  FO_length << endl;

	//read the data arrays for the decoupling information
	FO_surf* FOsurf_ptr = new FO_surf[FO_length];
	bool including_bulkpi = paraRdr->getVal("include_bulk_pi");
	read_decdat(FO_length, FOsurf_ptr, workingDirectory, including_bulkpi);
   
	//read the positions of the freeze out surface
	read_surfdat(FO_length, FOsurf_ptr, workingDirectory);
   
	//read the chemical potential on the freeze out surface
	int N_stableparticle;
	ifstream particletable("EOS/EOS_particletable.dat");
	particletable >> N_stableparticle;
	output << "Treating " << N_stableparticle << " stable particles!" << endl;
	double ** particle_mu = new double * [N_stableparticle];
	for (int i = 0; i < N_stableparticle; i++)
		particle_mu[i] = new double [FO_length];
	for (int i = 0; i < N_stableparticle; i++)
	for (int j = 0; j < FO_length; j++)
		particle_mu[i][j] = 0.0;
	if (N_stableparticle > 0)
		read_decdat_mu(FO_length, N_stableparticle, particle_mu, workingDirectory);

	//read particle resonance decay table
	particle_info *particle = new particle_info [Maxparticle];
	int Nparticle = read_resonance(particle);
	output << "read in total " << Nparticle << " particles!" << endl;
	output << "Calculating " << particle[particle_idx].name << endl;
	if (N_stableparticle > 0)
	{
		output << "EOS is partially chemical equilibrium " << endl;
		calculate_particle_mu(7, Nparticle, FOsurf_ptr, FO_length, particle, particle_mu);
	}
	else
	{
		output << "EOS is chemical equilibrium. " << endl;
		for(int j=0; j<Nparticle; j++)
		{
			particle[j].mu = 0.0e0;
			for(int i=0; i<FO_length; i++)
				FOsurf_ptr[i].particle_mu[j] = 0.0e0;
		}
	}
	//calculate (semi-analytic approximation of) pure thermal spectra for all particle species
	calculate_thermal_particle_yield(Nparticle, particle, FOsurf_ptr[0].Tdec);

	//use this to estimate resonance-decay contributions from each particles species to final state particle, here, pion(+),
	//as well as set effective branching ratios
	compute_total_contribution_percentages(particle_idx, Nparticle, particle);
	//sort all particles by importance of their percentage contributions, then compute resonance SVs for only contributions up to some threshold
	sw.Stop();
	output << "Read in data finished!" << endl;
	output << "Used " << sw.printTime() << " sec." << endl;

	//HBT calculations begin ...
	// Get chosen particles
	double threshold = 0.0;
	double net_fraction_resonance_contribution = 0.0;
	vector<int> chosen_resonance_indices;

	if ((int)(paraRdr->getVal("chosenParticlesMode")) == 0)				// calculate chosen resonances from threshold
	{
		threshold = paraRdr->getVal("CF_resonanceThreshold");
		output << "Working with threshold = " << threshold << endl;
		get_important_resonances(particle_idx, &chosen_resonance_indices, particle, Nparticle, threshold, net_fraction_resonance_contribution, output);
		get_all_descendants(&chosen_resonance_indices, particle, Nparticle, output);

		//chosen_resonance_indices.push_back(particle_idx);

		sort_by_mass(&chosen_resonance_indices, particle, Nparticle, output);

		for (int ii = 0; ii < (int)chosen_resonance_indices.size(); ii++)
			output << ii << "   " << chosen_resonance_indices[ii] << "   " << particle[chosen_resonance_indices[ii]].name
					<< "   ,   Gamma = " << particle[chosen_resonance_indices[ii]].width
					<< "   ,   pc = " << particle[chosen_resonance_indices[ii]].percent_contribution << endl;
	}
	else if ((int)(paraRdr->getVal("chosenParticlesMode")) == 1)		// read chosen resonances in from file
	{
		ifstream chosenParticlesStream("EOS/chosen_particles.dat");
		int len = 0;
		long tmp;
		while (chosenParticlesStream >> tmp)
		{
			chosen_resonance_indices.push_back(lookup_particle_id_from_monval(particle, Nparticle, tmp));
			len++;
		}

		chosenParticlesStream.close();
		output << "Read in " << len << " particles from EOS/chosen_particles.dat:" << endl;

		// sort chosen particles by mass
		sort_by_mass(&chosen_resonance_indices, particle, Nparticle, output);
		// erase the last element of sorted list, which will generally be the target particle
		//chosen_resonance_indices.pop_back();

		for (int ii = 0; ii < (int)chosen_resonance_indices.size(); ii++)
			output << ii << "   " << chosen_resonance_indices[ii] << "   " << particle[chosen_resonance_indices[ii]].name
					<< "   ,   Gamma = " << particle[chosen_resonance_indices[ii]].width
					<< "   ,   pc = " << particle[chosen_resonance_indices[ii]].percent_contribution << endl;

		net_fraction_resonance_contribution = 1.0;
	}
	else																// otherwise, you did something wrong!
	{
		cerr << "chosenParticlesMode = " << (int)(paraRdr->getVal("chosenParticlesMode")) << " not supported!  Exiting..." << endl;
		exit(1);
	}

//if (1) return (0);	//just for the moment

	// Create CorrelationFunction object
	CorrelationFunction correlation_function(paraRdr, &particle[particle_idx], particle, Nparticle, chosen_resonance_indices, particle_idx, output);

	// Set path and freeze-out surface information
	correlation_function.Set_path(workingDirectory);
	correlation_function.Set_FOsurf_ptr(FOsurf_ptr, FO_length);

	correlation_function.read_in_all_dN_dypTdpTdphi = false;
	correlation_function.output_all_dN_dypTdpTdphi = !(correlation_function.read_in_all_dN_dypTdpTdphi);

	correlation_function.fraction_of_resonances = net_fraction_resonance_contribution;
	output << "Using fraction_of_resonances = " << net_fraction_resonance_contribution << endl;

	//allows me to omit thermal pions easily, e.g.
	bool omit_specific_resonances = false;
	if (omit_specific_resonances)
	{
		vector<int> thermal_particles_to_omit;
		thermal_particles_to_omit.push_back(particle_idx);	//push back pion(+) to ignore thermal pions
		thermal_particles_to_omit.push_back(9);	//push back rho(+) to ignore thermal rhos (if doing phase space integration exactly)
		//double tmp = 0.0;
		//double threshhold_of_thermal_resonances_to_omit = 0.6;	//60%
		//get_important_resonances(particle_idx, &thermal_resonances_to_omit, particle, Nparticle, threshhold_of_thermal_resonances_to_omit, tmp, output);
		//get_all_descendants(&thermal_resonances_to_omit, particle, Nparticle, output);
		correlation_function.osr = thermal_particles_to_omit;
	}

//cheat...
/*
int test_pid = 1;	//rho+
double cos_res = 0.0, sin_res = 0.0;
double cos_norm = 0.0, sin_norm = 0.0;
double pT = 1.01002173, pphi = 0.0, p_Y = 0.0;
double qx = 0.15, qy = 0.0, qz = 0.0;
double m = particle[test_pid].mass;
double xi2 = m*m + pT*pT + 0.25 * (qx*qx + qy*qy + qz*qz);
double qt = sqrt(xi2 + qx*pT) - sqrt(xi2 - qx*pT);
correlation_function.Cal_dN_dypTdpTdphi_with_weights_function_approx(
	test_pid, pT, pphi, p_Y,
	0.0, 0.0, 0.0, 0.0, &cos_norm, &sin_norm);
correlation_function.Cal_dN_dypTdpTdphi_with_weights_function_approx(
	test_pid, pT, pphi, p_Y,
	qt, qx, qy, qz, &cos_res, &sin_res);
cout << "ANSWER: " << (cos_res*cos_res + sin_res*sin_res)/(cos_norm*cos_norm + sin_norm*sin_norm) << endl;
*/
/*
const int n_pts = 51;
double * dummy_pts_1 = new double [n_pts];
double * dummy_wts_1 = new double [n_pts];
double * dummy_pts_2 = new double [n_pts];
double * dummy_wts_2 = new double [n_pts];
double * dummy_pts_3 = new double [n_pts];
double * dummy_wts_3 = new double [n_pts];

double alpha = 7.0;

gauss_quadrature(n_pts, 6, 0.0, 0.0, 0.0, 0.5*alpha, dummy_pts_1, dummy_wts_1);
gauss_quadrature(n_pts, 6, 0.0, 0.0, 0.0, 1.0*alpha, dummy_pts_2, dummy_wts_2);
gauss_quadrature(n_pts, 6, 0.0, 0.0, 0.0, 2.0*alpha, dummy_pts_3, dummy_wts_3);

for (int i = 0; i < n_pts; ++i)
	cout << "GHQ: " << i
			<< "   " << dummy_pts_1[i] << "   " << dummy_wts_1[i]
			<< "   " << dummy_pts_2[i] << "   " << dummy_wts_2[i]
			<< "   " << dummy_pts_3[i] << "   " << dummy_wts_3[i]
			<< endl;

double integral1 = 0.0, integral2 = 0.0, integral3 = 0.0;
for (int i = 0; i < n_pts; ++i)
{
	integral1 += dummy_wts_1[i]*exp(-alpha*cosh(dummy_pts_1[i]));
	//integral2 += dummy_wts_2[i]*exp(-2.0*dummy_pts_2[i]*dummy_pts_2[i]);
	//integral3 += dummy_wts_3[i]*exp(-2.0*dummy_pts_3[i]*dummy_pts_3[i]);
}
cout << "CHECK: " << setw(30) << setprecision(20) << 2.0*gsl_sf_bessel_K0(alpha) << "   " << integral1 << endl;

*/


/*
		double ptr = 0.204426, phir = 0.0336253, spyr = -1.65126;
		double qt = -0.198469, qx = 0.0, qy = 0.0, qz = 0.0;
		double tempCS[4];
		tempCS[0] = 0.0, tempCS[1] = 0.0, tempCS[2] = 0.0, tempCS[3] = 0.0;
		correlation_function.Cal_dN_dypTdpTdphi_with_weights_function_approx(
				9, ptr, phir, spyr, qt, qx, qy, qz,
				&tempCS[0], &tempCS[1], &tempCS[2], &tempCS[3] );
		cout << "Check thermal routines: "
				<< ptr << "   " << phir << "   " << spyr << "   " 
				<< qt << "   " << qx << "   "
				<< qy << "   " << qz << endl
				<< "\t\t" << tempCS[0] << "   " << tempCS[1] << "   "
				<< tempCS[2] << "   " << tempCS[3] << endl;

		tempCS[0] = 0.0, tempCS[1] = 0.0, tempCS[2] = 0.0, tempCS[3] = 0.0;
		correlation_function.Cal_dN_dypTdpTdphi_with_weights_function_etas_integ(
				9, ptr, phir, spyr, qt, qx, qy, qz,
				&tempCS[0], &tempCS[1], &tempCS[2], &tempCS[3], 0.0 );
		cout << "\t\t" << tempCS[0] << "   " << tempCS[1] << "   "
				<< tempCS[2] << "   " << tempCS[3] << endl;

		tempCS[0] = 0.0, tempCS[1] = 0.0, tempCS[2] = 0.0, tempCS[3] = 0.0;
		correlation_function.Cal_dN_dypTdpTdphi_with_weights_function_etas_integ(
				9, ptr, phir, spyr, qt, qx, qy, qz,
				&tempCS[0], &tempCS[1], &tempCS[2], &tempCS[3], 1.0 );
		cout << "\t\t" << tempCS[0] << "   " << tempCS[1] << "   "
				<< tempCS[2] << "   " << tempCS[3] << endl;


if (1) exit(8);
*/

	////////////////////////////////////////////
	// Actual calculations start here...
	////////////////////////////////////////////

	//do calculations
	//decide whether to compute resonances or read them in
	if ((int)(paraRdr->getVal("calculate_CF_mode")) == 0)
	{
		sw.Stop();
		sw.Reset();
		output << "Calculating correlation function with all resonance decays (looping over qt and qz)..." << endl;

		int local_qtnpts = (int)(paraRdr->getVal("qtnpts"));
		int local_qznpts = (int)(paraRdr->getVal("qznpts"));

		correlation_function.Initialize_HDF_arrays();

		//looping in this way keeps *h5 files and total loaded memory of program small at any one time
		//#pragma omp parallel for
		for (int iqt = 0; iqt < (local_qtnpts+1)/2; ++iqt)
		for (int iqz = 0; iqz < (local_qznpts+1)/2; ++iqz)
		{
			//if (iqt + iqz > 0) continue;
			sw.Start();
			correlation_function.Fourier_transform_emission_function(iqt, iqz);
			correlation_function.Compute_phase_space_integrals(iqt, iqz);
			correlation_function.Set_target_moments(iqt, iqz);
			sw.Stop();
			output << "Completed this (iqt=" << iqt << ",iqz=" << iqz << ")-loop in " << sw.printTime() << " seconds." << endl;
			sw.Reset();
		}
	}

	//decide whether to compute correlation function or read it in
	if ((int)(paraRdr->getVal("calculate_CF_mode")) < 2)
	{
		correlation_function.Cal_correlationfunction(true);		//if we didn't compute resonance decays, must read them in from files
		correlation_function.Cal_correlationfunction(false);		//if we didn't compute resonance decays, must read them in from files
		output << "Finished calculating correlation function with all resonance decays..." << endl;
	}
	else
	{
		correlation_function.Set_correlation_function_q_pts();
		correlation_function.Allocate_CFvals();
		correlation_function.Read_in_correlationfunction();
		output << "Read in correlation function with all resonance decays." << endl;
	}

	//output results
	//decide what to output
	if ((int)(paraRdr->getVal("calculate_CF_mode")) < 2)
	{
		//included in Cal_correlationfunction() above...
		/*if ( (int)(paraRdr->getVal("qtnpts")) > 1 )
		{
			correlation_function.Output_correlationfunction();
			correlation_function.Output_correlationfunction(false);
		}*/
		if ((int)(paraRdr->getVal("calculate_CF_mode")) == 0)
		{
			correlation_function.Output_total_target_dN_dypTdpTdphi();
			correlation_function.Output_total_target_eiqx_dN_dypTdpTdphi();
			correlation_function.Output_chosen_resonances();
			correlation_function.Output_resonance_fraction();
		}
	}

	//if there's a full 3D grid to fit over, do the Gaussian fit and get the HBT radii too
	if ( (int)(paraRdr->getVal("qxnpts")) > 1
		&& (int)(paraRdr->getVal("qynpts")) > 1
		&& (int)(paraRdr->getVal("qznpts")) > 1 )
	{
		output << "Calculating HBT radii via Gaussian fit method..." << endl;
		correlation_function.Get_GF_HBTradii();
		correlation_function.Output_results(0);	//0 means do GF R2ij
		output << "Outputting lambdas..." << endl;
		correlation_function.Output_lambdas();
		output << "Finished calculating HBT radii via Gaussian fit method" << endl;
		/*output << "Calculating HBT radii via q-moments method..." << endl;
		correlation_function.Get_QM_HBTradii();
		correlation_function.Output_results(folderindex, 1);	//1 means do QM R2ij
		output << "Finished calculating HBT radii via q-moments method" << endl;*/

	 	/*if (FLESH_OUT_CF)
		{
			output << "Allocating fleshed out CFvals..." << endl;
			correlation_function.Allocate_fleshed_out_CF();
			for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
			for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
			{
				output << "Fleshing out ipt = " << ipt << ", ipphi = " << ipphi << "..." << endl;
				correlation_function.Flesh_out_CF(ipt, ipphi);
				correlation_function.Output_fleshed_out_correlationfunction(ipt, ipphi);
			}
		}*/
	}

	//do some clean up
	output << "Cleaning up..." << endl;
	for(int i = 0; i < N_stableparticle; ++i)
		delete [] particle_mu[i];
	delete [] particle_mu;

	delete [] FOsurf_ptr;
	output << "...finished cleaning up." << endl;

	sw.Stop();
	output << "Finished in " << sw.printTime() << " sec." << endl;
	sw_total.Stop();
	output << "Program totally finished in " << sw_total.printTime() << " sec." << endl;

	output << "/**********End of processing output**********/" << endl;

	output.close();

	//checkforfiles_PRfile(workingDirectory);

	finalize_PRfile(workingDirectory);

	return 0;
}

//End of file
