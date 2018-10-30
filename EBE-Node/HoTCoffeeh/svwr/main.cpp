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

#include "src/ParameterReader.h"
#include "src/Stopwatch.h"
#include "src/parameters.h"
#include "src/readindata.h"
#include "src/svwr.h"
#include "src/generate_processing_record.h"
#include "src/lib.h"
#include "main.h"

using namespace std;

int main(int argc, char *argv[])
{
   cout << endl
        << "            Source variances with resonances            " << endl
        << endl
        << "  Ver 1.0   ----- Christopher Plumberg, 07/2016   " << endl;
   cout << endl << "**********************************************************" << endl;
   display_logo(2); // Hail to the king~
   cout << endl << "**********************************************************" << endl << endl;
   
   // Read-in parameters
   ParameterReader *paraRdr = new ParameterReader;
   paraRdr->readFromFile("parameters.dat");
   paraRdr->readFromArguments(argc, argv);
   paraRdr->echo();

////////////////////////////////////////////////////////////////////
// below are Chris' set-up for computing the source variances
////////////////////////////////////////////////////////////////////

	Stopwatch sw;
	Stopwatch sw_total;
	sw_total.Start();
	sw.Start();

	string workingDirectory = "./results";
	initialize_PRfile(paraRdr, workingDirectory);

	ostringstream filename_stream;
	filename_stream << workingDirectory << "/Processing_record.txt";
	ofstream output(filename_stream.str().c_str(), ios::app);

	output << "/**********Processing output**********/" << endl;
	output << "entering folder: " << workingDirectory << endl;

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
	double** particle_mu = new double* [N_stableparticle];

	for (int i = 0; i < N_stableparticle; i++)
		particle_mu[i] = new double [FO_length];

	for (int i = 0; i < N_stableparticle; i++)
	for (int j = 0; j < FO_length; j++)
		particle_mu[i][j] = 0.0;
	if (N_stableparticle > 0)
	{
		//if(hydropara_ptr->IEOS==7)       //for s95p_PCE
		read_decdat_mu(FO_length, N_stableparticle, particle_mu, workingDirectory);
	}

	//read particle resonance decay table
	particle_info *particle = new particle_info [Maxparticle];
	int Nparticle=read_resonance(particle);
	output <<"read in total " << Nparticle << " particles!" << endl;
	output << "Calculating "<< particle[particle_idx].name << endl;
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
   double localy = 0.0e0;
   sw.Reset();
   sw.Start();
     if(fabs(localy) > 1e-16)
     {
        output << "Case of y != 0 not yet supported.  Exiting..." << endl;
        return 0;
     }

	// Get chosen particles
	double threshold = 0.0;
	double net_fraction_resonance_contribution = 0.0;
	vector<int> chosen_resonance_indices;

	if ((int)(paraRdr->getVal("chosenParticlesMode")) == 0)				// calculate chosen resonances from threshold
	{
		threshold = paraRdr->getVal("SV_resonanceThreshold");
		output << "Working with threshold = " << threshold << endl;
		get_important_resonances(particle_idx, &chosen_resonance_indices, particle, Nparticle, threshold, net_fraction_resonance_contribution, output);
		get_all_descendants(&chosen_resonance_indices, particle, Nparticle, output);
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
		while (!chosenParticlesStream.eof())
		{
			long tmp;
			chosenParticlesStream >> tmp;
			chosen_resonance_indices.push_back(lookup_particle_id_from_monval(particle, Nparticle, tmp));
			len++;
		}

		chosenParticlesStream.close();
		output << "Read in " << len << " particles from EOS/chosen_particles.dat:" << endl;

		// sort chosen particles by mass
		sort_by_mass(&chosen_resonance_indices, particle, Nparticle, output);
		for (int ii = 0; ii < (int)chosen_resonance_indices.size(); ii++)
			output << ii << "   " << chosen_resonance_indices[ii] << "   " << particle[chosen_resonance_indices[ii]].name
					<< "   ,   Gamma = " << particle[chosen_resonance_indices[ii]].width
					<< "   ,   pc = " << particle[chosen_resonance_indices[ii]].percent_contribution << endl;

		// erase the last element of sorted list, which will generally include the target particle
		chosen_resonance_indices.erase (chosen_resonance_indices.end());
		net_fraction_resonance_contribution = 1.0;
	}
	else																// otherwise, you did something wrong!
	{
		cerr << "chosenParticlesMode = " << (int)(paraRdr->getVal("chosenParticlesMode")) << " not supported!  Exiting..." << endl;
		exit(1);
	}

//if (1) return(0);

	// Create SourceVariances object
	SourceVariances Source_function(paraRdr, &particle[particle_idx], particle, Nparticle, chosen_resonance_indices, particle_idx, output);
	
	// Set path and freeze-out surface information
	Source_function.Set_path(workingDirectory);
	Source_function.Set_FOsurf_ptr(FOsurf_ptr, FO_length);

	//read in space-time moments from file (default: false)
	Source_function.read_in_all_dN_dypTdpTdphi = false;
	Source_function.output_all_dN_dypTdpTdphi = !(Source_function.read_in_all_dN_dypTdpTdphi);

	// Do source variances calculations
	output << "Calculating HBT radii via source variances method..." << endl;
	Source_function.Analyze_sourcefunction();

	Source_function.Output_emission_density(particle_idx);

	// Output most interesting results
	Source_function.Output_total_target_dN_dypTdpTdphi();
	Source_function.Output_chosen_resonances();
	Source_function.Output_results();
	output << "Finished calculating HBT radii via source variances method" << endl;



   sw.Stop();
   output << "Finished in " << sw.printTime() << " sec." << endl;
   sw_total.Stop();
   output << "Program totally finished in " << sw_total.printTime() << " sec." << endl;

   output << "/**********End of processing output**********/" << endl;

   output.close();

   //checkforfiles_PRfile(workingDirectory, folderindex, generatedcorrfuncs);

   finalize_PRfile(workingDirectory);

   return 0;
}
