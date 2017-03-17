#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<cstdlib>
#include<gsl/gsl_sf_bessel.h>
#include<vector>
#include<algorithm>

#include "readindata.h"
using namespace std;

double * stable_particle_monval;
int Nstable_particle;
double ** all_b_j_to_i;
bool output_mu_and_yield = false;

void read_hydropar(hydropara* hp, string localpath)
{
  ostringstream hydropar_stream;
  hydropar_stream << localpath << "/VISH2p1_tec.dat";
  
  cout << "read " << hydropar_stream.str().c_str() << "....";
  char dummy[256];
  char delim = '=';
  ifstream if_hydropar(hydropar_stream.str().c_str());
  
  if_hydropar.getline(dummy, 256, delim);
  if_hydropar >> hp->IEOS;
  if_hydropar.getline(dummy, 256);  //line 1
  if_hydropar.getline(dummy, 256, delim);
  if_hydropar >> hp->tau_0;
  if_hydropar.getline(dummy, 256);  //line 2
  if_hydropar.getline(dummy, 256, delim);
  if_hydropar >> hp->e_dec;
  if_hydropar.getline(dummy, 256);  //line 3
  if_hydropar.getline(dummy, 256, delim);
  if_hydropar >> hp->T_dec;
  if_hydropar.getline(dummy, 256);  //line 4
  if_hydropar.getline(dummy, 256);  //line 5
  if_hydropar.getline(dummy, 256);  //line 6
  if_hydropar.getline(dummy, 256);  //line 7
  if_hydropar.getline(dummy, 256);  //line 8
  if_hydropar.getline(dummy, 256, delim);
  if_hydropar >> hp->eta_s;
  if_hydropar.getline(dummy, 256);  //line 9
  if_hydropar.getline(dummy, 256);  //line 10
  if_hydropar.getline(dummy, 256);  //line 11
  if_hydropar.getline(dummy, 256, delim);
  if_hydropar >> hp->zeta_s;
  if_hydropar.getline(dummy, 256);  //line 12
  if_hydropar.getline(dummy, 256);  //line 13
  if_hydropar.getline(dummy, 256);  //line 14
  if_hydropar.getline(dummy, 256, delim);
  if_hydropar >> hp->dtau;
  if_hydropar.getline(dummy, 256);  //line 15
  if_hydropar.getline(dummy, 256, delim);
  if_hydropar >> hp->dx;
  if_hydropar.getline(dummy, 256);  //line 16
  if_hydropar.getline(dummy, 256, delim);
  if_hydropar >> hp->dy;
  if_hydropar.getline(dummy, 256);  //line 17
  if_hydropar.getline(dummy, 256, delim);
  if_hydropar >> hp->NLS;
  if_hydropar.getline(dummy, 256);  //line 18
  if_hydropar.getline(dummy, 256, delim);
  if_hydropar >> hp->NXD;
  if_hydropar.getline(dummy, 256);  //line 19
  if_hydropar.getline(dummy, 256, delim);
  if_hydropar >> hp->NYD;
  if_hydropar.getline(dummy, 256);  //line 20
  if_hydropar.getline(dummy, 256, delim);
  if_hydropar >> hp->NTauD;         
  if_hydropar.getline(dummy, 256);  //line 21
  if_hydropar.getline(dummy, 256);  //line 22
  if_hydropar.close();
  cout << "done!" << endl;

  cout << "IEOS= " << hp->IEOS << endl;
  cout << "tau_0= " << hp->tau_0 << endl;
  cout << "T_dec= " << hp->T_dec << endl;
  cout << "eta/s= " << hp->eta_s << endl;
  cout << "zeta/s= " << hp->zeta_s << endl;
  cout << "NXD= " << hp->NXD << endl;
  cout << "NYD= " << hp->NYD << endl;
  cout << "NTauD= " << hp->NTauD << endl;
  cout << "dtau= " << hp->dtau << endl;
  cout << "dx= " << hp->dx << endl;
  cout << "dy= " << hp->dy << endl;
  cout << "NLS= " << hp->NLS << endl;
  return;
}

int get_filelength(string filepath)
{
   int length=0; 
   char line[512];
   ostringstream filepath_stream;
   filepath_stream << filepath;
   ifstream infile(filepath_stream.str().c_str());
   //Determine the length of the file
   while (!infile.eof ())
   {
      infile.getline(line, 512);
      length++;
   }
   length = length-1;
   infile.close();
   return(length);
}

int get_filewidth(string filepath)
{
	ostringstream filepath_stream;
	filepath_stream << filepath;
	ifstream infile(filepath_stream.str().c_str());
	string line, temp;
	stringstream ss;
	int ncols=0;
	getline(infile, line);
	ss.clear();
	ss << line;
    
	while (ss >> temp)
		ncols++;

	return ncols;
}

void read_decdat(int length, FO_surf* surf_ptr, string localpath, bool include_bulk_pi /* = false*/)
{ 
  //cout<<"read in information on freeze out surface...";
  ostringstream decdat_stream;
  decdat_stream << localpath << "/decdat2.dat";
  ifstream decdat(decdat_stream.str().c_str());
  for(int i=0; i<length; i++)
  {
     decdat >> surf_ptr[i].tau;
     decdat >> surf_ptr[i].da0;
     decdat >> surf_ptr[i].da1;
     decdat >> surf_ptr[i].da2;
     decdat >> surf_ptr[i].vx;
     decdat >> surf_ptr[i].vy;
     decdat >> surf_ptr[i].Edec;
     decdat >> surf_ptr[i].Bn;
     decdat >> surf_ptr[i].Tdec;
     decdat >> surf_ptr[i].muB;
     decdat >> surf_ptr[i].muS;
     decdat >> surf_ptr[i].Pdec;
     decdat >> surf_ptr[i].pi33;
     decdat >> surf_ptr[i].pi00;
     decdat >> surf_ptr[i].pi01;
     decdat >> surf_ptr[i].pi02;
     decdat >> surf_ptr[i].pi11;
     decdat >> surf_ptr[i].pi12;
     decdat >> surf_ptr[i].pi22;
     if (include_bulk_pi) decdat >> surf_ptr[i].bulkPi;
     surf_ptr[i].gammaT = 1./sqrt(1.- surf_ptr[i].vx*surf_ptr[i].vx - surf_ptr[i].vy*surf_ptr[i].vy);
  }
  decdat.close();
  //cout<<"done"<<endl;
  return;
}

void read_surfdat(int length, FO_surf* surf_ptr, string localpath)
{
  //cout<<"read spacial positions of freeze out surface...";
  ostringstream surfdat_stream;
  double dummy;
  char rest_dummy[512];
  surfdat_stream << localpath << "/surface.dat";
  ifstream surfdat(surfdat_stream.str().c_str());
  for(int i=0; i<length; i++)
  {
     surfdat >> dummy >> dummy;
     surfdat >> surf_ptr[i].xpt;
     surfdat >> surf_ptr[i].ypt;
     surf_ptr[i].r = sqrt(surf_ptr[i].xpt*surf_ptr[i].xpt + surf_ptr[i].ypt*surf_ptr[i].ypt);
     surf_ptr[i].phi = atan2(surf_ptr[i].ypt, surf_ptr[i].xpt);
     surf_ptr[i].sin_phi = sin(surf_ptr[i].phi);
     surf_ptr[i].cos_phi = cos(surf_ptr[i].phi);
     surfdat.getline(rest_dummy, 512);
  }
  surfdat.close();
  //cout<<"done"<<endl;
  return;
}

void read_decdat_mu(int FO_length, int N_stable, double** particle_mu, string localpath)
{
  //cout<<" -- Read chemical potential for stable particles...";
  ostringstream decdat_mu_stream;
  double dummy;
  decdat_mu_stream << localpath << "/decdat_mu.dat";
  ifstream decdat_mu(decdat_mu_stream.str().c_str());

  //For backward compatibility: decdat_mu.dat can be one line or FO_length lines
  for(int j=0; j<FO_length; j++)
  {
    decdat_mu >> dummy;  //not used in the code plz ignore it

    if(decdat_mu.eof())
    {
      for(int k=j; k<FO_length; k++)
        for(int i=0; i<N_stable; i++)
           particle_mu[i][k]=particle_mu[i][j-1];
      break;
    }

    for(int i=0; i<N_stable; i++)
    {
       decdat_mu >> particle_mu[i][j];
    }
  }

  //cout<<"done" << endl;
  return;
}

int read_resonance(particle_info* particle)
{
   int Nparticle=0; 
   //cout << "Reading in particle resonance decay table...";
   ifstream resofile("/home/plumberg.1/HBTPlumberg/EOS/pdg.dat");
   int local_i = 0;
   int dummy_int;
   while (!resofile.eof())
   {
      resofile >> particle[local_i].monval;
      resofile >> particle[local_i].name;
      resofile >> particle[local_i].mass;
      resofile >> particle[local_i].width;
      resofile >> particle[local_i].gspin;	      //spin degeneracy
      resofile >> particle[local_i].baryon;
      resofile >> particle[local_i].strange;
      resofile >> particle[local_i].charm;
      resofile >> particle[local_i].bottom;
      resofile >> particle[local_i].gisospin;     //isospin degeneracy
      resofile >> particle[local_i].charge;
      resofile >> particle[local_i].decays;
      for (int j = 0; j < particle[local_i].decays; j++)
      {
         resofile >> dummy_int;
         resofile >> particle[local_i].decays_Npart[j];
         resofile >> particle[local_i].decays_branchratio[j];
         resofile >> particle[local_i].decays_part[j][0];
         resofile >> particle[local_i].decays_part[j][1];
         resofile >> particle[local_i].decays_part[j][2];
         resofile >> particle[local_i].decays_part[j][3];
         resofile >> particle[local_i].decays_part[j][4];
      }
      
      //decide whether particle is stable under strong interactions
      if(particle[local_i].decays_Npart[0] == 1) 	      
         particle[local_i].stable = 1;
      else
         particle[local_i].stable = 0;

      //add anti-particle entry
      if(particle[local_i].baryon == 1)
      {
         local_i++;
         particle[local_i].monval = -particle[local_i-1].monval;
         ostringstream antiname;
         antiname << "Anti-" << particle[local_i-1].name;
         particle[local_i].name = antiname.str();
         particle[local_i].mass = particle[local_i-1].mass;
         particle[local_i].width = particle[local_i-1].width;
         particle[local_i].gspin = particle[local_i-1].gspin;
         particle[local_i].baryon = -particle[local_i-1].baryon;
         particle[local_i].strange = -particle[local_i-1].strange;
         particle[local_i].charm = -particle[local_i-1].charm;
         particle[local_i].bottom = -particle[local_i-1].bottom;
         particle[local_i].gisospin = particle[local_i-1].gisospin;
         particle[local_i].charge = -particle[local_i-1].charge;
         particle[local_i].decays = particle[local_i-1].decays;
         particle[local_i].stable = particle[local_i-1].stable;
         for (int j = 0; j < particle[local_i].decays; j++)
         {
            particle[local_i].decays_Npart[j]=particle[local_i-1].decays_Npart[j];
            particle[local_i].decays_branchratio[j]=particle[local_i-1].decays_branchratio[j];
            //for (int k=0; k< Maxdecaypart; k++)						//commented out by Chris Plumberg - 06/08/2015
            //   particle[local_i].decays_part[j][k]=particle[local_i-1].decays_part[j][k];
            for (int k=0; k< Maxdecaypart; k++)							//replaced with following k-loop
            {
               int idx = 0;  
               for(int ii=0; ii < local_i; ii++) // find the index for decay particle
               {
                  if(particle[local_i-1].decays_part[j][k] == particle[ii].monval)
                  {
                     idx = ii;
                     break;
                  }
               }
               if(idx == local_i-1 && particle[local_i-1].stable == 0)  // check
               {
                  cout << "Error: can not find decay particle index for anti-baryon!" << endl;
                  cout << "particle monval : " << particle[local_i-1].decays_part[j][k] << endl;
                  exit(1);
               }
               if(particle[idx].baryon == 0 && particle[idx].charge == 0 && particle[idx].strange == 0)
                  particle[local_i].decays_part[j][k]= particle[local_i-1].decays_part[j][k];
               else
                  particle[local_i].decays_part[j][k]= -particle[local_i-1].decays_part[j][k];
            } // end of k-loop
         }
       }
       local_i++;	// Add one to the counting variable "i" for the meson/baryon
   }
   resofile.close();
   Nparticle=local_i-1; //take account the final fake one
   for(int i=0; i < Nparticle; i++)
   {
      if(particle[i].baryon==0)
         particle[i].sign=-1;
      else
         particle[i].sign=1;
   }
   //cout << "done! Antiparticles are added!" << endl;
   return(Nparticle);
}

void calculate_particle_mu(int IEOS, int Nparticle, FO_surf* FOsurf_ptr, int FO_length, particle_info* particle, double** particle_mu)
{
	Nstable_particle = set_stable_particle_monval();
      for(int i=0; i<Nstable_particle; i++)
         for(int j=0; j<Nparticle; j++)
            if(particle[j].monval == stable_particle_monval[i])
            {
               particle[j].stable = 1;
               for(int k=0; k<FO_length; k++)
                   FOsurf_ptr[k].particle_mu[j] = particle_mu[i][k];
               break;
            }

      for(int i=0; i < Nparticle ; i++)
      {
         if(particle[i].stable==0)
         {
            for(int j=0; j < particle[i].decays; j++)
            {
               for(int k=0; k < abs(particle[i].decays_Npart[j]); k++)
               {
                  for(int l=0; l < Nparticle; l++)
                  {
                     if(particle[i].decays_part[j][k] == particle[l].monval)
                     {
                        for(int m=0; m<FO_length; m++)
                          FOsurf_ptr[m].particle_mu[i] += particle[i].decays_branchratio[j]*FOsurf_ptr[m].particle_mu[l];
                        break;
                     }
                     if(l==Nparticle-1)
                        cout<<"warning: can not find particle" <<  particle[i].name << endl;
                  }
               }
            }
         }
      }
	for(int i = 0; i < Nparticle; i++)	//set particle.mus for easy access --> assumes mu is constant along FO surface
		particle[i].mu = FOsurf_ptr[0].particle_mu[i];
   return;
}

void calculate_thermal_particle_yield(int Nparticle, particle_info* particle, double Temperature)
{
	double one_by_Tconv = 1./Temperature;
	//double * all_particle_fugacities = new double [Maxparticle];
	//set_to_zero(all_particle_thermal, Nparticle);
	for (int j = 0; j < Nparticle; j++)
	{
		double yield = 0.0, fugacity = 0.0;
		double gj = particle[j].gspin;
		double mj = particle[j].mass;
		if (mj == 0)
		{
			fugacity = 0.0;
			particle[j].thermal_yield = 0.0;
			continue;
		}
		double pm = -particle[j].sign;
		fugacity = exp(one_by_Tconv * particle[j].mu);
		for (int k = 1; k <= 10; k++)
			yield += double(pow(pm, k+1))*pow(fugacity, (double)k)*gsl_sf_bessel_Kn(2, double(k)*mj*one_by_Tconv)/double(k);
		particle[j].thermal_yield = yield*gj*mj*mj/(2.*M_PI*M_PI);
	}
	//**********************************************************************
	if (output_mu_and_yield)
	{
		ostringstream output_stream;
		output_stream << "check_mu_and_yield.dat";
		ofstream output(output_stream.str().c_str());
		for (int j = 0; j < Nparticle; j++)
			output << particle[j].name << "   " << particle[j].mu << "   " << particle[j].thermal_yield << endl;
		output.close();
	}
	//**********************************************************************
	
	return;
}

void compute_total_contribution_percentages(int stable_particle_idx, int Nparticle, particle_info* particle)
{
	double denominator = 0.0, temp;
	all_b_j_to_i = new double * [Nparticle];
	for (int i = 0; i < Nparticle; i++)
	{
		particle[i].percent_contribution = 0.0;
		all_b_j_to_i[i] = new double [Nparticle];
		for (int j = 0; j < Nparticle; j++)
			all_b_j_to_i[i][j] = 0.0;
	}
	for (int i = 0; i < Nparticle; i++)
	{
		temp = b_j_to_i(particle, Nparticle, i, stable_particle_idx);
		particle[i].percent_contribution = temp * particle[i].thermal_yield;
		particle[i].effective_branchratio = temp;
		denominator += temp * particle[i].thermal_yield;
	}
	for (int i = 0; i < Nparticle; i++)
		particle[i].percent_contribution /= 0.01 * denominator;	//0.01 factor makes it a percentage
	
	return;
}

void calculate_pTdep_thermal_particle_yield(int Nparticle, particle_info* particle, double Temperature, double * pTarray, int npT, double ** resultsArray)
{
	double one_by_Tconv = 1./Temperature;
	for (int j = 0; j < Nparticle; j++)
	{
		double fugacity = exp(one_by_Tconv * particle[j].mu);
		double gj = particle[j].gspin;
		double mj = particle[j].mass;
		double pm = -particle[j].sign;
		if (mj == 0)
			continue;
		for (int i = 0; i < npT; i++)
		{
			double pT = pTarray[i];
			double mTj = sqrt(mj*mj + pT*pT);
			double yield = 0.0;
			for (int k = 1; k <= 10; k++)
				yield += double(pow(pm, k+1))*pow(fugacity, (double)k)*gsl_sf_bessel_Kn(1, double(k)*mTj*one_by_Tconv)
											*gsl_sf_bessel_I0 (double(k)*mTj*one_by_Tconv);
			//resultsArray[j][i] = yield*gj*mTj*mTj/(M_PI*M_PI);
			resultsArray[j][i] = yield*gj*mTj/(M_PI*M_PI);
		}
	}
	
	return;
}

void compute_total_pTdep_contribution_percentages(int stable_particle_idx, int Nparticle, particle_info * particle, int npT,
						double ** pTdepThermalYields, vector<int> * chosen_resonance_indices_ptr, double * pTdep_Cumulative_Resonance_Contributions)
{
	for (int j = 0; j < npT; j++)
	{
		double denominator = 0.0, temp;
		pTdep_Cumulative_Resonance_Contributions[j] = 0.0;
		
		// get total pT-dependent contributions to stable particle (e.g., pions) from all resonances
		for (int i = 0; i < Nparticle; i++)
			denominator += all_b_j_to_i[i][stable_particle_idx] * pTdepThermalYields[i][j];

		// get total pT-dependent contributions to stable particle (e.g., pions) from just the chosen resonances
		for (int i = 0; i < (*chosen_resonance_indices_ptr).size(); i++)
		{
			int idx = (*chosen_resonance_indices_ptr)[i];
			pTdep_Cumulative_Resonance_Contributions[j] += all_b_j_to_i[idx][stable_particle_idx] * pTdepThermalYields[idx][j] / denominator;
		}
	}
	
	return;
}


double b_j_to_i(particle_info * particle, int Nparticle, int j, int i, int verbose_monval /* = 0*/)
{
	double result = 0.0;
	particle_info parent = particle[j];
	particle_info target = particle[i];
	int verbose = 0;
	if (parent.monval == verbose_monval) verbose = 1;
	if (verbose > 0) cout << "Currently looking at decay chains of " << parent.name << " to " << target.name << ":" << endl;
	if ( (parent.decays == 1) && (parent.decays_Npart[0] == 1) )
	{
		if (verbose > 0) cout << "   --> " << parent.name << " is stable, so moving on!" << endl;
		return(0.0);	// if parent is already stable, return 0
	}
	else
	{
		if (verbose > 0) cout << "   --> " << parent.name << " is unstable, so continuing!" << endl;
		//particle[j].stable == 0;	//just in case
	}
	
	// added this to try to save some time and simplify debugging output
	if (fabs(all_b_j_to_i[j][i]) > 1.e-12)
	{
		if (verbose > 0) cout << "   --> already calculated!  Recycling stored value and moving on!" << endl;
		if (verbose > 0) cout << "   --> " << parent.name << "-->" << target.name << ": using b_j_to_i = " << all_b_j_to_i[j][i] << endl;
		return all_b_j_to_i[j][i];
	}
	
	for (int k = 0; k < parent.decays; k++)
	{
		int nki = count_targets(parent.decays_part[k], &target);	// number of target decay particles in kth decay channel
		double bk = parent.decays_branchratio[k];			// branching ratio for kth decay channel
		int nks = count_stable(particle, Nparticle, parent.decays_part[k]);			// number of stable decay particles in kth decay channel
		int nktot = abs(parent.decays_Npart[k]);				// total number of decay particles in kth decay channel
		if (verbose > 0) cout << " - " << parent.name << "(monval = " << parent.monval << "): decay channel " << k + 1 << " of " << parent.decays << endl;
		if (verbose > 0) cout << "   --> nki = " << nki << ", nks = " << nks << ", nktot = " << nktot << endl;
		if ((nki == 0) && (nks == nktot))
		{
			if (verbose > 0) cout << "   --> found no " << target.name << "s or unstable particles, so moving on!" << endl;
			continue;			// if kth decay channel contains no target particles or other particles which might decay to some
		}
		if (nki != 0) result += double(nki)*bk;				// if kth decay channel contains target particles
		parent.decays_effective_branchratio[k] = double(nki)*bk;
		particle[j].decays_effective_branchratio[k] = double(nki)*bk;
		if (nks != nktot)						// if there are unstable particles
		{
			if (verbose > 0) cout << "   --> found some unstable particles!" << endl;
			for (int ipart = 0; ipart < nktot; ipart++)
			{		// apply this same process recursively to all unstable daughter resonances
				if ( !is_stable( particle, Nparticle, parent.decays_part[k][ipart] ) )
				{
					int decay_particle_idx = lookup_particle_id_from_monval(particle, Nparticle, parent.decays_part[k][ipart]);
					if (verbose > 0) cout << "   --> now working on unstable particle (" << particle[decay_particle_idx].name << ") with monval = "
						<< parent.decays_part[k][ipart] << endl << endl;
					double temp_bj2i = b_j_to_i(particle, Nparticle, decay_particle_idx, i);
					if (verbose > 0) cout << "   --> " << parent.name << "-->" << target.name << ": using b_j_to_i = " << temp_bj2i << endl;
					result += bk * temp_bj2i;
					parent.decays_effective_branchratio[k] += bk * temp_bj2i;
					particle[j].decays_effective_branchratio[k] += bk * temp_bj2i;
				}
			}
		}
		
	}
	
	all_b_j_to_i[j][i] = result;
	if (verbose > 0) cout << "***FINAL***: " << parent.name << "-->" << target.name << ": using b_j_to_i = " << all_b_j_to_i[j][i] << endl;
	return (result);
}

int lookup_particle_id_from_monval(particle_info * particle, int Nparticle, int monval)
{
	for (int ipart = 0; ipart < Nparticle; ipart++)
	{
		//cout << "lookup(" << monval << "): " << particle[ipart].name << "   " << particle[ipart].monval << endl;
		if (monval == particle[ipart].monval) return ipart;
	}
	cerr << "monval = " << monval << endl;
	cerr << "Only available monvals are:" << endl;
	for (int ipart = 0; ipart < Nparticle; ipart++)
		cerr << particle[ipart].name << "   " << particle[ipart].monval << endl;
	cerr << "Could not find monval in PDG table!  Aborting..." << endl;
	exit(1);
}

int count_targets(int * decay_channel_particles, particle_info * i)
{
	int count = 0;
	for (int idcp = 0; idcp < Maxdecaypart; idcp++)
		if (i->monval == decay_channel_particles[idcp]) count++;
	return(count);
}

int count_stable(particle_info * particle, int Nparticle, int * decay_channel_particles)
{
	int count = 0;
	for (int idcp = 0; idcp < Maxdecaypart; idcp++)
	{
		if (decay_channel_particles[idcp] == 0) continue;
		if (is_stable(particle, Nparticle, decay_channel_particles[idcp])) count++;
	}
	return(count);
}

bool is_stable(particle_info * particle, int Nparticle, int monval)
{
	int local_idx = lookup_particle_id_from_monval(particle, Nparticle, monval);
	return(particle[local_idx].decays_Npart[0] == 1);
}

int set_stable_particle_monval()
{
   int local_Nstable_particle;
   int Idummy;
   char cdummy[256];
   //if(IEOS!=7)
   //{
   //   cout << "Error! IEOS = "<<IEOS << " is not support for PCE!"<<endl;
   //   exit(1);
   //}
   //else
   //{
      //cout << "Reading particle table and calculating chemical potential for particles...";
      //ifstream particletable("EOS/EOS_particletable.dat");
      ifstream particletable("/home/plumberg.1/HBTPlumberg/EOS/EOS_particletable.dat");
      particletable >> local_Nstable_particle;
      stable_particle_monval = new double [local_Nstable_particle];
      for(int i=0; i<local_Nstable_particle; i++)
      {
          particletable >> Idummy >> stable_particle_monval[i];
          particletable.getline(cdummy, 256);
      }
      particletable.close();
   //}
	return(local_Nstable_particle);
}

void print_particle_stability(particle_info * particle, int Nparticle)
{
	for (int ipart = 0; ipart < Nparticle; ipart++)
		cout << particle[ipart].name << "   " << particle[ipart].monval << "   " << particle[ipart].mass << "   " << particle[ipart].stable << endl;
	return;
}

int get_number_of_decay_channels(vector<int> chosen_resonances, particle_info * particle)
{
	int count = 0, total_number_of_decays = 0;
	//cerr << "get_number_of_decay_channels(): (int)chosen_resonances.size() = " << (int)chosen_resonances.size() << endl;
	for (int icr = 0; icr < (int)chosen_resonances.size(); icr++)
	{
		//cerr << "get_number_of_decay_channels(): accessing icr = " << icr << ", chosen_resonances[" << icr << "] = " << chosen_resonances[icr] << endl;
		total_number_of_decays = particle[chosen_resonances[icr]].decays;
		//cerr << "get_number_of_decay_channels(): total_number_of_decays = " << total_number_of_decays << endl;
		//for (int idecay = 0; idecay < total_number_of_decays; idecay++)
		//	count += (particle[chosen_resonances[icr]].decays_effective_branchratio[idecay] > 1e-12) ? 1 : 0;
		count += total_number_of_decays;
		//cerr << "get_number_of_decay_channels(): count = " << count << endl;
	}
	return count;
}



void get_important_resonances(int chosen_target_particle_idx, vector<int> * chosen_resonance_indices_ptr, particle_info * particle, int Nparticle, double threshold, double &running_total_percentage, std::ofstream& output)
{
	//**********************************************************************************
	//SELECT RESONANCES TO INCLUDE IN SOURCE VARIANCES CALCULATIONS
	//sort resonances by importance and loop over all resonances needed to achieve a certain minimum percentage of total decay pions
	//double threshold = 0.6;	//include only enough of the most important resonances to account for fixed fraction of total resonance-decay pion(+)s
					//threshold = 1.0 means include all resonance-decay pion(+)s,
					//threshold = 0.0 means include none of them
	vector<double> percent_contributions;
	for (int i = 0; i < Nparticle; i++)
		percent_contributions.push_back(particle[i].percent_contribution);
	vector<size_t> sorted_resonance_indices = ordered(percent_contributions);
	reverse(sorted_resonance_indices.begin(), sorted_resonance_indices.end());
	//vector<int> chosen_resonance_indices_ptr;
	//double running_total_percentage = 0.0;
	int count = 0;
	if (threshold < 1e-12)
	{
		output << "No resonances included." << endl;
	}
	else if (fabs(1. - threshold) < 1e-12)
	{
		count = Nparticle;	//			if (sorted_resonance_indices[ii - 1] == chosen_target_particle_idx)
		for (int ii = 1; ii <= count; ii++)
			(*chosen_resonance_indices_ptr).push_back(sorted_resonance_indices[ii - 1]);
		running_total_percentage = 1.0;
		output << "All resonances included." << endl;
	}
	else
	{
		while (running_total_percentage <= threshold)
		{
			running_total_percentage += 0.01 * particle[sorted_resonance_indices[count]].percent_contribution;
			(*chosen_resonance_indices_ptr).push_back(sorted_resonance_indices[count]);
			count++;
		}
		if ((*chosen_resonance_indices_ptr).size() == 0)
		{
			output << "No resonances included!  Choose a higher threshold!" << endl;
			exit(1);
		}
		else
		{
			output << "Including the following " << count << " resonances (accounting for " << 100.*running_total_percentage
				<< "%, threshold = " << 100.*threshold << "%): " << endl;
			for (int ii = 1; ii <= count; ii++)
				output << "\t" << ii << ": " << particle[sorted_resonance_indices[ii - 1]].name << endl;
		}
	}
	//END OF CODE SELECTING INCLUDED RESONANCES
	//**********************************************************************************
	
	return;
}

void get_all_descendants(vector<int> * chosen_resonance_indices_ptr, particle_info * particle, int Nparticle, std::ofstream& output)
{
	// This function ensures that no particles are missing from the chosen_resonances vector,
	// even if a truncated (threshold < 1.0) version is used
	int amount_of_output = 1;
	bool no_missing_descendants = false;
	int original_size = (int)(*chosen_resonance_indices_ptr).size();
	while (!no_missing_descendants)
	{
		int old_size = (int)(*chosen_resonance_indices_ptr).size();
		int new_size = old_size;
		for (int icr = 0; icr < old_size; icr++)
		{
			particle_info resonance = particle[(*chosen_resonance_indices_ptr)[icr]];
			if (amount_of_output > 1) output << "Currently looking at resonance " << resonance.name << endl;
			if (resonance.stable == 1 && resonance.decays_Npart[0] == 1)
				continue;
			int number_of_decays = resonance.decays;
			for (int k = 0; k < number_of_decays; k++)
			{
				if (amount_of_output > 1) output << resonance.name << ": currently looking at decay #" << k << endl;
				int nb = abs(resonance.decays_Npart[k]);
				for (int l = 0; l < nb; l++)
				{
					int pid = lookup_particle_id_from_monval(particle, Nparticle, resonance.decays_part[k][l]);
					if (amount_of_output > 1) output << resonance.name << ", decay #" << k
								<< ": currently looking at daughter " << particle[pid].name << endl;
					//if this decay particle is not in the list and has large effective br, include it
					bool decay_particle_is_not_in_list = ( find ((*chosen_resonance_indices_ptr).begin(), (*chosen_resonance_indices_ptr).end(), pid)
															==
															(*chosen_resonance_indices_ptr).end() );
					bool br_tot_is_not_small = ( particle[pid].effective_branchratio >= 1.e-12 );
					if (decay_particle_is_not_in_list)
						if (amount_of_output > 1) output << resonance.name << ", decay #" << k
							<< ", daughter " << particle[pid].name << ": daughter is not in list!" << endl;
					if (br_tot_is_not_small)
						if (amount_of_output > 1) output << resonance.name << ", decay #" << k
							<< ", daughter " << particle[pid].name << ": effective br. is not small!" << endl;
					if (decay_particle_is_not_in_list && br_tot_is_not_small)
					{
						(*chosen_resonance_indices_ptr).push_back(pid);
						if (amount_of_output > 1) output << resonance.name << ", decay #" << k
							<< ", daughter " << particle[pid].name << ": adding daughter to list!" << endl;
						new_size++;
					}
				}
			}
		}
		no_missing_descendants = ( old_size == new_size );
		if (new_size == original_size)
			if (amount_of_output > 0) output << "get_all_descendants(): No new particles added!" << endl;
		else
			if (amount_of_output > 0) output << "get_all_descendants(): " << new_size - original_size << " new particles added!" << endl;
	}
	return;
}

void sort_by_mass(vector<int> * chosen_resonance_indices_ptr, particle_info * particle, int Nparticle, std::ofstream& output)
{	//with the heaviest first
	output << "sort_by_mass(): Sorting by mass..." << endl;
	int number_of_chosen_particles = (int)(*chosen_resonance_indices_ptr).size();
	for (int m = 0; m < number_of_chosen_particles; m++)
	for (int n = 0; n < number_of_chosen_particles - m - 1; n++)
	if (particle[(*chosen_resonance_indices_ptr)[n]].mass < particle[(*chosen_resonance_indices_ptr)[n + 1]].mass)
	{
		// swap them
		int particle_idx = (*chosen_resonance_indices_ptr)[n + 1];
		(*chosen_resonance_indices_ptr)[n + 1] = (*chosen_resonance_indices_ptr)[n];
		(*chosen_resonance_indices_ptr)[n] = particle_idx;
	}
	output << "sort_by_mass(): ...finished!" << endl;
	return;
}

//End of file
