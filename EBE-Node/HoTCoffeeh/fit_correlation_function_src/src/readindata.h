#ifndef READINDATA_H
#define READINDATA_H

#include "parameters.h"
#include "sorter.h"
#include<fstream>
#include<vector>
using namespace std;

typedef struct
{
   int IEOS;            //EOS selector
   double b;            //Impact parameter
   double tau_0;        //thermalization time
   double e_dec;        //kinetic decoupling energy density
   double T_dec;
   double eta_s;        //specific shear viscosity
   double zeta_s;       //specific bulk viscosity
   int NXD;             //number of skipped points along x direction
   int NYD;             //number of skipped points along y direction
   int NTauD;           //number of skipped points along tau direction
   double dtau;         //lattice spacing in tau direction
   double dx;           //lattice spacing in x direction
   double dy;           //lattice spacing in y direction
   int NLS;             //number of lattice size (assume square lattice NX=NY)
}hydropara;

typedef struct
{
   double tau, xpt, ypt, r, phi, sin_phi, cos_phi;
   double da0, da1, da2;
   double vx, vy, gammaT;
   double Edec, Tdec, Pdec;
   double Bn, muB, muS;
   double pi00, pi01, pi02, pi11, pi12, pi22, pi33, bulkPi;
   double particle_mu[Maxparticle];
//08/05/2013: added r and phi to correspond to xpt and ypt
}FO_surf;

typedef struct
{
  int monval;			// Montecarlo number according PDG
  string name;
  double mass;
  double width;
  int gspin;			// spin degeneracy
  int baryon;
  int strange;
  int charm;
  int bottom;
  int gisospin;			// isospin degeneracy
  int charge;
  int decays;			// amount of decays listed for this resonance
  int stable;			// defines whether this particle is considered as stable
  int decays_Npart[Maxdecaychannel];
  double decays_branchratio[Maxdecaychannel];
  int decays_part[Maxdecaychannel][Maxdecaypart];
  double mu;
  double thermal_yield;
  double percent_contribution;	//used to compute approximate percentage contribution of particle to net yield of some daughter particle
  double effective_branchratio;	//the effective branching ratio of generic resonance to some specified daughter particle
				//N.B. - effective branching ratio may be larger than 1
  double decays_effective_branchratio[Maxdecaychannel];
				//similar to effective_branchratio, but specific to each decay channel
  int sign;       //Bose-Einstein or Dirac-Fermi statistics
}particle_info;

void read_hydropar(hydropara* hp, string localpath);
int get_filelength(string filepath);
int get_filewidth(string filepath);
void read_decdat(int length, FO_surf* surf_ptr, string localpath, bool include_bulk_pi = false);
void read_surfdat(int length, FO_surf* surf_ptr, string localpath);
void read_decdat_mu(int FO_length, int N_stable, double** particle_mu, string localpath);
int read_resonance(particle_info* particle);
void calculate_particle_mu(int IEOS, int Nparticle, FO_surf* FOsurf_ptr, int FO_length, particle_info* particle, double** particle_mu);
void compute_total_contribution_percentages(int stable_particle_idx, int Nparticle, particle_info* particle);
void calculate_thermal_particle_yield(int Nparticle, particle_info* particle, double Temperature);
double b_j_to_i(particle_info * particle, int Nparticle, int j, int i, int verbose_monval = 0);
int count_targets(int * decay_channel_particles, particle_info * i);
int count_stable(particle_info * all_particles, int Nparticle, int * decay_channel_particles);
bool is_stable(particle_info * all_particles, int Nparticle, int monval);
int set_stable_particle_monval();
int lookup_particle_id_from_monval(particle_info * all_particles, int Nparticle, int monval);
void print_particle_stability(particle_info * all_particles, int Nparticle);
int get_number_of_decay_channels(vector<int> chosen_resonances, particle_info * particle);
void get_important_resonances(int chosen_target_particle_idx, vector<int> * chosen_resonance_indices_ptr, particle_info * particle, int Nparticle, double threshold, double &running_total_percentage, std::ofstream& output);
void get_all_descendants(vector<int> * chosen_resonance_indices_ptr, particle_info * particle, int Nparticle, std::ofstream& output);
void sort_by_mass(vector<int> * chosen_resonance_indices_ptr, particle_info * particle, int Nparticle, std::ofstream& output);



void calculate_pTdep_thermal_particle_yield(int Nparticle, particle_info* particle, double Temperature, double * pTarray, int npT, double ** resultsArray);
void compute_total_pTdep_contribution_percentages(int stable_particle_idx, int Nparticle, particle_info * particle, int npT,
						double ** pTdepThermalYields, vector<int> * chosen_resonance_indices_ptr, double * pTdep_Cumulative_Resonance_Contributions);


#endif
