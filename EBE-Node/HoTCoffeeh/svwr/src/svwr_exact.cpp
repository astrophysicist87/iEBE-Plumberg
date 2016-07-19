#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>
#include<time.h>

#include "svwr.h"
#include "Arsenal.h"
#include "Stopwatch.h"
#include "gauss_quadrature.h"

using namespace std;

double SourceVariances::Cal_wfi_dN_dypTdpTdphi_function(FO_surf* FOsurf_ptr, int local_pid, double pT, double pphi, int wfi)
{
	double result = 0.0;

	switch (wfi)
	{
		case 0:
			result = Cal_dN_dypTdpTdphi_function(FOsurf_ptr, local_pid, pT, pphi);
			break;
		case 1:
			result = Cal_xs_dN_dypTdpTdphi_function(FOsurf_ptr, local_pid, pT, pphi);
			break;
		case 2:
			result = Cal_xs2_dN_dypTdpTdphi_function(FOsurf_ptr, local_pid, pT, pphi);
			break;
		case 3:
			result = Cal_xo_dN_dypTdpTdphi_function(FOsurf_ptr, local_pid, pT, pphi);
			break;
		case 4:
			result = Cal_xo2_dN_dypTdpTdphi_function(FOsurf_ptr, local_pid, pT, pphi);
			break;
		case 5:
			//result = Cal_xl_dN_dypTdpTdphi_function(FOsurf_ptr, local_pid, pT, pphi);
			result = 0.0;
			break;
		case 6:
			result = Cal_xl2_dN_dypTdpTdphi_function(FOsurf_ptr, local_pid, pT, pphi);
			break;
		case 7:
			result = Cal_t_dN_dypTdpTdphi_function(FOsurf_ptr, local_pid, pT, pphi);
			break;
		case 8:
			result = Cal_t2_dN_dypTdpTdphi_function(FOsurf_ptr, local_pid, pT, pphi);
			break;
		case 9:
			result = Cal_xsxo_dN_dypTdpTdphi_function(FOsurf_ptr, local_pid, pT, pphi);
			break;
		case 10:
			//result = Cal_xsxl_dN_dypTdpTdphi_function(FOsurf_ptr, local_pid, pT, pphi);
			result = 0.0;
			break;
		case 11:
			result = Cal_xst_dN_dypTdpTdphi_function(FOsurf_ptr, local_pid, pT, pphi);
			break;
		case 12:
			//result = Cal_xoxl_dN_dypTdpTdphi_function(FOsurf_ptr, local_pid, pT, pphi);
			result = 0.0;
			break;
		case 13:
			result = Cal_xot_dN_dypTdpTdphi_function(FOsurf_ptr, local_pid, pT, pphi);
			break;
		case 14:
			//result = Cal_xlt_dN_dypTdpTdphi_function(FOsurf_ptr, local_pid, pT, pphi);
			result = 0.0;
			break;
		default:
			break;
	}

	return (result);
}

double SourceVariances::Cal_dN_dypTdpTdphi_function(FO_surf* FOsurf_ptr, int local_pid, double pT, double pphi)
{
	// set particle information
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;
	double deltaf_prefactor = 0.;
	if (INCLUDE_DELTA_F)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 2.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double dN_dypTdpTdphi = 0.0;

	for(int isurf=0; isurf<FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double r = surf->r;
		double sin_temp_phi = surf->sin_phi;
		double cos_temp_phi = surf->cos_phi;

		double vx = surf->vx;
		double vy = surf->vy;
		double gammaT = surf->gammaT;

		double da0 = surf->da0;
		double da1 = surf->da1;
		double da2 = surf->da2;

		double pi00 = surf->pi00;
		double pi01 = surf->pi01;
		double pi02 = surf->pi02;
		double pi11 = surf->pi11;
		double pi12 = surf->pi12;
		double pi22 = surf->pi22;
		double pi33 = surf->pi33;

		for(int ieta=0; ieta < eta_s_npts; ++ieta)
		{
			double p0 = sqrt(pT*pT+localmass*localmass)*cosh(SP_p_y - eta_s[ieta]);
			double pz = sqrt(pT*pT+localmass*localmass)*sinh(SP_p_y - eta_s[ieta]);

			//now get distribution function, emission function, etc.
			double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions

			//viscous corrections
			double deltaf = 0.;
			//if (INCLUDE_DELTA_F)
				deltaf = deltaf_prefactor * (1. - sign*f0)
							* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

			//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
			double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

			//ignore points where delta f is large or emission function goes negative from pdsigma
			if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
				S_p = 0.0;

			double S_p_withweight = S_p*tau*eta_s_weight[ieta];

			double weight_factor = eta_even_factor;

			dN_dypTdpTdphi += weight_factor * S_p_withweight;
		}
	}

	return dN_dypTdpTdphi;
}

double SourceVariances::Cal_xs_dN_dypTdpTdphi_function(FO_surf* FOsurf_ptr, int local_pid, double pT, double pphi)
{
	// set particle information
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;
	double deltaf_prefactor = 0.;
	if (INCLUDE_DELTA_F)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 2.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double dN_dypTdpTdphi = 0.0;

	for(int isurf=0; isurf<FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double r = surf->r;
		double sin_temp_phi = surf->sin_phi;
		double cos_temp_phi = surf->cos_phi;

		double vx = surf->vx;
		double vy = surf->vy;
		double gammaT = surf->gammaT;

		double da0 = surf->da0;
		double da1 = surf->da1;
		double da2 = surf->da2;

		double pi00 = surf->pi00;
		double pi01 = surf->pi01;
		double pi02 = surf->pi02;
		double pi11 = surf->pi11;
		double pi12 = surf->pi12;
		double pi22 = surf->pi22;
		double pi33 = surf->pi33;

		double z2 = r * (sin_temp_phi * cos_pphi - cos_temp_phi * sin_pphi);

		for(int ieta=0; ieta < eta_s_npts; ++ieta)
		{
			double p0 = sqrt(pT*pT+localmass*localmass)*cosh(SP_p_y - eta_s[ieta]);
			double pz = sqrt(pT*pT+localmass*localmass)*sinh(SP_p_y - eta_s[ieta]);

			//now get distribution function, emission function, etc.
			double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
	
			//viscous corrections
			double deltaf = 0.;
			//if (INCLUDE_DELTA_F)
				deltaf = deltaf_prefactor * (1. - sign*f0)
							* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

			//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
			double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

			//ignore points where delta f is large or emission function goes negative from pdsigma
			if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
				S_p = 0.0;

			double S_p_withweight = S_p*tau*eta_s_weight[ieta];

			double weight_factor = eta_even_factor*z2;

			dN_dypTdpTdphi += weight_factor * S_p_withweight;
		}
	}

	return dN_dypTdpTdphi;
}

double SourceVariances::Cal_xs2_dN_dypTdpTdphi_function(FO_surf* FOsurf_ptr, int local_pid, double pT, double pphi)
{
	// set particle information
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;
	double deltaf_prefactor = 0.;
	if (INCLUDE_DELTA_F)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 2.0;
	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double dN_dypTdpTdphi = 0.0;

	for(int isurf=0; isurf<FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double r = surf->r;
		double sin_temp_phi = surf->sin_phi;
		double cos_temp_phi = surf->cos_phi;

		double vx = surf->vx;
		double vy = surf->vy;
		double gammaT = surf->gammaT;

		double da0 = surf->da0;
		double da1 = surf->da1;
		double da2 = surf->da2;

		double pi00 = surf->pi00;
		double pi01 = surf->pi01;
		double pi02 = surf->pi02;
		double pi11 = surf->pi11;
		double pi12 = surf->pi12;
		double pi22 = surf->pi22;
		double pi33 = surf->pi33;

		double z2 = r * (sin_temp_phi * cos_pphi - cos_temp_phi * sin_pphi);

		for(int ieta=0; ieta < eta_s_npts; ++ieta)
		{
			double p0 = sqrt(pT*pT+localmass*localmass)*cosh(SP_p_y - eta_s[ieta]);
			double pz = sqrt(pT*pT+localmass*localmass)*sinh(SP_p_y - eta_s[ieta]);

			//now get distribution function, emission function, etc.
			double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
	
			//viscous corrections
			double deltaf = 0.;
			//if (INCLUDE_DELTA_F)
				deltaf = deltaf_prefactor * (1. - sign*f0)
							* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

			//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
			double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

			//ignore points where delta f is large or emission function goes negative from pdsigma
			if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
				S_p = 0.0;

			double S_p_withweight = S_p*tau*eta_s_weight[ieta];

			double weight_factor = eta_even_factor*z2*z2;	//2

			dN_dypTdpTdphi += weight_factor * S_p_withweight;
		}
	}

	return dN_dypTdpTdphi;
}

double SourceVariances::Cal_xo_dN_dypTdpTdphi_function(FO_surf* FOsurf_ptr, int local_pid, double pT, double pphi)
{
	// set particle information
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;
	double deltaf_prefactor = 0.;
	if (INCLUDE_DELTA_F)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 2.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double dN_dypTdpTdphi = 0.0;

	for(int isurf=0; isurf<FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double r = surf->r;
		double sin_temp_phi = surf->sin_phi;
		double cos_temp_phi = surf->cos_phi;

		double vx = surf->vx;
		double vy = surf->vy;
		double gammaT = surf->gammaT;

		double da0 = surf->da0;
		double da1 = surf->da1;
		double da2 = surf->da2;

		double pi00 = surf->pi00;
		double pi01 = surf->pi01;
		double pi02 = surf->pi02;
		double pi11 = surf->pi11;
		double pi12 = surf->pi12;
		double pi22 = surf->pi22;
		double pi33 = surf->pi33;

		double z1 = r * (cos_temp_phi * cos_pphi + sin_temp_phi * sin_pphi);

		for(int ieta=0; ieta < eta_s_npts; ++ieta)
		{
			double p0 = sqrt(pT*pT+localmass*localmass)*cosh(SP_p_y - eta_s[ieta]);
			double pz = sqrt(pT*pT+localmass*localmass)*sinh(SP_p_y - eta_s[ieta]);

			//now get distribution function, emission function, etc.
			double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
	
			//viscous corrections
			double deltaf = 0.;
			//if (INCLUDE_DELTA_F)
				deltaf = deltaf_prefactor * (1. - sign*f0)
							* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

			//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
			double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

			//ignore points where delta f is large or emission function goes negative from pdsigma
			if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
				S_p = 0.0;

			double S_p_withweight = S_p*tau*eta_s_weight[ieta];

			double weight_factor = eta_even_factor*z1;	//3

			dN_dypTdpTdphi += weight_factor * S_p_withweight;
		}
	}

	return dN_dypTdpTdphi;
}

double SourceVariances::Cal_xo2_dN_dypTdpTdphi_function(FO_surf* FOsurf_ptr, int local_pid, double pT, double pphi)
{
	// set particle information
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;
	double deltaf_prefactor = 0.;
	if (INCLUDE_DELTA_F)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 2.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double dN_dypTdpTdphi = 0.0;

	for(int isurf=0; isurf<FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double r = surf->r;
		double sin_temp_phi = surf->sin_phi;
		double cos_temp_phi = surf->cos_phi;

		double vx = surf->vx;
		double vy = surf->vy;
		double gammaT = surf->gammaT;

		double da0 = surf->da0;
		double da1 = surf->da1;
		double da2 = surf->da2;

		double pi00 = surf->pi00;
		double pi01 = surf->pi01;
		double pi02 = surf->pi02;
		double pi11 = surf->pi11;
		double pi12 = surf->pi12;
		double pi22 = surf->pi22;
		double pi33 = surf->pi33;

		double z1 = r * (cos_temp_phi * cos_pphi + sin_temp_phi * sin_pphi);

		for(int ieta=0; ieta < eta_s_npts; ++ieta)
		{
			double p0 = sqrt(pT*pT+localmass*localmass)*cosh(SP_p_y - eta_s[ieta]);
			double pz = sqrt(pT*pT+localmass*localmass)*sinh(SP_p_y - eta_s[ieta]);

			//now get distribution function, emission function, etc.
			double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
	
			//viscous corrections
			double deltaf = 0.;
			//if (INCLUDE_DELTA_F)
				deltaf = deltaf_prefactor * (1. - sign*f0)
							* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

			//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
			double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

			//ignore points where delta f is large or emission function goes negative from pdsigma
			if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
				S_p = 0.0;

			double S_p_withweight = S_p*tau*eta_s_weight[ieta];

			double weight_factor = eta_even_factor*z1*z1;	//4

			dN_dypTdpTdphi += weight_factor * S_p_withweight;
		}
	}

	return dN_dypTdpTdphi;
}

double SourceVariances::Cal_xl_dN_dypTdpTdphi_function(FO_surf* FOsurf_ptr, int local_pid, double pT, double pphi)
{
	// set particle information
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;
	double deltaf_prefactor = 0.;
	if (INCLUDE_DELTA_F)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 2.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double dN_dypTdpTdphi = 0.0;

	for(int isurf=0; isurf<FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double r = surf->r;
		double sin_temp_phi = surf->sin_phi;
		double cos_temp_phi = surf->cos_phi;

		double vx = surf->vx;
		double vy = surf->vy;
		double gammaT = surf->gammaT;

		double da0 = surf->da0;
		double da1 = surf->da1;
		double da2 = surf->da2;

		double pi00 = surf->pi00;
		double pi01 = surf->pi01;
		double pi02 = surf->pi02;
		double pi11 = surf->pi11;
		double pi12 = surf->pi12;
		double pi22 = surf->pi22;
		double pi33 = surf->pi33;

		for(int ieta=0; ieta < eta_s_npts; ++ieta)
		{
			double p0 = sqrt(pT*pT+localmass*localmass)*cosh(SP_p_y - eta_s[ieta]);
			double pz = sqrt(pT*pT+localmass*localmass)*sinh(SP_p_y - eta_s[ieta]);

			//now get distribution function, emission function, etc.
			double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
	
			//viscous corrections
			double deltaf = 0.;
			//if (INCLUDE_DELTA_F)
				deltaf = deltaf_prefactor * (1. - sign*f0)
							* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

			//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
			double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

			//ignore points where delta f is large or emission function goes negative from pdsigma
			if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
				S_p = 0.0;

			double S_p_withweight = S_p*tau*eta_s_weight[ieta];
			double z3 = tau*sh_eta_s[ieta];

			double weight_factor = 0.0;	//5

			dN_dypTdpTdphi += weight_factor * S_p_withweight;
		}
	}

	return dN_dypTdpTdphi;
}

double SourceVariances::Cal_xl2_dN_dypTdpTdphi_function(FO_surf* FOsurf_ptr, int local_pid, double pT, double pphi)
{
	// set particle information
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;
	double deltaf_prefactor = 0.;
	if (INCLUDE_DELTA_F)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 2.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double dN_dypTdpTdphi = 0.0;

	for(int isurf=0; isurf<FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double r = surf->r;
		double sin_temp_phi = surf->sin_phi;
		double cos_temp_phi = surf->cos_phi;

		double vx = surf->vx;
		double vy = surf->vy;
		double gammaT = surf->gammaT;

		double da0 = surf->da0;
		double da1 = surf->da1;
		double da2 = surf->da2;

		double pi00 = surf->pi00;
		double pi01 = surf->pi01;
		double pi02 = surf->pi02;
		double pi11 = surf->pi11;
		double pi12 = surf->pi12;
		double pi22 = surf->pi22;
		double pi33 = surf->pi33;

		for(int ieta=0; ieta < eta_s_npts; ++ieta)
		{
			double p0 = sqrt(pT*pT+localmass*localmass)*cosh(SP_p_y - eta_s[ieta]);
			double pz = sqrt(pT*pT+localmass*localmass)*sinh(SP_p_y - eta_s[ieta]);

			//now get distribution function, emission function, etc.
			double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
	
			//viscous corrections
			double deltaf = 0.;
			//if (INCLUDE_DELTA_F)
				deltaf = deltaf_prefactor * (1. - sign*f0)
							* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

			//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
			double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

			//ignore points where delta f is large or emission function goes negative from pdsigma
			if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
				S_p = 0.0;

			double S_p_withweight = S_p*tau*eta_s_weight[ieta];
			double z3 = tau*sh_eta_s[ieta];

			double weight_factor = eta_even_factor*z3*z3;	//6

			dN_dypTdpTdphi += weight_factor * S_p_withweight;
		}
	}

	return dN_dypTdpTdphi;
}

double SourceVariances::Cal_t_dN_dypTdpTdphi_function(FO_surf* FOsurf_ptr, int local_pid, double pT, double pphi)
{
	// set particle information
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;
	double deltaf_prefactor = 0.;
	if (INCLUDE_DELTA_F)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 2.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double dN_dypTdpTdphi = 0.0;

	for(int isurf=0; isurf<FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double r = surf->r;
		double sin_temp_phi = surf->sin_phi;
		double cos_temp_phi = surf->cos_phi;

		double vx = surf->vx;
		double vy = surf->vy;
		double gammaT = surf->gammaT;

		double da0 = surf->da0;
		double da1 = surf->da1;
		double da2 = surf->da2;

		double pi00 = surf->pi00;
		double pi01 = surf->pi01;
		double pi02 = surf->pi02;
		double pi11 = surf->pi11;
		double pi12 = surf->pi12;
		double pi22 = surf->pi22;
		double pi33 = surf->pi33;

		for(int ieta=0; ieta < eta_s_npts; ++ieta)
		{
			double p0 = sqrt(pT*pT+localmass*localmass)*cosh(SP_p_y - eta_s[ieta]);
			double pz = sqrt(pT*pT+localmass*localmass)*sinh(SP_p_y - eta_s[ieta]);

			//now get distribution function, emission function, etc.
			double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
	
			//viscous corrections
			double deltaf = 0.;
			//if (INCLUDE_DELTA_F)
				deltaf = deltaf_prefactor * (1. - sign*f0)
							* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

			//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
			double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

			//ignore points where delta f is large or emission function goes negative from pdsigma
			if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
				S_p = 0.0;

			double S_p_withweight = S_p*tau*eta_s_weight[ieta];
			double z0 = tau*ch_eta_s[ieta];

			double weight_factor = eta_even_factor*z0;	//7

			dN_dypTdpTdphi += weight_factor * S_p_withweight;
		}
	}

	return dN_dypTdpTdphi;
}

double SourceVariances::Cal_t2_dN_dypTdpTdphi_function(FO_surf* FOsurf_ptr, int local_pid, double pT, double pphi)
{
	// set particle information
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;
	double deltaf_prefactor = 0.;
	if (INCLUDE_DELTA_F)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 2.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double dN_dypTdpTdphi = 0.0;

	for(int isurf=0; isurf<FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double r = surf->r;
		double sin_temp_phi = surf->sin_phi;
		double cos_temp_phi = surf->cos_phi;

		double vx = surf->vx;
		double vy = surf->vy;
		double gammaT = surf->gammaT;

		double da0 = surf->da0;
		double da1 = surf->da1;
		double da2 = surf->da2;

		double pi00 = surf->pi00;
		double pi01 = surf->pi01;
		double pi02 = surf->pi02;
		double pi11 = surf->pi11;
		double pi12 = surf->pi12;
		double pi22 = surf->pi22;
		double pi33 = surf->pi33;

		for(int ieta=0; ieta < eta_s_npts; ++ieta)
		{
			double p0 = sqrt(pT*pT+localmass*localmass)*cosh(SP_p_y - eta_s[ieta]);
			double pz = sqrt(pT*pT+localmass*localmass)*sinh(SP_p_y - eta_s[ieta]);

			//now get distribution function, emission function, etc.
			double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
	
			//viscous corrections
			double deltaf = 0.;
			//if (INCLUDE_DELTA_F)
				deltaf = deltaf_prefactor * (1. - sign*f0)
							* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

			//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
			double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

			//ignore points where delta f is large or emission function goes negative from pdsigma
			if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
				S_p = 0.0;

			double S_p_withweight = S_p*tau*eta_s_weight[ieta];
			double z0 = tau*ch_eta_s[ieta];

			double weight_factor = eta_even_factor*z0*z0;	//8

			dN_dypTdpTdphi += weight_factor * S_p_withweight;
		}
	}

	return dN_dypTdpTdphi;
}

double SourceVariances::Cal_xsxo_dN_dypTdpTdphi_function(FO_surf* FOsurf_ptr, int local_pid, double pT, double pphi)
{
	// set particle information
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;
	double deltaf_prefactor = 0.;
	if (INCLUDE_DELTA_F)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 2.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double dN_dypTdpTdphi = 0.0;

	for(int isurf=0; isurf<FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double r = surf->r;
		double sin_temp_phi = surf->sin_phi;
		double cos_temp_phi = surf->cos_phi;

		double vx = surf->vx;
		double vy = surf->vy;
		double gammaT = surf->gammaT;

		double da0 = surf->da0;
		double da1 = surf->da1;
		double da2 = surf->da2;

		double pi00 = surf->pi00;
		double pi01 = surf->pi01;
		double pi02 = surf->pi02;
		double pi11 = surf->pi11;
		double pi12 = surf->pi12;
		double pi22 = surf->pi22;
		double pi33 = surf->pi33;

		double z1 = r * (cos_temp_phi * cos_pphi + sin_temp_phi * sin_pphi);
		double z2 = r * (sin_temp_phi * cos_pphi - cos_temp_phi * sin_pphi);

		for(int ieta=0; ieta < eta_s_npts; ++ieta)
		{
			double p0 = sqrt(pT*pT+localmass*localmass)*cosh(SP_p_y - eta_s[ieta]);
			double pz = sqrt(pT*pT+localmass*localmass)*sinh(SP_p_y - eta_s[ieta]);

			//now get distribution function, emission function, etc.
			double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
	
			//viscous corrections
			double deltaf = 0.;
			//if (INCLUDE_DELTA_F)
				deltaf = deltaf_prefactor * (1. - sign*f0)
							* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

			//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
			double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

			//ignore points where delta f is large or emission function goes negative from pdsigma
			if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
				S_p = 0.0;

			double S_p_withweight = S_p*tau*eta_s_weight[ieta];

			double weight_factor = eta_even_factor*z2*z1;	//9

			dN_dypTdpTdphi += weight_factor * S_p_withweight;
		}
	}

	return dN_dypTdpTdphi;
}

double SourceVariances::Cal_xsxl_dN_dypTdpTdphi_function(FO_surf* FOsurf_ptr, int local_pid, double pT, double pphi)
{
	// set particle information
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;
	double deltaf_prefactor = 0.;
	if (INCLUDE_DELTA_F)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 2.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double dN_dypTdpTdphi = 0.0;

	for(int isurf=0; isurf<FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double r = surf->r;
		double sin_temp_phi = surf->sin_phi;
		double cos_temp_phi = surf->cos_phi;

		double vx = surf->vx;
		double vy = surf->vy;
		double gammaT = surf->gammaT;

		double da0 = surf->da0;
		double da1 = surf->da1;
		double da2 = surf->da2;

		double pi00 = surf->pi00;
		double pi01 = surf->pi01;
		double pi02 = surf->pi02;
		double pi11 = surf->pi11;
		double pi12 = surf->pi12;
		double pi22 = surf->pi22;
		double pi33 = surf->pi33;

		double z2 = r * (sin_temp_phi * cos_pphi - cos_temp_phi * sin_pphi);

		for(int ieta=0; ieta < eta_s_npts; ++ieta)
		{
			double p0 = sqrt(pT*pT+localmass*localmass)*cosh(SP_p_y - eta_s[ieta]);
			double pz = sqrt(pT*pT+localmass*localmass)*sinh(SP_p_y - eta_s[ieta]);

			//now get distribution function, emission function, etc.
			double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
	
			//viscous corrections
			double deltaf = 0.;
			//if (INCLUDE_DELTA_F)
				deltaf = deltaf_prefactor * (1. - sign*f0)
							* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

			//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
			double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

			//ignore points where delta f is large or emission function goes negative from pdsigma
			if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
				S_p = 0.0;

			double S_p_withweight = S_p*tau*eta_s_weight[ieta];
		//double z3 = tau*sh_eta_s[ieta];

			double weight_factor = 0.0;

			dN_dypTdpTdphi += weight_factor * S_p_withweight;
		}
	}

	return dN_dypTdpTdphi;
}

double SourceVariances::Cal_xst_dN_dypTdpTdphi_function(FO_surf* FOsurf_ptr, int local_pid, double pT, double pphi)
{
	// set particle information
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;
	double deltaf_prefactor = 0.;
	if (INCLUDE_DELTA_F)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 2.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double dN_dypTdpTdphi = 0.0;

	for(int isurf=0; isurf<FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double r = surf->r;
		double sin_temp_phi = surf->sin_phi;
		double cos_temp_phi = surf->cos_phi;

		double vx = surf->vx;
		double vy = surf->vy;
		double gammaT = surf->gammaT;

		double da0 = surf->da0;
		double da1 = surf->da1;
		double da2 = surf->da2;

		double pi00 = surf->pi00;
		double pi01 = surf->pi01;
		double pi02 = surf->pi02;
		double pi11 = surf->pi11;
		double pi12 = surf->pi12;
		double pi22 = surf->pi22;
		double pi33 = surf->pi33;

		double z2 = r * (sin_temp_phi * cos_pphi - cos_temp_phi * sin_pphi);

		for(int ieta=0; ieta < eta_s_npts; ++ieta)
		{
			double p0 = sqrt(pT*pT+localmass*localmass)*cosh(SP_p_y - eta_s[ieta]);
			double pz = sqrt(pT*pT+localmass*localmass)*sinh(SP_p_y - eta_s[ieta]);

			//now get distribution function, emission function, etc.
			double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
	
			//viscous corrections
			double deltaf = 0.;
			//if (INCLUDE_DELTA_F)
				deltaf = deltaf_prefactor * (1. - sign*f0)
							* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

			//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
			double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

			//ignore points where delta f is large or emission function goes negative from pdsigma
			if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
				S_p = 0.0;

			double S_p_withweight = S_p*tau*eta_s_weight[ieta];
			double z0 = tau*ch_eta_s[ieta];

			double weight_factor = eta_even_factor*z2*z0;	//11

			dN_dypTdpTdphi += weight_factor * S_p_withweight;
		}
	}

	return dN_dypTdpTdphi;
}

double SourceVariances::Cal_xoxl_dN_dypTdpTdphi_function(FO_surf* FOsurf_ptr, int local_pid, double pT, double pphi)
{
	// set particle information
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;
	double deltaf_prefactor = 0.;
	if (INCLUDE_DELTA_F)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 2.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double dN_dypTdpTdphi = 0.0;

	for(int isurf=0; isurf<FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double r = surf->r;
		double sin_temp_phi = surf->sin_phi;
		double cos_temp_phi = surf->cos_phi;

		double vx = surf->vx;
		double vy = surf->vy;
		double gammaT = surf->gammaT;

		double da0 = surf->da0;
		double da1 = surf->da1;
		double da2 = surf->da2;

		double pi00 = surf->pi00;
		double pi01 = surf->pi01;
		double pi02 = surf->pi02;
		double pi11 = surf->pi11;
		double pi12 = surf->pi12;
		double pi22 = surf->pi22;
		double pi33 = surf->pi33;

		for(int ieta=0; ieta < eta_s_npts; ++ieta)
		{
			double p0 = sqrt(pT*pT+localmass*localmass)*cosh(SP_p_y - eta_s[ieta]);
			double pz = sqrt(pT*pT+localmass*localmass)*sinh(SP_p_y - eta_s[ieta]);

			//now get distribution function, emission function, etc.
			double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
	
			//viscous corrections
			double deltaf = 0.;
			//if (INCLUDE_DELTA_F)
				deltaf = deltaf_prefactor * (1. - sign*f0)
							* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

			//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
			double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

			//ignore points where delta f is large or emission function goes negative from pdsigma
			if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
				S_p = 0.0;

			double S_p_withweight = S_p*tau*eta_s_weight[ieta];
		//double z1 = r * (cos_temp_phi * cos_pphi + sin_temp_phi * sin_pphi);
		//double z3 = tau*sh_eta_s[ieta];


			double weight_factor = 0.0;

			dN_dypTdpTdphi += weight_factor * S_p_withweight;
		}
	}

	return dN_dypTdpTdphi;
}

double SourceVariances::Cal_xot_dN_dypTdpTdphi_function(FO_surf* FOsurf_ptr, int local_pid, double pT, double pphi)
{
	// set particle information
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;
	double deltaf_prefactor = 0.;
	if (INCLUDE_DELTA_F)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 2.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double dN_dypTdpTdphi = 0.0;

	for(int isurf=0; isurf<FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double r = surf->r;
		double sin_temp_phi = surf->sin_phi;
		double cos_temp_phi = surf->cos_phi;

		double vx = surf->vx;
		double vy = surf->vy;
		double gammaT = surf->gammaT;

		double da0 = surf->da0;
		double da1 = surf->da1;
		double da2 = surf->da2;

		double pi00 = surf->pi00;
		double pi01 = surf->pi01;
		double pi02 = surf->pi02;
		double pi11 = surf->pi11;
		double pi12 = surf->pi12;
		double pi22 = surf->pi22;
		double pi33 = surf->pi33;

		double z1 = r * (cos_temp_phi * cos_pphi + sin_temp_phi * sin_pphi);

		for(int ieta=0; ieta < eta_s_npts; ++ieta)
		{
			double p0 = sqrt(pT*pT+localmass*localmass)*cosh(SP_p_y - eta_s[ieta]);
			double pz = sqrt(pT*pT+localmass*localmass)*sinh(SP_p_y - eta_s[ieta]);

			//now get distribution function, emission function, etc.
			double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
	
			//viscous corrections
			double deltaf = 0.;
			//if (INCLUDE_DELTA_F)
				deltaf = deltaf_prefactor * (1. - sign*f0)
							* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

			//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
			double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

			//ignore points where delta f is large or emission function goes negative from pdsigma
			if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
				S_p = 0.0;

			double S_p_withweight = S_p*tau*eta_s_weight[ieta];
			double z0 = tau*ch_eta_s[ieta];

			double weight_factor = eta_even_factor*z1*z0;

			dN_dypTdpTdphi += weight_factor * S_p_withweight;
		}
	}

	return dN_dypTdpTdphi;
}

double SourceVariances::Cal_xlt_dN_dypTdpTdphi_function(FO_surf* FOsurf_ptr, int local_pid, double pT, double pphi)
{
	// set particle information
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;
	double deltaf_prefactor = 0.;
	if (INCLUDE_DELTA_F)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 2.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double dN_dypTdpTdphi = 0.0;

	for(int isurf=0; isurf<FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double r = surf->r;
		double sin_temp_phi = surf->sin_phi;
		double cos_temp_phi = surf->cos_phi;

		double vx = surf->vx;
		double vy = surf->vy;
		double gammaT = surf->gammaT;

		double da0 = surf->da0;
		double da1 = surf->da1;
		double da2 = surf->da2;

		double pi00 = surf->pi00;
		double pi01 = surf->pi01;
		double pi02 = surf->pi02;
		double pi11 = surf->pi11;
		double pi12 = surf->pi12;
		double pi22 = surf->pi22;
		double pi33 = surf->pi33;

		for(int ieta=0; ieta < eta_s_npts; ++ieta)
		{
			double p0 = sqrt(pT*pT+localmass*localmass)*cosh(SP_p_y - eta_s[ieta]);
			double pz = sqrt(pT*pT+localmass*localmass)*sinh(SP_p_y - eta_s[ieta]);

			//now get distribution function, emission function, etc.
			double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
	
			//viscous corrections
			double deltaf = 0.;
			//if (INCLUDE_DELTA_F)
				deltaf = deltaf_prefactor * (1. - sign*f0)
							* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

			//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
			double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

			//ignore points where delta f is large or emission function goes negative from pdsigma
			if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
				S_p = 0.0;

			double S_p_withweight = S_p*tau*eta_s_weight[ieta];
		//double z0 = tau*ch_eta_s[ieta];
		//double z3 = tau*sh_eta_s[ieta];

			double weight_factor = 0.0;

			dN_dypTdpTdphi += weight_factor * S_p_withweight;
		}
	}

	return dN_dypTdpTdphi;
}
