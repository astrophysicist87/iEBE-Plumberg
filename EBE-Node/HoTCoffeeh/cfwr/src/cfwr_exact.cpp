#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<set>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<queue>
#include<map>
#include<cstdlib>
#include<numeric>

#include "cfwr.h"
#include "gauss_quadrature.h"
#include "bessel.h"

using namespace std;

const std::complex<double> i(0, 1);

void CorrelationFunction::Cal_dN_dypTdpTdphi_with_rapdep_thermal(int local_pid)
{
	n_pY_pts = 101;
	SP_Del_pY = new double [n_pY_pts];
	double * dummy = new double [n_pY_pts];
	gauss_quadrature(n_pY_pts, 1, 0.0, 0.0, -4.0, 4.0, SP_Del_pY, dummy);

	for (int ipt = 0; ipt < n_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	{
			*global_out_stream_ptr << "Made it to " << ipt << "   " << ipphi << endl;

			for (int iqt = 0; iqt < qtnpts; ++iqt)
			for (int iqx = 0; iqx < qxnpts; ++iqx)
			for (int iqy = 0; iqy < qynpts; ++iqy)
			for (int iqz = 0; iqz < qznpts; ++iqz)
			{
				//if ( (iqt != (qtnpts-1)/2) || (iqx != (qxnpts-1)/2) || (iqy != (qynpts-1)/2) || (iqz != (qznpts-1)/2 ) )
				//	continue;

				Cal_dN_dypTdpTdphi_with_weights_function_approx(local_pid, SP_pT[ipt], SP_pphi[ipphi], SP_Del_pY,
																qt_pts[iqt], qx_pts[iqx], qy_pts[iqy], qz_pts[iqz], cos_vs_pY, sin_vs_pY);

				for (int ipY = 0; ipY < n_pY_pts; ++ipY)
				{
					//if (ipY != 1) continue;
					double tmp_cos = 0.0, tmp_sin = 0.0, tmp_cos_BA = 0.0, tmp_sin_BA = 0.0;
					Cal_dN_dypTdpTdphi_with_weights_function(local_pid, SP_pT[ipt], SP_pphi[ipphi], SP_Del_pY[ipY],
																qt_pts[iqt], qx_pts[iqx], qy_pts[iqy], qz_pts[iqz], &tmp_cos, &tmp_sin);
					Cal_dN_dypTdpTdphi_with_weights_function_BoltzmannApproximation(local_pid, SP_pT[ipt], SP_pphi[ipphi], SP_Del_pY[ipY],
																qt_pts[iqt], qx_pts[iqx], qy_pts[iqy], qz_pts[iqz], &tmp_cos_BA, &tmp_sin_BA);
					printf("FT THERMAL RHOS: %15.8le  %15.8le  %15.8le  %15.8le  %15.8le  %15.8le  %15.8le  %15.8le  %15.8le  %15.8le  %15.8le  %15.8le  %15.8le\n",
							SP_pT[ipt], SP_pphi[ipphi], SP_Del_pY[ipY], qt_pts[iqt], qx_pts[iqx], qy_pts[iqy], qz_pts[iqz],
							tmp_cos, tmp_sin, tmp_cos_BA, tmp_sin_BA, cos_vs_pY[ipY], sin_vs_pY[ipY]);
				}
		}
		if (1) exit(0);
	}

	if (1) exit(0);
	return;
}

void CorrelationFunction::Cal_dN_dypTdpTdphi_with_weights_function(int local_pid, double pT, double pphi, double p_y,
					double qt, double qx, double qy, double qz, double * cosqx_dN_dypTdpTdphi, double * sinqx_dN_dypTdpTdphi)
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
	if (use_delta_f)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));
	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 2.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;
	double ch_pY = cosh(p_y);
	double sh_pY = sinh(p_y);

	double mT = sqrt(pT*pT+localmass*localmass);

	*cosqx_dN_dypTdpTdphi = 0.0;
	*sinqx_dN_dypTdpTdphi = 0.0;

	//double local_tmp_cos = 0.0, local_tmp_sin = 0.0;
	for (int isurf = 0; isurf < FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double xpt = surf->xpt;
		double ypt = surf->ypt;

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

		for (int ieta = 0; ieta < eta_s_npts; ++ieta)
		{
			double p0 = mT*ch_Delta_eta_s[ieta];
			double pz = mT*sh_Delta_eta_s[ieta];

			double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
			//viscous corrections
			double deltaf = 0.;
			if (use_delta_f)
				deltaf = deltaf_prefactor * (1. - sign*f0)
							* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

			//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
			double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

			//ignore points where delta f is large or emission function goes negative from pdsigma
			/*if ( (1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol) )
			{
				S_p = 0.0;
				continue;
			}*/

			double tpt = tau*(ch_pY*ch_Delta_eta_s[ieta] - sh_pY*sh_Delta_eta_s[ieta]);
			double zpt = tau*(sh_pY*ch_Delta_eta_s[ieta] - ch_pY*sh_Delta_eta_s[ieta]);
			double arg = tpt*qt-(xpt*qx+ypt*qy+zpt*qz);
			*cosqx_dN_dypTdpTdphi += cos(arg/hbarC)*S_p*tau*Delta_eta_s_weight[ieta];
			*sinqx_dN_dypTdpTdphi += sin(arg/hbarC)*S_p*tau*Delta_eta_s_weight[ieta];
			//local_tmp_cos += cos(arg/hbarC)*S_p*tau*Delta_eta_s_weight[ieta];
			//local_tmp_sin += sin(arg/hbarC)*S_p*tau*Delta_eta_s_weight[ieta];
		}
			//printf("RAPIDITY DEPENDENCE: %15.8le  %15.8le  %15.8le  %15.8le  %15.8le  %15.8le  %15.8le  %15.8le  %15.8le  %15.8le\n",
			//		pT, pphi, p_y, qt, qx, qy, qz, Delta_eta_s[ieta], local_tmp_cos, local_tmp_sin);
	}

	return;
}

void CorrelationFunction::Cal_dN_dypTdpTdphi_with_weights_function_BoltzmannApproximation(int local_pid, double pT, double pphi, double p_y,
					double qt, double qx, double qy, double qz, double * cosqx_dN_dypTdpTdphi, double * sinqx_dN_dypTdpTdphi)
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
	if (use_delta_f)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));
	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 2.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;
	double ch_pY = cosh(p_y);
	double sh_pY = sinh(p_y);

	double mT = sqrt(pT*pT+localmass*localmass);

	*cosqx_dN_dypTdpTdphi = 0.0;
	*sinqx_dN_dypTdpTdphi = 0.0;

	//double local_tmp_cos = 0.0, local_tmp_sin = 0.0;
	for (int isurf = 0; isurf < FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double xpt = surf->xpt;
		double ypt = surf->ypt;

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

		for (int ieta = 0; ieta < eta_s_npts; ++ieta)
		{
			double p0 = mT*ch_Delta_eta_s[ieta];
			double pz = mT*sh_Delta_eta_s[ieta];

			//double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
			double f0 = exp( - one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) );	//thermal equilibrium distributions
			//viscous corrections
			double deltaf = 0.;
			if (use_delta_f)
				deltaf = deltaf_prefactor * (1. - sign*f0)
							* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

			//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
			double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

			double tpt = tau*(ch_pY*ch_Delta_eta_s[ieta] - sh_pY*sh_Delta_eta_s[ieta]);
			double zpt = tau*(sh_pY*ch_Delta_eta_s[ieta] - ch_pY*sh_Delta_eta_s[ieta]);
			double arg = tpt*qt-(xpt*qx+ypt*qy+zpt*qz);
			*cosqx_dN_dypTdpTdphi += cos(arg/hbarC)*S_p*tau*Delta_eta_s_weight[ieta];
			*sinqx_dN_dypTdpTdphi += sin(arg/hbarC)*S_p*tau*Delta_eta_s_weight[ieta];
		}
	}

	return;
}

void CorrelationFunction::Cal_dN_dypTdpTdphi_with_weights_function_approx(int local_pid, double pT, double pphi, double * p_Y_pts,
					double qt, double qx, double qy, double qz,
					double * cosqx_dN_dypTdpTdphi, double * sinqx_dN_dypTdpTdphi)
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
	if (use_delta_f)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));
	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 0.0, eta_even_factor = 2.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double mT = sqrt(pT*pT+localmass*localmass);
	//double tmp_cos_ta_num

	for (int isurf = 0; isurf < FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double xpt = surf->xpt;
		double ypt = surf->ypt;

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

		double transverse_arg = -(xpt*qx+ypt*qy)/hbarC;
		double cos_ta = cos(transverse_arg);
		double sin_ta = sin(transverse_arg);

		double A = tau*prefactor*mT*da0;
		double B = tau*prefactor*(px*da1 + py*da2);
		double C = deltaf_prefactor;
		
		double a = mT*mT*(pi00 + pi33);
		double b = -2.0*mT*(px*pi01 + py*pi02);
		double c = px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 - mT*mT*pi33;
		
		double alpha = one_by_Tdec*gammaT*mT;
		double transverse_f0 = exp( one_by_Tdec*(gammaT*(px*vx + py*vy) + mu) );

		for (int ipY = 0; ipY < n_pY_pts; ++ipY)
		{
			double local_pY = SP_Del_pY[ipY];
			double ch_pY = cosh(local_pY);
			double sh_pY = sinh(local_pY);
			double beta = (tau / hbarC) * ( qt*ch_pY - qz*sh_pY );
			double gamma = (tau / hbarC) * ( qz*ch_pY - qt*sh_pY );

			complex<double> I0_a_b_g, I1_a_b_g, I2_a_b_g, I3_a_b_g;
			complex<double> I0_2a_b_g, I1_2a_b_g, I2_2a_b_g, I3_2a_b_g;
			I(alpha, beta, gamma, I0_a_b_g, I1_a_b_g, I2_a_b_g, I3_a_b_g);
			I(2.0*alpha, beta, gamma, I0_2a_b_g, I1_2a_b_g, I2_2a_b_g, I3_2a_b_g);

			complex<double> term1 = transverse_f0 * (A*I1_a_b_g + B*I0_a_b_g);
			complex<double> term2 = C * transverse_f0 * ( A*a*I3_a_b_g + (B*a+b*A)*I2_a_b_g + (B*b+c*A)*I1_a_b_g + B*c*I0_a_b_g );
			complex<double> term3 = -sign * C * transverse_f0 * transverse_f0 * ( A*a*I3_2a_b_g + (B*a+b*A)*I2_2a_b_g + (B*b+c*A)*I1_2a_b_g + B*c*I0_2a_b_g );

			complex<double> eiqx_S_x_K = term1 + term2 + term3;

			double cos_qx_S_x_K = eiqx_S_x_K.real();
			double sin_qx_S_x_K = eiqx_S_x_K.imag();

			cosqx_dN_dypTdpTdphi[ipY] += cos_ta * cos_qx_S_x_K + sin_ta * sin_qx_S_x_K;
			sinqx_dN_dypTdpTdphi[ipY] += cos_ta * sin_qx_S_x_K - sin_ta * cos_qx_S_x_K;
		}
	}

	return;
}

//End of file
