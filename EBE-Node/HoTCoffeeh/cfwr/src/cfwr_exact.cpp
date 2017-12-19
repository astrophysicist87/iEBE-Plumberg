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
#include "cfwr_lib.h"
#include "gauss_quadrature.h"
#include "bessel.h"

using namespace std;

const std::complex<double> i(0, 1);

inline void Iexact(double alpha, double beta, double gamma, complex<double> & I0, complex<double> & I1, complex<double> & I2, complex<double> & I3)
{
	complex<double> ci0, ci1, ck0, ck1, ci0p, ci1p, ck0p, ck1p;
	complex<double> z0 = alpha - i*beta;
	complex<double> z0sq = pow(z0, 2.0);
	double gsq = gamma*gamma;
	complex<double> z = sqrt(z0sq + gsq);
	int errorCode = bessf::cbessik01(z, ci0, ci1, ck0, ck1, ci0p, ci1p, ck0p, ck1p);
	
	I0 = 2.0*ck0;
	I1 = 2.0*z0*ck1 / z;
	I2 = 2.0*z0sq*ck0 / (z*z)
			+ 2.0*(z0sq - gsq)*ck1 / pow(z, 3.0);
	I3 = 2.0*z0*( ( pow(z0, 4.0) - 2.0* z0sq*gsq - 3.0 * pow(gamma, 4.0) ) * ck0 / z
						+ (-6.0*gsq + z0sq*(2.0 + z0sq + gsq)) * ck1
				) / pow(z,5.0);

	return;
}

void CorrelationFunction::Cal_dN_dypTdpTdphi_with_weights_function_approx(int local_pid, double pT, double pphi, double p_Y,
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

	//double pY_shift = - double(abs(qz)>1.e-10) * asinh(qz / sqrt(abs(qt*qt-qz*qz) + 1.e-100));
	double pY_shift = 0.5 * log(abs((qt+qz + 1.e-100)/(qt-qz + 1.e-100)));
	//double pY_shift = 0.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double mT = sqrt(pT*pT+localmass*localmass);
	*cosqx_dN_dypTdpTdphi = 0.0;
	*sinqx_dN_dypTdpTdphi = 0.0;

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

		double transverse_arg = (xpt*qx+ypt*qy)/hbarC;
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

		double ch_pY = cosh(p_Y+pY_shift);
		double sh_pY = sinh(p_Y+pY_shift);
		double beta = (tau / hbarC) * ( qt*ch_pY - qz*sh_pY );
		double gamma = (tau / hbarC) * ( qz*ch_pY - qt*sh_pY );

		complex<double> I0_a_b_g, I1_a_b_g, I2_a_b_g, I3_a_b_g;
		complex<double> I0_2a_b_g, I1_2a_b_g, I2_2a_b_g, I3_2a_b_g;
		Iexact(alpha, beta, gamma, I0_a_b_g, I1_a_b_g, I2_a_b_g, I3_a_b_g);
		Iexact(2.0*alpha, beta, gamma, I0_2a_b_g, I1_2a_b_g, I2_2a_b_g, I3_2a_b_g);

		complex<double> term1 = transverse_f0 * (A*I1_a_b_g + B*I0_a_b_g);
		complex<double> term2 = C * transverse_f0 * ( A*a*I3_a_b_g + (B*a+b*A)*I2_a_b_g + (B*b+c*A)*I1_a_b_g + B*c*I0_a_b_g );
		complex<double> term3 = -sign * C * transverse_f0 * transverse_f0 * ( A*a*I3_2a_b_g + (B*a+b*A)*I2_2a_b_g + (B*b+c*A)*I1_2a_b_g + B*c*I0_2a_b_g );

		complex<double> eiqx_S_x_K = term1 + term2 + term3;

		double cos_qx_S_x_K = eiqx_S_x_K.real();
		double sin_qx_S_x_K = eiqx_S_x_K.imag();

		*cosqx_dN_dypTdpTdphi += cos_ta * cos_qx_S_x_K + sin_ta * sin_qx_S_x_K;
		*sinqx_dN_dypTdpTdphi += cos_ta * sin_qx_S_x_K - sin_ta * cos_qx_S_x_K;

//if (abs(p_Y)<1.e-10) cout << "CHECK(direct): " << isurf << "   " << qt << "   " << qx << "   " << qy << "   " << qz << "   " << pT << "   " << pphi << "   "
//		<< cos_qx_S_x_K << "   " << sin_qx_S_x_K << "   " << cos_ta << "   " << sin_ta << endl;

	}

	return;
}



void CorrelationFunction::Cal_dN_dypTdpTdphi_with_weights_function_approx(int local_pid, double pT, double pphi, double Del_p_Y,
					double qt, double qx, double qy, double qz,
					double * cosLcosT_dN_dypTdpTdphi, double * cosLsinT_dN_dypTdpTdphi,
					double * sinLcosT_dN_dypTdpTdphi, double * sinLsinT_dN_dypTdpTdphi)
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

	double pY_shift = 0.5 * log(abs((qt+qz + 1.e-100)/(qt-qz + 1.e-100)));
	//double pY_shift = 0.0;

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double mT = sqrt(pT*pT+localmass*localmass);
	*cosLcosT_dN_dypTdpTdphi = 0.0;
	*cosLsinT_dN_dypTdpTdphi = 0.0;
	*sinLcosT_dN_dypTdpTdphi = 0.0;
	*sinLsinT_dN_dypTdpTdphi = 0.0;

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

		double transverse_arg = (xpt*qx+ypt*qy)/hbarC;
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

		double ch_pY = cosh(Del_p_Y+pY_shift);
		double sh_pY = sinh(Del_p_Y+pY_shift);
		double beta = (tau / hbarC) * ( qt*ch_pY - qz*sh_pY );
		double gamma = (tau / hbarC) * ( qz*ch_pY - qt*sh_pY );

		complex<double> I0_a_b_g, I1_a_b_g, I2_a_b_g, I3_a_b_g;
		complex<double> I0_2a_b_g, I1_2a_b_g, I2_2a_b_g, I3_2a_b_g;
		Iexact(alpha, beta, gamma, I0_a_b_g, I1_a_b_g, I2_a_b_g, I3_a_b_g);
		Iexact(2.0*alpha, beta, gamma, I0_2a_b_g, I1_2a_b_g, I2_2a_b_g, I3_2a_b_g);

		complex<double> term1 = transverse_f0 * (A*I1_a_b_g + B*I0_a_b_g);
		complex<double> term2 = C * transverse_f0 * ( A*a*I3_a_b_g + (B*a+b*A)*I2_a_b_g + (B*b+c*A)*I1_a_b_g + B*c*I0_a_b_g );
		complex<double> term3 = -sign * C * transverse_f0 * transverse_f0 * ( A*a*I3_2a_b_g + (B*a+b*A)*I2_2a_b_g + (B*b+c*A)*I1_2a_b_g + B*c*I0_2a_b_g );

		complex<double> eiqx_S_x_K = term1 + term2 + term3;

		double cos_qx_S_x_K = eiqx_S_x_K.real();
		double sin_qx_S_x_K = eiqx_S_x_K.imag();

		*cosLcosT_dN_dypTdpTdphi += cos_ta * cos_qx_S_x_K;
		*cosLsinT_dN_dypTdpTdphi -= sin_ta * cos_qx_S_x_K;	//for consistency with thermal calculation
		*sinLcosT_dN_dypTdpTdphi += cos_ta * sin_qx_S_x_K;
		*sinLsinT_dN_dypTdpTdphi += sin_ta * sin_qx_S_x_K;
	}

	return;
}

//define some parameters for the exact emission function
const double Rad = 5.0, Del_tau = 1.0, tau0 = 5.0, etaf = 0.6;

inline double Hfactor(double r, double tau)
{
	return (
			exp( -r*r/(2.0*Rad*Rad) - (tau-tau0)*(tau-tau0)/(2.0*Del_tau*Del_tau) ) / (M_PI*Del_tau)
			);
}

inline double eta_t(double r)
{
	return ( etaf*r/Rad );
}

void CorrelationFunction::Cal_dN_dypTdpTdphi_no_weights_toy(int local_pid)
{
	//space-time integration grid
	const int n_tau_pts = 51;
	const int n_r_pts = 81;
	const int n_phi_pts = 31;
	double * tau_pts = new double [n_tau_pts];
	double * tau_wts = new double [n_tau_pts];
	double * x_pts = new double [n_r_pts];
	double * x_wts = new double [n_r_pts];
	double * phi_pts = new double [n_phi_pts];
	double * phi_wts = new double [n_phi_pts];
	gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, 0.0, 10.0, tau_pts, tau_wts);
	gauss_quadrature(n_r_pts, 1, 0.0, 0.0, -1.0, 1.0, x_pts, x_wts);
	gauss_quadrature(n_phi_pts, 1, 0.0, 0.0, 0.0, 2.0*M_PI, phi_pts, phi_wts);

	// set particle information
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;

	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	{
//if (ipT > 0 || ipphi > 0) continue;
		double pT = SP_pT[ipT];
		double pphi = SP_pphi[ipphi];
		double mT = sqrt(pT*pT+localmass*localmass);

		double spectra_at_pTpphi = 0.0;

		for (int ir = 0; ir < n_r_pts; ++ir)
		{

			//set r-point inside pT loop, since optimal distribution of integration points
			// depends on value of MT
			double rmin = 0.0, rmax = 5.0 * Rad / sqrt(1.0 + mT*one_by_Tdec*etaf*etaf);
			double hw = 0.5 * (rmax - rmin), cen = 0.5 * (rmax + rmin);
			double rpt = cen + hw * x_pts[ir];

			double ch_eta_t = cosh(eta_t(rpt));
			double sh_eta_t = sinh(eta_t(rpt));

			double alpha = mT*ch_eta_t*one_by_Tdec;
			complex<double> I0, I1, I2, I3;
			Iexact(alpha, 0.0, 0.0, I0, I1, I2, I3);

			for (int itau = 0; itau < n_tau_pts; ++itau)
			{
				double tau = tau_pts[itau];
				double local_H = Hfactor(rpt, tau);

				for (int iphi = 0; iphi < n_phi_pts; ++iphi)
				{
					double phipt = phi_pts[iphi];
//if (local_pid != target_particle_id)
//	cout << "SPECTRA: " << local_pid << "   " << ipT << "   " << ipphi << "   " << ir << "   " << itau << "   " << iphi << "   ";
// << "   " << ch_eta_t << "   " << sh_eta_t << "   ";
//	cout << tau_wts[itau] << "   " << r_wts[ir] << "   " << phi_wts[iphi] << "   " << mT << "   " << tau << "   " << rpt << "   ";
//	cout << prefactor << "   " << local_H << "   " << one_by_Tdec << "   " << pT << "   " << phipt << "   " << pphi << "   " << "   " << cos(phipt - pphi)
//			<< "   " << exp( one_by_Tdec*pT*sh_eta_t*cos(phipt - pphi) ) << "   ";
					double S_p_with_weight = tau_wts[itau]*hw*x_wts[ir]*phi_wts[iphi]*mT*tau*rpt*prefactor
												*local_H*exp( one_by_Tdec*pT*sh_eta_t*cos(phipt - pphi) );
//if (local_pid != target_particle_id)
//	cout << S_p_with_weight << "   " << I1.real() << "   " << I1.imag() << endl;

					spectra_at_pTpphi += S_p_with_weight*I1.real();
				}
			}
		}

		//update spectra
		spectra[local_pid][ipT][ipphi] = spectra_at_pTpphi;
		thermal_spectra[local_pid][ipT][ipphi] = spectra_at_pTpphi;
		log_spectra[local_pid][ipT][ipphi] = log(abs(spectra_at_pTpphi) + 1.e-100);
		sign_spectra[local_pid][ipT][ipphi] = sgn(spectra_at_pTpphi);
//		cout << "FINAL SPECTRA: " << all_particles[local_pid].name << "   " << SP_pT[ipT] << "   " << SP_pphi[ipphi] << "   " << spectra[local_pid][ipT][ipphi] << endl;
	}

	//clean up
	delete [] tau_pts;
	delete [] tau_wts;
	delete [] x_pts;
	delete [] x_wts;
	delete [] phi_pts;
	delete [] phi_wts;

	return;
}

void CorrelationFunction::Cal_dN_dypTdpTdphi_with_weights_toy(int local_pid, int iqt, int iqz, int ipY, double * moments_to_update)
{
	//space-time integration grid
	const int n_tau_pts = 51;
	const int n_r_pts = 81;
	const int n_phi_pts = 31;
	double * tau_pts = new double [n_tau_pts];
	double * tau_wts = new double [n_tau_pts];
	double * x_pts = new double [n_r_pts];
	double * x_wts = new double [n_r_pts];
	double * phi_pts = new double [n_phi_pts];
	double * phi_wts = new double [n_phi_pts];
	gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, 0.0, 10.0, tau_pts, tau_wts);
	gauss_quadrature(n_r_pts, 1, 0.0, 0.0, -1.0, 1.0, x_pts, x_wts);
	gauss_quadrature(n_phi_pts, 1, 0.0, 0.0, 0.0, 2.0*M_PI, phi_pts, phi_wts);
	/*for (int itau = 0; itau < n_tau_pts; itau++)
		cout << itau << "   " << tau_pts[itau] << endl;
	for (int ir = 0; ir < n_r_pts; ir++)
		cout << ir << "   " << r_pts[ir] << endl;
	for (int iphi = 0; iphi < n_phi_pts; iphi++)
		cout << iphi << "   " << phi_pts[iphi] << endl;
	if (1) exit (1);*/

	//q-point info
	double qt_loc = qt_pts[iqt], qz_loc = qz_pts[iqz];

	//zero arrays
	if (local_pid == target_particle_id)
	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int ii = 0; ii < 4; ++ii)
		moments_to_update[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, ii)] = 0.0;
	else
	for (int ipT = 0; ipT < n_pT_pts; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int ii = 0; ii < 4; ++ii)
		moments_to_update[fixQTQZ_indexer(ipT, ipphi, ipY, iqx, iqy, ii)] = 0.0;


	// set particle information
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;

	double ch_pY = ch_SP_pY[ipY];
	double sh_pY = sh_SP_pY[ipY];

	if (local_pid == target_particle_id)
	{
		for (int ir = 0; ir < n_r_pts; ++ir)
		for (int itau = 0; itau < n_tau_pts; ++itau)
		{
			double tau = tau_pts[itau];
			double beta = tau * hbarCm1 * ( qt_loc*ch_pY - qz_loc*sh_pY );
			double gamma = tau * hbarCm1 * ( qz_loc*ch_pY - qt_loc*sh_pY );

//cout << "CHECKbg: " << beta << "   " << gamma << "   " << tau << "   " << hbarCm1 << "   " << qt_loc << "   " << ch_pY << "   " << qz_loc << "   " << sh_pY << endl;

			for (int ipT = 0; ipT < n_pT_pts; ++ipT)
			{
//if (ipT > 0) continue;
				double pT = SP_pT[ipT];
				double mT = sqrt(pT*pT+localmass*localmass);

				//set r-point inside pT loop, since optimal distribution of integration points
				// depends on value of MT
				double rmin = 0.0, rmax = 7.5 * Rad / sqrt(1.0 + mT*one_by_Tdec*etaf*etaf);
				double hw = 0.5 * (rmax - rmin), cen = 0.5 * (rmax + rmin);
				double rpt = cen + hw * x_pts[ir];

				double local_H = Hfactor(rpt, tau);
				double ch_eta_t = cosh(eta_t(rpt));
				double sh_eta_t = sinh(eta_t(rpt));

				double alpha = mT * one_by_Tdec * ch_eta_t;
				complex<double> I0, I1, I2, I3;
				Iexact(alpha, beta, gamma, I0, I1, I2, I3);
				double cos_phi_L = I1.real();
				double sin_phi_L = I1.imag();

				for (int iphi = 0; iphi < n_phi_pts; ++iphi)
				{
					double phipt = phi_pts[iphi];
					double xpt = rpt*cos(phipt), ypt = rpt*sin(phipt);

					for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
					{
//if (ipphi > 0) continue;
						double pphi = SP_pphi[ipphi];

//if (local_pid != target_particle_id)
//	cout << "MOMENTS: " << local_pid << "   " << ir << "   " << itau << "   " << iphi << "   ";
						double S_p_with_weight = tau_wts[itau]*hw*x_wts[ir]*phi_wts[iphi]*mT*tau*rpt*prefactor*local_H*exp( one_by_Tdec*pT*sh_eta_t*cos(phipt - pphi) );

//if (local_pid != target_particle_id)
//	cout << S_p_with_weight << "   " << I1.real() << "   " << I1.imag() << endl;

						for (int iqx = 0; iqx < qxnpts; ++iqx)
						for (int iqy = 0; iqy < qynpts; ++iqy)
						{
							double cosAx = cos(qx_pts[iqx]*xpt/hbarC), sinAx = sin(qx_pts[iqx]*xpt/hbarC);
							double cosAy = cos(qy_pts[iqy]*ypt/hbarC), sinAy = sin(qy_pts[iqy]*ypt/hbarC);
							double cos_trans_Fourier = cosAx*cosAy - sinAx*sinAy;	//==cos(qx x + qy y)
							double sin_trans_Fourier = sinAx*cosAy + cosAx*sinAy;	//==sin(qx x + qy y)
							moments_to_update[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, 0)] += S_p_with_weight * cos_phi_L * cos_trans_Fourier;
							moments_to_update[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, 1)] -= S_p_with_weight * cos_phi_L * sin_trans_Fourier;
							moments_to_update[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, 2)] += S_p_with_weight * sin_phi_L * cos_trans_Fourier;
							moments_to_update[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, 3)] += S_p_with_weight * sin_phi_L * sin_trans_Fourier;
						}
					}
				}
			}
		}
	}
	else
	{
		for (int itau = 0; itau < n_tau_pts; ++itau)
		for (int ir = 0; ir < n_r_pts; ++ir)
		{
			double tau = tau_pts[itau];
			double beta = tau * hbarCm1 * ( qt_loc*ch_pY - qz_loc*sh_pY );
			double gamma = tau * hbarCm1 * ( qz_loc*ch_pY - qt_loc*sh_pY );
//cout << "CHECKbg: " << beta << "   " << gamma << "   " << tau << "   " << hbarCm1 << "   " << qt_loc << "   " << ch_pY << "   " << qz_loc << "   " << sh_pY << endl;

			for (int ipT = 0; ipT < n_pT_pts; ++ipT)
			{
				double pT = SP_pT[ipT];
				double mT = sqrt(pT*pT+localmass*localmass);

				//set r-point inside pT loop, since optimal distribution of integration points
				// depends on value of MT
				double rmin = 0.0, rmax = 7.5 * Rad / sqrt(1.0 + mT*one_by_Tdec*etaf*etaf);
				double hw = 0.5 * (rmax - rmin), cen = 0.5 * (rmax + rmin);
				double rpt = cen + hw * x_pts[ir];

				double local_H = Hfactor(rpt, tau);
				double ch_eta_t = cosh(eta_t(rpt));
				double sh_eta_t = sinh(eta_t(rpt));

				double alpha = mT * one_by_Tdec * ch_eta_t;
				complex<double> I0, I1, I2, I3;
				Iexact(alpha, beta, gamma, I0, I1, I2, I3);
				double cos_phi_L = I1.real();
				double sin_phi_L = I1.imag();

				for (int iphi = 0; iphi < n_phi_pts; ++iphi)
				{
					double phipt = phi_pts[iphi];
					double xpt = rpt*cos(phipt), ypt = rpt*sin(phipt);

					for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
					{
						double pphi = SP_pphi[ipphi];

	//if (ipY==ipY0)
	//	cout << "MOMENTS: " << local_pid << "   " << iqt << "   " << ipT << "   " << ipphi << "   " << ir << "   " << itau << "   " << iphi << "   ";
						double S_p_with_weight = tau_wts[itau]*hw*x_wts[ir]*phi_wts[iphi]*mT*tau*rpt*prefactor*local_H*exp( one_by_Tdec*pT*sh_eta_t*cos(phipt - pphi) );
	//if (ipY==ipY0)
	//	cout << S_p_with_weight << "   " << I1.real() << "   " << I1.imag() << endl;

						for (int iqx = 0; iqx < qxnpts; ++iqx)
						for (int iqy = 0; iqy < qynpts; ++iqy)
						{
							double cosAx = cos(qx_pts[iqx]*xpt/hbarC), sinAx = sin(qx_pts[iqx]*xpt/hbarC);
							double cosAy = cos(qy_pts[iqy]*ypt/hbarC), sinAy = sin(qy_pts[iqy]*ypt/hbarC);
							double cos_trans_Fourier = cosAx*cosAy - sinAx*sinAy;	//==cos(qx x + qy y)
							double sin_trans_Fourier = sinAx*cosAy + cosAx*sinAy;	//==sin(qx x + qy y)
							moments_to_update[fixQTQZ_indexer(ipT, ipphi, ipY, iqx, iqy, 0)] += S_p_with_weight * cos_phi_L * cos_trans_Fourier;
							moments_to_update[fixQTQZ_indexer(ipT, ipphi, ipY, iqx, iqy, 1)] -= S_p_with_weight * cos_phi_L * sin_trans_Fourier;
							moments_to_update[fixQTQZ_indexer(ipT, ipphi, ipY, iqx, iqy, 2)] += S_p_with_weight * sin_phi_L * cos_trans_Fourier;
							moments_to_update[fixQTQZ_indexer(ipT, ipphi, ipY, iqx, iqy, 3)] += S_p_with_weight * sin_phi_L * sin_trans_Fourier;
						}
					}
				}
			}
		}
	}

	/*if (ipY==ipY0)
	{
		if (local_pid==target_particle_id)
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
			cout << "FINAL MOMENTS: " << all_particles[local_pid].name << "   " << iqt << "   " << ipY << "   " << SP_pT[ipT] << "   " << SP_pphi[ipphi]
					<< "   " << moments_to_update[indexer(ipT, ipphi, iqt, iqx, iqy, iqz, 0)] << endl;
		else
		for (int ipT = 0; ipT < n_pT_pts; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_pts; ++ipphi)
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
			cout << "FINAL MOMENTS: " << all_particles[local_pid].name << "   " << iqt << "   " << ipY << "   " << SP_pT[ipT] << "   " << SP_pphi[ipphi]
					<< "   " << moments_to_update[fixQTQZ_indexer(ipT, ipphi, ipY, iqx, iqy, 0)] << endl;
	}*/

	//clean up
	delete [] tau_pts;
	delete [] tau_wts;
	delete [] x_pts;
	delete [] x_wts;
	delete [] phi_pts;
	delete [] phi_wts;

	return;
}

//End of file
