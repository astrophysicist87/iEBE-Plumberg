#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <complex>

#include "HBT_event_generator.h"
#include "Arsenal.h"
#include "Stopwatch.h"

using namespace std;

void HBT_event_generator::initialize_all(
	ParameterReader * paraRdr_in,
	const vector<EventRecord> & allEvents_in )
{
	// Load parameters
	paraRdr = paraRdr_in;

	// Copy in records of all events
	allEvents = allEvents_in;

	//Set header info
	//Define various grid sizes
	// - SP momentum points at which to evaluate correlation function
	n_pT_pts 		= paraRdr->getVal("n_pT_pts");
	pT_min 			= paraRdr->getVal("pTmin");
	pT_max 			= paraRdr->getVal("pTmax");
	n_pphi_pts 		= paraRdr->getVal("n_pphi_pts");
	pphi_min 		= -M_PI+1.e-10;
	pphi_max 		= M_PI-1.e-10;
	n_pY_pts 		= paraRdr->getVal("n_pY_pts");
	pY_min 			= paraRdr->getVal("pYmin");
	pY_max 			= paraRdr->getVal("pYmax");
	// - pair momenta points at which to interpolate HBT results
	n_KT_pts 		= paraRdr->getVal("n_KT_pts");
	KT_min 			= paraRdr->getVal("KTmin");
	KT_max 			= paraRdr->getVal("KTmax");
	n_Kphi_pts 		= paraRdr->getVal("n_Kphi_pts");
	Kphi_min 		= -M_PI+1.e-10;
	Kphi_max 		= M_PI-1.e-10;
	n_KL_pts 		= paraRdr->getVal("n_KL_pts");
	KL_min 			= paraRdr->getVal("KLmin");
	KL_max 			= paraRdr->getVal("KLmax");
	// - relative momentum points at which to evaluate
	//   correlation function
	n_qo_pts 		= paraRdr->getVal("n_qo_pts");
	n_qs_pts 		= paraRdr->getVal("n_qs_pts");
	n_ql_pts 		= paraRdr->getVal("n_ql_pts");
	// - step size in q directions
	delta_qo 		= paraRdr->getVal("delta_qo");
	delta_qs 		= paraRdr->getVal("delta_qs");
	delta_ql 		= paraRdr->getVal("delta_ql");
	// - minimum value in each q direction
	init_qo 		= -0.5*double(n_qo_pts-1)*delta_qo;
	init_qs 		= -0.5*double(n_qs_pts-1)*delta_qs;
	init_ql 		= -0.5*double(n_ql_pts-1)*delta_ql;

	// - number of points to use when fleshing out correlation
	//   function in each direction
	//new_nqopts 		= ( n_qo_pts > 1 ) ? new_nqpts : 1;
	//new_nqspts 		= ( n_qs_pts > 1 ) ? new_nqpts : 1;
	//new_nqlpts 		= ( n_ql_pts > 1 ) ? new_nqpts : 1;

	n_qo_bins 		= n_qo_pts - 1;
	n_qs_bins 		= n_qs_pts - 1;
	n_ql_bins 		= n_ql_pts - 1;

	n_pT_bins 		= n_pT_pts  - 1;
	n_pphi_bins 	= n_pphi_pts  - 1;
	n_pY_bins 		= n_pY_pts  - 1;

	n_KT_bins 		= n_KT_pts - 1;
	n_Kphi_bins 	= n_Kphi_pts - 1;
	n_KL_bins 		= n_KL_pts - 1;

	pT_pts 			= vector<double> (n_pT_pts);
	pphi_pts 		= vector<double> (n_pphi_pts);
	pY_pts 			= vector<double> (n_pY_pts);

	KT_pts 			= vector<double> (n_KT_pts);
	Kphi_pts 		= vector<double> (n_Kphi_pts);
	KL_pts 			= vector<double> (n_KL_pts);

	dN_pTdpTdpphidpY = vector<double> (n_pT_bins*n_pphi_bins*n_pY_pts);

	linspace(pT_pts, pT_min, pT_max);
	linspace(pphi_pts, pphi_min, pphi_max);
	linspace(pY_pts, pY_min, pY_max);

	linspace(KT_pts, KT_min, KT_max);
	linspace(Kphi_pts, Kphi_min, Kphi_max);
	linspace(KL_pts, KL_min, KL_max);

	pT_bin_width 	= pT_pts[1]-pT_pts[0];
	pphi_bin_width 	= pphi_pts[1]-pphi_pts[0];
	pY_bin_width 	= pY_pts[1]-pY_pts[0];

	KT_bin_width 	= KT_pts[1]-KT_pts[0];
	Kphi_bin_width 	= Kphi_pts[1]-Kphi_pts[0];
	KL_bin_width 	= KL_pts[1]-KL_pts[0];

	numerator 		= vector<double> (n_KT_bins*n_Kphi_bins*n_KL_bins*n_qo_bins*n_qs_bins*n_ql_bins);

	return;
}

HBT_event_generator::~HBT_event_generator()
{
	//clear everything

	return;
}


/*
void HBT_event_generator::Compute_numerator()
{
	double sumOverEvents = 0.0;

	// Sum over all events
	for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	{
		EventRecord event = allEvents[iEvent];

		int KPhaseSpaceSize = n_KT_bins*n_Kphi_bins*n_KL_bins;
		int qPhaseSpaceSize = n_qo_bins*n_qs_bins*n_ql_bins;
		vector<complex<double> > sumOverParticles1 (KPhaseSpaceSize*qPhaseSpaceSize);
		vector<complex<double> > sumOverParticles2 (KPhaseSpaceSize);
		
		for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
		{
			ParticleRecord p = event.particles[iParticle];

			int iKT 	= bin_function(p.pT, KT_pts);
			int iKphi 	= bin_function(p.pphi, Kphi_pts);
			int iKL 	= bin_function(p.pz, KL_pts);

			// Check if particle is out of range
			if ( 	   iKT < 0 or iKT >= n_KT_bins
					or iKL < 0 or iKL >= n_KL_bins
					//or iqo   < 0 or iqo   >= n_qo_bins
					//or iqs   < 0 or iqs   >= n_qs_bins
					//or iql   < 0 or iql   >= n_ql_bins
				) continue;

			//periodicity
			if ( iKphi < 0 )
				err << "Warning: iKphi = " << iKphi << " < 0!" << endl;
			else if ( iKphi > n_Kphi_bins )
				err << "Warning: iKphi = " << iKphi << " > n_Kphi_bins!" << endl;

			// Get particle mass (can also read this in)
			double m = sqrt(p.E*p.E - p.px*p.px - p.py*p.py - p.pz*p.pz);

			// Get bin centers
			double KT = 0.5*(KT_pts[iKT]+KT_pts[iKT+1]);
			double Kphi = 0.5*(Kphi_pts[iKphi]+Kphi_pts[iKphi+1]);
			double cKphi = cos(Kphi), sKphi = sin(Kphi);
			double KL = 0.5*(KL_pts[iKL]+KL_pts[iKL+1]);

			// Rotate from xyz to osl
			double xo = p.x * cKphi + p.y * sKphi;
			double xs = -p.x * sKphi + p.y * cKphi;
			double xl = p.z;

			// If not, include it in the sums
			// Loop over q bins
			for (int iqo = 0; iqo < n_qo_bins; iqo++)
			for (int iqs = 0; iqs < n_qs_bins; iqs++)
			for (int iql = 0; iql < n_ql_bins; iql++)
			{
				double qo = 0.5*(qo_pts[iqo]+qo_pts[iqo+1]);
				double qs = 0.5*(qs_pts[iqs]+qs_pts[iqs+1]);
				double ql = 0.5*(ql_pts[iql]+ql_pts[iql+1]);
				double q0 = get_q0(m, qo, qs, ql, KT, KL);

				double arg = q0*p.t - (qo*xo+qs*xs+ql*xl);
				complex<double> phase_factor = exp(i*arg);

				int loc_index = indexer(iKT, iKphi, iKL, iqo, iqs, iql);
				sumOverParticles1[loc_index] += phase_factor;
			}

			sumOverParticles2[indexer(iKT, iKphi, iKL)] += 1.0;



			//cout << "Check: "
			//		<< p.eventID << "   " << p.particleID << "   "
			//		<< p.E << "   " << p.px << "   "
			//		<< p.py << "   " << p.pz << "   "
			//		<< p.t << "   " << p.x << "   "
			//		<< p.y << "   " << p.z << endl;
		}
		
		// Normalize sums appropriately
		for (int iKT = 0; iKT < n_KT_bins; iKT++)
		for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
		for (int iKL = 0; iKL < n_KL_bins; iKL++)
		{
			double KT_bin_center = 0.5*(KT_pts[iKT]+KT_pts[iKT+1]);
			sumOverParticles2[indexer(iKT, iKphi, iKL)]
				/= KT_bin_center*KT_bin_width
					*Kphi_bin_width*KL_bin_width
					*KT_bin_center*KT_bin_width
					*Kphi_bin_width*KL_bin_width;

			complex<double> sum2 = sumOverParticles2[indexer(iKT, iKphi, iKL)];

			for (int iqo = 0; iqo < n_qo_bins; iqo++)
			for (int iqs = 0; iqs < n_qs_bins; iqs++)
			for (int iql = 0; iql < n_ql_bins; iql++)
			{
				sumOverParticles1[indexer(iKT, iKphi, iKL, iqo, iqs, iql)]
					/= KT_bin_center*KT_bin_width
						*Kphi_bin_width*KL_bin_width
						*delta_qo*delta_qs*delta_ql;
				complex<double> sum1 = sumOverParticles1[indexer(iKT, iKphi, iKL)];

				numerator[indexer(iKT, iKphi, iKL, iqo, iqs, iql)]
					+= sum*conj(sum1) - sum2;
			}

		}


	}

	return;
}*/


/*
void HBT_event_generator::Compute_denominator()
{
	double sumOverEvents = 0.0;

	// Sum over all events
	for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	{
		EventRecord event = allEvents[iEvent];

		double sumOverParticles1 = 0.0;
		double sumOverParticles2 = 0.0;
		
		for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
		{
			ParticleRecord p = event.particles[iParticle];

		}
	}

	return;
}
*/










//End of file
