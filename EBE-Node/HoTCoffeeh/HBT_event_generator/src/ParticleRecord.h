#ifndef PARTICLERECORD_H
#define PARTICLERECORD_H

#include <vector>

typedef struct
{
	int eventID;		//which event did this particle come from
	//int pdgID;			//particle type for accessing info (e.g., mass)
	int particleID; 	//to distinguish it from other particles in the same event

	//vector<double> x;	//where was this particle produced?
	double t, x, y, z, r, eta, phi, tau;
	//vector<double> p;	//what momentum was this particle produced with?
	double E, px, py, pz, pT, pY, pphi;
	
} ParticleRecord;

#endif
