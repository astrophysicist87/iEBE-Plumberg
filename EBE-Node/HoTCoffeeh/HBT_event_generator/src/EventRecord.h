#ifndef EVENTRECORD_H
#define EVENTRECORD_H

#include <vector>

#include "ParticleRecord.h"

typedef struct
{
	int eventID;	//to distinguish it from other events in the same ensemble
	int Nparticles;				//may not need this

	vector<int> particleIDs;	//may not need this
	vector<ParticleRecord> particles;
	
} EventRecord;

#endif
