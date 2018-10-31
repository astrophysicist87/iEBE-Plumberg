#ifndef MAIN_H
#define MAIN_H

#include <string>
#include <sstream>
#include <fstream>

#include "src/EventRecord.h"
#include "src/ParticleRecord.h"

using namespace std;

//this is just to give this file a reason to exist for the moment...
const double plumbergtest = 0.;


// function to read in catalogue of event files and return number of files to read
int read_file_catalogue(string catalogue_name, vector<string> & allFilenames)
{
	ifstream catalogue_in(catalogue_name.c_str());

	string filename;
	while (getline(catalogue_in, filename))
		allFilenames.push_back(filename);

	return ( allFilenames.size() );
}



// function to read in a file containing some number of events
void read_in_file(string filename, vector<EventRecord> & eventsInFile)
{
	ifstream infile(filename.c_str());

	int count = 0;
	string line;
	int previous_eventID = -1, current_eventID = -1;
	
	EventRecord event;

	while (getline(infile, line))
	{
		istringstream iss(line);

		ParticleRecord particle;
		int eventID, particleID;
		double E, px, py, pz;
		double t, x, y, z;

		if ( !( iss >> eventID
					>> particleID
					>> E >> px >> py >> pz
					>> t >> x >> y >> z
			 ) ) { break; }

		particle.eventID 	= eventID;
		particle.particleID = particleID;
		particle.E 			= E;
		particle.px 		= px;
		particle.py 		= py;
		particle.pz 		= pz;
		particle.t 			= t;
		particle.x 			= x;
		particle.y 			= y;
		particle.z 			= z;

		// Decide what to do with new particle
		// if on first iteration
		if (count == 0)
		{
			// initialize previous eventID
			// and current eventID
			previous_eventID = eventID;
			current_eventID = eventID;

			// push particle to event
			event.particles.push_back(particle);
		}
		// otherwise...
		else
		{
			current_eventID = eventID;

			// if newest particle does not
			// correspond to a new event
			if (current_eventID == previous_eventID)
			{
				event.particles.push_back(particle);
			}
			// if newest particle corresponds
			// to a new event
			else
			{
				// push event to eventsInFile
				eventsInFile.push_back(event);

				// reset event
				event = EventRecord();

				// push new particle to new event
				event.particles.push_back(particle);
			}
		}

		previous_eventID = current_eventID;
		++count;
	}

	// push final event to eventsInFile
	eventsInFile.push_back(event);

	infile.close();

	return;
}

#endif
