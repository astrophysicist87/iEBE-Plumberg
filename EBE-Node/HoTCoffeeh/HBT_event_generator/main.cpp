#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <sys/time.h>

#include "src/Stopwatch.h"
//#include "src/HBT_event_generator.h"
#include "src/ParameterReader.h"
#include "src/EventRecord.h"
#include "src/ParticleRecord.h"
#include "main.h"

using namespace std;

int main(int argc, char *argv[])
{
	// Display intro
	cout << endl
			<< "              HBT event generator              " << endl
			<< endl
			<< "  Ver 1.0   ----- Christopher Plumberg, 10/2018" << endl;
	//cout << endl << "**********************************************************" << endl;
	//display_logo(2); // Hail to the king~
	//cout << endl << "**********************************************************" << endl << endl;
   

	// Read-in parameters
	/*ParameterReader * paraRdr = new ParameterReader;
	paraRdr->readFromFile("parameters.dat");
	paraRdr->readFromArguments(argc, argv);
	paraRdr->echo();*/


	// Start timing
	Stopwatch sw;
	sw.Start();


	// Specify path
    string workingDirectory = "./results";


	// Load pdg.dat and any relevant resonance information
	//int particle_idx = 1;	//pion^+


	// Specify files containing all position-momentum information
	// from which to construct HBT correlation function
	vector<string> all_file_names;
	read_file_catalogue("./catalogue.dat", all_file_names);


	// Read in the files
	vector<EventRecord> allEvents, eventsInFile;

	allEvents.clear();
	for (int iFile = 0; iFile < all_file_names.size(); ++iFile)
	{
		// Reset
		eventsInFile.clear();


		// Read in this file
		read_in_file(all_file_names[iFile], eventsInFile);


		// Concatenate these events to allEvents vector
		allEvents.insert( allEvents.end(),
							eventsInFile.begin(),
							eventsInFile.end() );
	}

	// check that everything was read in correctly
	/*for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	{
		EventRecord event = allEvents[iEvent];
		
		for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
		{
			ParticleRecord p = event.particles[iParticle];
			cout << "Check: "
					<< p.eventID << "   " << p.particleID << "   "
					<< p.E << "   " << p.px << "   "
					<< p.py << "   " << p.pz << "   "
					<< p.t << "   " << p.x << "   "
					<< p.y << "   " << p.z << endl;
		}
	}*/
	

	// Create HBT_event_generator object from allEvents
	HBT_event_generator HBT_event_ensemble(paraRdr, allEvents);


	// Compute correlation function
	HBT_event_ensemble.Compute_correlation_function();


	// Output correlation function
	HBT_event_ensemble.Output_correlation_function();


	// Print out run-time
	sw.Stop();
	cout 	<< "Finished everything in "
			<< sw.printTime() << " seconds." << endl;

	// Wrap it up!
	return (0);
}

//End of file
