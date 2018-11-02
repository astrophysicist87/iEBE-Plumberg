#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <sys/time.h>

#include "src/Stopwatch.h"
#include "src/HBT_event_generator.h"
#include "src/ParameterReader.h"
#include "src/EventRecord.h"
#include "src/ParticleRecord.h"
#include "main.h"

using namespace std;

int main(int argc, char *argv[])
{
	// Display intro
	//cout << endl
	//		<< "              HBT event generator              " << endl
	//		<< endl
	//		<< "  Ver 1.0   ----- Christopher Plumberg, 10/2018" << endl;
	//cout << endl << "**********************************************************" << endl;
	//display_logo(2); // Hail to the king~
	//cout << endl << "**********************************************************" << endl << endl;
   

	// Read-in parameters
	ParameterReader * paraRdr = new ParameterReader;
	paraRdr->readFromFile("./parameters.dat");
	paraRdr->readFromArguments(argc, argv);
	//paraRdr->echo();


	// Start timing
	Stopwatch sw;
	sw.Start();


	// Specify files containing all position-momentum information
	// from which to construct HBT correlation function
	vector<string> all_file_names;
	read_file_catalogue("./catalogue.dat", all_file_names);


	// Vector to hold all event information
	vector<EventRecord> allEvents;


	// Read in the files
	get_all_events(all_file_names, allEvents);


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


	// Check single-particle spectra first
	HBT_event_ensemble.Compute_spectra();


	// Numerator and denominator in definition
	// of correlation function
	HBT_event_ensemble.Compute_numerator();
	HBT_event_ensemble.Compute_denominator();


	// Compute correlation function itself
	HBT_event_ensemble.Compute_correlation_function();


	/*
	// Output correlation function
	HBT_event_ensemble.Output_correlation_function();
	*/


	// Print out run-time
	sw.Stop();
	//cout 	<< "Finished everything in "
	//		<< sw.printTime() << " seconds." << endl;


	// Wrap it up!
	return (0);
}

//End of file
