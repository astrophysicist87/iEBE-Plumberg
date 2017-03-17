#ifndef GENERATE_PROCESSING_RECORD_H
#define GENERATE_PROCESSING_RECORD_H

#include <string>
#include <sstream>
#include <time.h>

using namespace std;

#include "lib.h"
#include "parameters.h"
#include "ParameterReader.h"

string truestring = "true";
string falsestring = "false";

inline string return_boolean_string(bool test){return (test ? truestring : falsestring);}

void initialize_PRfile(ParameterReader* paraRdr, string currentworkingdirectory, string PRfilename = "Processing_record.txt")
{
	ostringstream filename_stream;
	filename_stream << currentworkingdirectory << "/" << PRfilename;
	ofstream output(filename_stream.str().c_str(), ios::app);

	stringstream out;

	time_t now = time(0);
	tm *ltm = localtime(&now);

	char* dt = ctime(&now);

	out << setfill('0') << setw(2) << 1 + ltm->tm_mon << "/" << setw(2) << ltm->tm_mday << "/" << setw(4) << 1900 + ltm->tm_year;

	string date = out.str();


	output << "/***************************************************/" << endl;
	output << "/****************" << PRfilename << "**************/" << endl;
	output << "/***************************************************/" << endl;

	output << "Beginning timestamp: " << dt << endl << endl;

	output << "Presets:" << endl;
	output << "   - GROUPING_PARTICLES: " << return_boolean_string(paraRdr->getVal("grouping_particles")) << endl;
	output << "   - PARTICLE_DIFF_TOLERANCE: " << paraRdr->getVal("particle_diff_tolerance") << endl;
	//output << "   - TRUNCATE_COOPER_FRYE: " << return_boolean_string(paraRdr->getVal("truncate_cooper_frye")) << endl;
	//output << "   - SPACETIME_MOMENTS_ONLY: " << return_boolean_string(paraRdr->getVal("spacetime_moments_only")) << endl;
	//output << "   - INCLUDE_SOURCE_VARIANCES: " << return_boolean_string(paraRdr->getVal("include_source_variances")) << endl;
	output << "   - INCLUDE_DELTA_F: " << return_boolean_string(paraRdr->getVal("include_delta_f")) << endl;
	output << "   - INCLUDE_BULK_PI: " << return_boolean_string(paraRdr->getVal("include_bulk_pi")) << endl;
	output << "   - IGNORE_LONG_LIVED_RESONANCES: " << return_boolean_string(paraRdr->getVal("ignore_long_lived_resonances")) << endl;
	output << "   - USE_PLANE_PSI_ORDER: " << paraRdr->getVal("use_plane_psi_order") << endl;

	output << endl;

	output << "General initializations:" << endl;
	output << "   - Spatial rapidity information:" << endl;
	output << "      --> eta_s_npts: " << eta_s_npts << endl;
	output << "      --> eta_s_i: " << eta_s_i << endl;
	output << "      --> eta_s_f: " << eta_s_f << endl;

	output << "   - Single particle momentum information:" << endl;
	output << "      --> n_pT_pts: " << paraRdr->getVal("SV_npT") << endl;
	output << "      --> n_pphi_pts: " << paraRdr->getVal("SV_npphi") << endl;
	output << "      --> SP_pT_min: " << SP_pT_min << endl;
	output << "      --> SP_pT_max: " << SP_pT_max << endl;
	output << "      --> SP_pphi_min: " << SP_pphi_min << endl;
	output << "      --> SP_pphi_max: " << SP_pphi_max << endl;

	output << "   - Pair momentum information:" << endl;
	output << "      --> nKT: " << paraRdr->getVal("nKT") << endl;
	output << "      --> KT_min: " << paraRdr->getVal("KTmin") << endl;
	output << "      --> KT_max: " << paraRdr->getVal("KTmax") << endl;
	output << "      --> nKphi: " << paraRdr->getVal("nKphi") << endl;
	output << "      --> Kphi_min: " << Kphi_min << endl;
	output << "      --> Kphi_max: " << Kphi_max << endl;

	output << "   - Phase-space integral information:" << endl;
	output << "      --> n_s_pts: " << n_s_pts << endl;
	output << "      --> n_v_pts: " << n_v_pts << endl;
	output << "      --> n_zeta_pts: " << n_zeta_pts << endl;

	output << "   - HBT Fourier information:" << endl;
	output << "      --> n_order: " << paraRdr->getVal("n_order") << endl;

	output << "   - Miscellaneous information:" << endl;
	output << "      --> CWD: " << currentworkingdirectory << endl;
	output << "      --> tol: " << paraRdr->getVal("tolerance") << endl;
	output << "      --> flagneg: " << paraRdr->getVal("flag_negative_S") << endl;
	if ( paraRdr->getVal("ignore_long_lived_resonances") )
		output << "      --> max_lifetime (fm/c): " << paraRdr->getVal("max_lifetime") << endl;
	else
		output << "      --> max_lifetime (fm/c): 10000000000.0" << endl;

	output << "/***************************************************/" << endl << endl;

	output.close();

	return;
}

/*void checkforfiles_PRfile(string currentworkingdirectory, int folderindex, bool corrfuncsgenerated, string PRfilename = "Processing_record.txt")
{
	ostringstream filename_stream;
	filename_stream << currentworkingdirectory << "/" << PRfilename;
	ofstream output(filename_stream.str().c_str(), ios::app);

	stringstream HBTSVfilename, HBTcfsSVfilename, planepsifilename;
	stringstream HBTGFfilename, HBTcfsGFfilename;
	HBTSVfilename << currentworkingdirectory << "/HBTradii_ev" << folderindex << ".dat" << endl;
	HBTcfsSVfilename << currentworkingdirectory << "/HBTradii_cfs_ev" << folderindex << ".dat" << endl;
	HBTGFfilename << currentworkingdirectory << "/HBTradii_GF_ev" << folderindex << ".dat" << endl;
	HBTcfsGFfilename << currentworkingdirectory << "/HBTradii_GF_cfs_ev" << folderindex << ".dat" << endl;
	planepsifilename << currentworkingdirectory << "/plane_psi_ev" << folderindex << ".dat" << endl;

	string HBTSVexists = fexists(HBTSVfilename.str().c_str()) ? truestring : falsestring;
	string HBTcfsSVexists = fexists(HBTcfsSVfilename.str().c_str()) ? truestring : falsestring;
	string HBTGFexists = fexists(HBTGFfilename.str().c_str()) ? truestring : falsestring;
	string HBTcfsGFexists = fexists(HBTcfsGFfilename.str().c_str()) ? truestring : falsestring;
	string planepsiexists = fexists(planepsifilename.str().c_str()) ? truestring : falsestring;
	string corrfuncsexist = corrfuncsgenerated ? truestring : falsestring;*/

	//output << "/***************************************************/" << endl
	/*	<< "HBT output files (source variances method):" << endl
		<< "      --> Generated " << HBTSVfilename.str() << ": " << HBTSVexists << endl
		<< "      --> Generated " << HBTcfsSVfilename.str() << ": " << HBTcfsSVexists << endl
		<< "      --> Generated " << planepsifilename.str() << ": " << planepsiexists << endl;

	output << "HBT output files (Gaussian fit method):" << endl
		<< "      --> Generated " << HBTGFfilename.str() << ": " << HBTGFexists << endl
		<< "      --> Generated " << HBTcfsGFfilename.str() << ": " << HBTcfsGFexists << endl;

	output << "Correlation function output files:" << endl
		<< "      --> Generated correlation functions in corrfuncs_ev" << folderindex << ".zip: " << corrfuncsexist << endl;

	output.close();

	return;
}*/

void finalize_PRfile(string currentworkingdirectory, string PRfilename = "Processing_record.txt")
{
	ostringstream filename_stream;
	filename_stream << currentworkingdirectory << "/" << PRfilename;
	ofstream output(filename_stream.str().c_str(), ios::app);

	stringstream out;

	time_t now = time(0);
	tm *ltm = localtime(&now);

	char* dt = ctime(&now);

	out << setfill('0') << setw(2) << 1 + ltm->tm_mon << "/" << setw(2) << ltm->tm_mday << "/" << setw(4) << 1900 + ltm->tm_year;

	string date = out.str();

	output << "Ending timestamp: " << dt << endl;

	output << "/***************************************************/" << endl;
	output << "/*********End of processing for " << date << "**********/" << endl;
	output << "/***************************************************/" << endl;

	output.close();

	return;
}

#endif
