#ifndef GENERATE_PROCESSING_RECORD_H
#define GENERATE_PROCESSING_RECORD_H

#include <string>
#include <sstream>
#include <time.h>

using namespace std;

#include "plumberglib.h"
#include "parameters.h"

string truestring = "true";
string falsestring = "false";

inline string return_boolean_string(bool test){return (test ? truestring : falsestring);}

void initialize_PRfile(string currentworkingdirectory, string PRfilename = "Processing_record.txt")
{
	ostringstream filename_stream;
	filename_stream << currentworkingdirectory << "/" << PRfilename;
	ofstream output(filename_stream.str().c_str(), ios::app);

	string lambdaflagstring = lambdaflag ? truestring : falsestring;

	time_t now = time(0);
	tm *ltm = localtime(&now);

	char* dt = ctime(&now);

	output << "/***************************************************/" << endl;
	output << "/****************" << PRfilename << "**************/" << endl;
	output << "/***************************************************/" << endl;

	output << "Beginning timestamp: " << dt << endl << endl;

	output << "Presets:" << endl;
	output << "   - USE_LAMBDA: " << return_boolean_string(USE_LAMBDA) << endl;
	output << "   - IGNORE_LONG_LIVED_RESONANCES: " << return_boolean_string(IGNORE_LONG_LIVED_RESONANCES) << endl;
	output << "   - FIT_WITH_PROJECTED_CFVALS: " << return_boolean_string(FIT_WITH_PROJECTED_CFVALS) << endl;
	output << "   - FLESH_OUT_CF: " << return_boolean_string(FLESH_OUT_CF) << endl;
	output << "   - REGULATE_CF: " << return_boolean_string(REGULATE_CF) << endl;
	output << endl;

	output << "General initializations:" << endl;
	output << "   - Spatial rapidity information:" << endl;
	output << "      --> eta_s_npts: " << eta_s_npts << endl;
	output << "      --> eta_s_i: " << eta_s_i << endl;
	output << "      --> eta_s_f: " << eta_s_f << endl;

	output << "   - Single-particle momentum information:" << endl;
	//output << "      --> n_interp_pT_pts: " << n_interp_pT_pts << endl;
	//output << "      --> n_interp_pphi_pts: " << n_interp_pphi_pts << endl;
	output << "      --> interp_pT_min: " << interp_pT_min << endl;
	output << "      --> interp_pT_max: " << interp_pT_max << endl;
	output << "      --> interp_pphi_min: " << interp_pphi_min << endl;
	output << "      --> interp_pphi_max: " << interp_pphi_max << endl;

	output << "   - Phase-space integral information:" << endl;
	output << "      --> s_npts: " << s_npts << endl;
	output << "      --> v_npts: " << v_npts << endl;
	output << "      --> zeta_npts: " << zeta_npts << endl;

	output << "   - Relative momentum information:" << endl;
	//output << "      --> qtnpts: " << qtnpts << endl;
	//output << "      --> qxnpts: " << qxnpts << endl;
	//output << "      --> qynpts: " << qynpts << endl;
	//output << "      --> qznpts: " << qznpts << endl;
	output << "      --> delta_qt: " << delta_qt << endl;
	output << "      --> delta_qx: " << delta_qx << endl;
	output << "      --> delta_qy: " << delta_qy << endl;
	output << "      --> delta_qz: " << delta_qz << endl;
	//output << "      --> init_qt: " << init_qt << endl;
	//output << "      --> init_qx: " << init_qx << endl;
	//output << "      --> init_qy: " << init_qy << endl;
	//output << "      --> init_qz: " << init_qz << endl;

	output << "   - Pair momentum information:" << endl;
	output << "      --> n_localp_T: " << n_localp_T << endl;
	output << "      --> localp_T_min: " << localp_T_min << endl;
	output << "      --> localp_T_max: " << localp_T_max << endl;
	output << "      --> n_localp_phi: " << n_localp_phi << endl;
	output << "      --> localp_phi_min: " << localp_phi_min << endl;
	output << "      --> localp_phi_max: " << localp_phi_max << endl;

	output << "   - Correlation function information:" << endl;
	output << "      --> corrfuncdim: " << corrfuncdim << endl;
	output << "      --> lambdaflag: " << lambdaflagstring << endl;

	output << "   - HBT Fourier information:" << endl;
	output << "      --> n_order: " << n_order << endl;

	output << "   - Miscellaneous information:" << endl;
	output << "      --> CWD: " << currentworkingdirectory << endl;

	output << "/***************************************************/" << endl << endl;

	output.close();

	return;
}

void checkforfiles_PRfile(string currentworkingdirectory, int folderindex, bool corrfuncsgenerated, string PRfilename = "Processing_record.txt")
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
	string corrfuncsexist = corrfuncsgenerated ? truestring : falsestring;

	output << "/***************************************************/" << endl
		<< "HBT output files (source variances method):" << endl
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
}

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
