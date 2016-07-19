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

	time_t now = time(0);
	tm *ltm = localtime(&now);

	char* dt = ctime(&now);

	output << "/***************************************************/" << endl;
	output << "/****************" << PRfilename << "**************/" << endl;
	output << "/***************************************************/" << endl;

	output << "Beginning timestamp: " << dt << endl << endl;

	output << "Presets:" << endl;
	output << "   - GROUPING_PARTICLES: " << return_boolean_string(paraRdr->getVal("grouping_particles")) << endl;
	output << "   - PARTICLE_DIFF_TOLERANCE: " << paraRdr->getVal("particle_diff_tolerance") << endl;
	output << "   - INCLUDE_DELTA_F: " << return_boolean_string(paraRdr->getVal("include_delta_f")) << endl;
	output << "   - USE_LAMBDA: " << return_boolean_string(paraRdr->getVal("use_lambda")) << endl;
	output << "   - USE_EXTRAPOLATION: " << return_boolean_string(paraRdr->getVal("use_extrapolation")) << endl;
	output << "   - COMPUTE_RESONANCE_DECAYS: " << return_boolean_string(paraRdr->getVal("compute_resonance_decays")) << endl;
	output << "   - IGNORE_LONG_LIVED_RESONANCES: " << return_boolean_string(paraRdr->getVal("ignore_long_lived_resonances")) << endl;
	output << "   - FIT_WITH_PROJECTED_CFVALS: " << return_boolean_string(paraRdr->getVal("fit_with_projected_cfvals")) << endl;
	output << "   - FLESH_OUT_CF: " << return_boolean_string(paraRdr->getVal("flesh_out_cf")) << endl;

	output << endl;

	output << "General initializations:" << endl;
	output << "   - Spatial rapidity information:" << endl;
	output << "      --> eta_s_npts: " << eta_s_npts << endl;
	output << "      --> eta_s_i: " << eta_s_i << endl;
	output << "      --> eta_s_f: " << eta_s_f << endl;

	output << "   - Single-particle momentum information:" << endl;
	output << "      --> n_interp_pT_pts: " << n_interp_pT_pts << endl;
	output << "      --> n_interp_pphi_pts: " << n_interp_pphi_pts << endl;
	output << "      --> interp_pT_min: " << interp_pT_min << endl;
	output << "      --> interp_pT_max: " << interp_pT_max << endl;
	output << "      --> interp_pphi_min: " << interp_pphi_min << endl;
	output << "      --> interp_pphi_max: " << interp_pphi_max << endl;

	output << "   - Phase-space integral information:" << endl;
	output << "      --> s_npts: " << s_npts << endl;
	output << "      --> v_npts: " << v_npts << endl;
	output << "      --> zeta_npts: " << zeta_npts << endl;

	output << "   - Relative momentum information:" << endl;
	output << "      --> qtnpts: " << qtnpts << endl;
	output << "      --> qxnpts: " << qxnpts << endl;
	output << "      --> qynpts: " << qynpts << endl;
	output << "      --> qznpts: " << qznpts << endl;
	output << "      --> delta_qt: " << delta_qt << endl;
	output << "      --> delta_qx: " << delta_qx << endl;
	output << "      --> delta_qy: " << delta_qy << endl;
	output << "      --> delta_qz: " << delta_qz << endl;
	output << "      --> init_qt: " << init_qt << endl;
	output << "      --> init_qx: " << init_qx << endl;
	output << "      --> init_qy: " << init_qy << endl;
	output << "      --> init_qz: " << init_qz << endl;

	output << "   - Pair momentum information:" << endl;
	output << "      --> n_localp_T: " << n_localp_T << endl;
	output << "      --> localp_T_min: " << localp_T_min << endl;
	output << "      --> localp_T_max: " << localp_T_max << endl;
	output << "      --> n_localp_phi: " << n_localp_phi << endl;
	output << "      --> localp_phi_min: " << localp_phi_min << endl;
	output << "      --> localp_phi_max: " << localp_phi_max << endl;

	output << "   - HBT Fourier information:" << endl;
	output << "      --> n_order: " << paraRdr->getVal("n_order") << endl;
	output << endl;
	output << "   - Miscellaneous information:" << endl;
	output << "      --> CWD: " << currentworkingdirectory << endl;
	output << "      --> tol: " << paraRdr->getVal("tolerance") << endl;
	output << "      --> flagneg: " << paraRdr->getVal("flag_negative_S") << endl;
	if ( paraRdr->getVal("do_all_decay_channels") )
		output << "      --> max_lifetime (fm/c): " << paraRdr->getVal("max_lifetime") << endl;
	else
		output << "      --> max_lifetime (fm/c): 10000000000.0" << endl;

	output << "/***************************************************/" << endl << endl;

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
