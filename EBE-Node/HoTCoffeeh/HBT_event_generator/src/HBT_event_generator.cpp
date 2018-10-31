#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>

#include "HBT_event_generator.h"
#include "Stopwatch.h"

using namespace std;

HBT_event_generator::HBT_event_generator(
	ParameterReader * paraRdr_in,
	vector<EventRecord> &allEvents,
	ostream& out_stream /*= std::cout*/,
	ostream& err_stream /*= std::cerr*/)
{
	paraRdr = paraRdr_in;

	//FIT_WITH_PROJECTED_CFVALS
	//	= paraRdr->getVal("fit_with_projected_cfvals");
	//FLESH_OUT_CF
	//	= paraRdr->getVal("flesh_out_cf");

	//Set header info
	//Define various grid sizes
	// - SP momentum points at which to evaluate correlation function
	n_pT_pts 		= paraRdr->getVal("CF_npT");
	n_pphi_pts 		= paraRdr->getVal("CF_npphi");
	n_pY_pts 		= paraRdr->getVal("CF_npY");
	// - pair momenta points at which to interpolate HBT results
	nKT 			= paraRdr->getVal("nKT");
	nKphi 			= paraRdr->getVal("nKphi");
	KT_min 			= paraRdr->getVal("KTmin");
	KT_max 			= paraRdr->getVal("KTmax");
	// - relative momentum points at which to evaluate
	//   correlation function
	qonpts 			= paraRdr->getVal("qonpts");
	qsnpts 			= paraRdr->getVal("qsnpts");
	qlnpts 			= paraRdr->getVal("qlnpts");
	// - step size in q directions
	delta_qo 		= paraRdr->getVal("delta_qo");
	delta_qs 		= paraRdr->getVal("delta_qs");
	delta_ql 		= paraRdr->getVal("delta_ql");
	// - minimum value in each q direction
	init_qo 		= -0.5*double(qonpts-1)*delta_qo;
	init_qs 		= -0.5*double(qsnpts-1)*delta_qs;
	init_ql 		= -0.5*double(qlnpts-1)*delta_ql;

	// - number of points to use when fleshing out correlation
	//   function in each direction
	new_nqopts 		= ( qonpts > 1 ) ? new_nqpts : 1;
	new_nqspts 		= ( qsnpts > 1 ) ? new_nqpts : 1;
	new_nqlpts 		= ( qlnpts > 1 ) ? new_nqpts : 1;

	qonbins = qonpts - 1;
	qsnbins = qsnpts - 1;
	qlnbins = qlnpts - 1;

	//set ofstream for output file
	out 			= out_stream;
	err 			= err_stream;
	
	return;
}

HBT_event_generator::~HBT_event_generator()
{
	//clear everything

	return;
}

//End of file
