#ifndef LIB_H
#define LIB_H

#include <string>
#include <fstream>
#include <ctime>
#include <vector>
#include <algorithm>
#include <iostream>

#include <limits.h>

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

using namespace std;

#define COLUMN_INDEX_TO_SORT_BY	1	// duh

////////////////////////////////////////////////////////////////////////////////
// lib.h (original)
////////////////////////////////////////////////////////////////////////////////

bool fexists(const char *filename);
std::string get_selfpath();
int get_folder_index (string& str);

////////////////////////////////////////////////////////////////////////////////
// sorter.h
////////////////////////////////////////////////////////////////////////////////

//N.B. - template stuff must appear in header file (not cpp file)!!!

template <typename T>
vector<size_t> ordered(vector<T> const& values, int lt_or_gt = 0)
{
	using namespace boost::phoenix;
	using namespace boost::phoenix::arg_names;

	vector<size_t> indices(values.size());
	int i = 0;
	transform(values.begin(), values.end(), indices.begin(), ref(i)++);
	if (lt_or_gt == 0)
		sort(indices.begin(), indices.end(), ref(values)[arg1] < ref(values)[arg2]);
	else
		sort(indices.begin(), indices.end(), ref(values)[arg1] > ref(values)[arg2]);
	return indices;
}

template <typename T>
vector<size_t> partial_ordered(vector<T> const& values, int cutoff, int lt_or_gt = 0)
{
	using namespace boost::phoenix;
	using namespace boost::phoenix::arg_names;

	vector<size_t> indices(values.size());
	int i = 0;
	transform(values.begin(), values.end(), indices.begin(), ref(i)++);
	if (lt_or_gt == 0)
		partial_sort(indices.begin(), indices.begin() + cutoff, indices.end(), ref(values)[arg1] < ref(values)[arg2]);
	else
		partial_sort(indices.begin(), indices.begin() + cutoff, indices.end(), ref(values)[arg1] > ref(values)[arg2]);
	return indices;
}


template <typename T>
bool cmp(const vector<T>& a, const vector<T>& b)
{
    return a[COLUMN_INDEX_TO_SORT_BY] < b[COLUMN_INDEX_TO_SORT_BY];
}

template <typename T>
void sort_by_column(vector< vector<T> > * values)
{
	std::stable_sort((*values).begin(), (*values).end(), cmp<T>);

	return;
}


/*USAGE: debugger(__LINE__, __FILE__);*/
void debugger(int cln, const char* cfn);
void print_now();

////////////////////////////////////////////////////////////////////////////////
// stringsplit.h
////////////////////////////////////////////////////////////////////////////////

class splitstring : public string {
    vector<string> flds;
public:
    splitstring(char *s) : string(s) { };
    vector<string>& split(char delim, int rep=0);
};

////////////////////////////////////////////////////////////////////////////////
// misc. string converter
////////////////////////////////////////////////////////////////////////////////

namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}


#endif
