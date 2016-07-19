#include <string>
#include <fstream>
#include <ctime>
#include <vector>
#include <algorithm>
#include <iostream>

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

#include "lib.h"

#define COLUMN_INDEX_TO_SORT_BY	1	// duh

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// lib.h (original)
////////////////////////////////////////////////////////////////////////////////

bool fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile;
}

std::string get_selfpath() {
    char buff[PATH_MAX];
    ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff)-1);
    if (len != -1) {
      buff[len] = '\0';
      std::string temp = std::string(buff);
      unsigned found = temp.find_last_of("/\\");
      return temp.substr(0,found);
    } else {
     return "ERROR";
    }
}

int get_folder_index (string& str)
{
  unsigned found = str.find_last_of("-");
  return atoi((str.substr(found+1)).c_str());
}

////////////////////////////////////////////////////////////////////////////////
// sorter.h
////////////////////////////////////////////////////////////////////////////////

/*template <typename T>
vector<size_t> ordered(vector<T> values, int lt_or_gt)
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
vector<size_t> partial_ordered(vector<T> const& values, int cutoff, int lt_or_gt)
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
*/

/*USAGE: debugger(__LINE__, __FILE__)*/
void debugger(int cln, const char* cfn)
{
	cerr << "You made it to " << cfn << ":" << cln << "!" << endl;
	return;
}

void print_now()
{
	time_t ctt = time(0);
	cerr << asctime(localtime(&ctt)) << endl;
}

////////////////////////////////////////////////////////////////////////////////
// stringsplit.h
////////////////////////////////////////////////////////////////////////////////

// split: receives a char delimiter; returns a vector of strings
// By default ignores repeated delimiters, unless argument rep == 1.
vector<string>& splitstring::split(char delim, int rep)
{
    if (!flds.empty()) flds.clear();  // empty vector if necessary
    string work = data();
    string buf = "";
    int i = 0;
    while (i < work.length()) {
        if (work[i] != delim)
            buf += work[i];
        else if (rep == 1) {
            flds.push_back(buf);
            buf = "";
        } else if (buf.length() > 0) {
            flds.push_back(buf);
            buf = "";
        }
        i++;
    }
    if (!buf.empty())
        flds.push_back(buf);
    return flds;
}

//End of file
