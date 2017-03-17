#ifndef PLUMBERGLIB_H
#define PLUMBERGLIB_H

#include<string>
#include<fstream>
#include<ctime>

using namespace std;

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

#endif
