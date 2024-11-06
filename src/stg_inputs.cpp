#include "stg_inputs.h"
#include <fstream>
#include <iostream>
#include <string>

// Define the inputs
double STGInputs::delta0;
double STGInputs::ub0;
double STGInputs::lx, STGInputs::ly, STGInputs::lz;
double STGInputs::xinflow, STGInputs::yinflow, STGInputs::zinflow;
bool STGInputs::lGenerateSEM;
double STGInputs::xlinf;
double STGInputs::qinfd, STGInputs::rhoinf;

void STGInputs::read_input() {
  std::fstream f;
  std::string str;
  std::string fname{"stg.inp"};
  f.open(fname, std::ios::in);
  if (f.is_open()) {
    while (getline(f, str)) {
      if (str.front() == '#')
        continue;
      for (std::size_t i = 0; i < str.size(); ++i)
        if (str[i] == ' ') {
          str.erase(i, 1);
          --i;
        }
      std::size_t found = str.find("=");
      std::string str_key = str.substr(0, found);
      std::string str_value = str.substr(found + 1, str.size());
      // if (str_key.compare("N")== 0) N = std::stoi(str_value);i
      if (str_key.compare("delta0") == 0)
        STGInputs::delta0 = std::stold(str_value);
      // if (str_key.compare("sigma")== 0) sigma = std::stod(str_value);
      if (str_key.compare("ub0") == 0)
        STGInputs::ub0 = std::stod(str_value);
      if (str_key.compare("lx") == 0)
        STGInputs::lx = std::stod(str_value);
      if (str_key.compare("ly") == 0)
        STGInputs::ly = std::stod(str_value);
      if (str_key.compare("lz") == 0)
        lz = std::stod(str_value);
      if (str_key.compare("xinflow") == 0)
        STGInputs::xinflow = std::stod(str_value);
      if (str_key.compare("yinflow") == 0)
        STGInputs::yinflow = std::stod(str_value);
      if (str_key.compare("zinflow") == 0)
        STGInputs::zinflow = std::stod(str_value);
      if (str_key.compare("lGenerateSEM") == 0)
        STGInputs::lGenerateSEM = (str_value.compare("true") == 0);
      // if (str_key.compare("lVSTG")== 0) lVSTG = (str_value.compare("true") ==
      STGInputs::lGenerateSEM = true; 
      // 0);
    }
    f.close();
  }
}
