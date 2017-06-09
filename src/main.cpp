/* main.cpp (MySQLDynMetId)
*
* Copyright (C) <2017>  Giuseppe Marco Randazzo
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <algorithm>
#include <iostream>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <cctype>
#include <cstdarg>

#include "lcmsannotate.h"
#include "stroperation.h"

// Datastructure for adducts
struct ADDUCT{
  ADDUCT(std::string name_, std::string ms_): name(name_), ms(atof(ms_.c_str())){}
  std::string name;
  double ms;
};

//Data structure for feature
struct FEATURE{
  FEATURE(std::string mass_, std::string tr_, std::string origname_): mass(mass_), tr(tr_), origname(origname_){}
  std::string mass;
  std::string tr;
  std::string origname;
};

//missing string printf
//this is safe and convenient but not exactly efficient
inline std::string format(const char* fmt, ...){
    int size = 512;
    char* buffer = 0;
    buffer = new char[size];
    va_list vl;
    va_start(vl, fmt);
    int nsize = vsnprintf(buffer, size, fmt, vl);
    if(size<=nsize){ //fail delete buffer and try again
        delete[] buffer;
        buffer = 0;
        buffer = new char[nsize+1]; //+1 for /0
        nsize = vsnprintf(buffer, size, fmt, vl);
    }
    std::string ret(buffer);
    va_end(vl);
    delete[] buffer;
    return ret;
}

void PrintRes(std::string adductname, std::vector<std::string> r)
{
  //std::cout << "___________FOUND___________ " << std::endl;
  for(size_t i = 0; i < r.size(); i++)
    std::cout << "Adduct: " << adductname << ";" << r[i] << std::endl;
  //std::cout << "___________________________ " << std::endl;
}

std::string Annotation2JSON(std::string adductname, std::string str){
  std::string res;
  std::vector<std::string> v = strsplit(str, ';');
  res.append("{\n");
  res.append((std::string)"adduct: " + (std::string)"\"" + adductname + (std::string)"\",\n");
  for(size_t i = 0; i < v.size()-1; i++){
    if(v[i].find("link") != std::string::npos){
      res.append((std::string)"link: " +(std::string) "\"" + trim(purgestring(v[i], "link:")) + (std::string)"\",\n" );
    }
    else{
      std::vector<std::string> a = strsplit(v[i], ':');
      res.append(trim(a[0]) + (std::string)": " + (std::string)"\"" + trim(a[1]) + (std::string)"\",\n");
    }
  }

  if(v[v.size()-1].find("link") != std::string::npos){
    res.append("link: " + (std::string)"\"" + trim(purgestring(v[v.size()-1], (std::string)"link:")) + "\"\n" );
  }
  else{
    std::vector<std::string> a = strsplit(v[v.size()-1], ':');
    res.append(trim(a[0]) + (std::string)": " + (std::string)"\"" + trim(a[1]) + (std::string)"\"\n");
  }
  res.append("}");
  return res;
}

int main(int argc, char **argv)
{
  if(argc >= 16){
    /*Initialize the database*/

    LCMSAnnotate *lcmsann = new LCMSAnnotate;
    lcmsann->init(argv[1], argv[2], argv[3], argv[4], argv[5]);

    if(argc == 17){
      //lcmsann->setRTLinearAligner(f_trcorr, qline);
      //db->setRTLinearAligner(std::stof(argv[12]), std::stof(argv[13]));
    }

    std::vector<ADDUCT> adductlst;
    std::vector<FEATURE> featlst;
    std::string inpstr;
    std::vector<std::string> r;
    std::string line;

    std::ifstream f_featlst(argv[6]);
    if(f_featlst.is_open()){
      while(getline(f_featlst, line)){
        std::vector<std::string> v = strsplit(line, '_'); // Progenesis support
        if(v.size() == 2){
          removeCharsFromString(v[1], "m/zn");
          featlst.push_back(FEATURE(v[1], v[0], line));
        }
        else{ //TODO: add xcms support
          v = strsplit(line, '@'); // Mass hunter support
          if(v.size() == 2){
            removeCharsFromString(v[1], "m/zn");
            featlst.push_back(FEATURE(v[0], v[1], line));
          }
          else{
            continue;
          }
        }
      }
      f_featlst.close();
    }
    else std::cout << "Unable to open file m/z tr list" << std::endl;

    std::ifstream faddlst(argv[8]);
    if(faddlst.is_open()){
      while(getline(faddlst, line)){
        std::vector<std::string> v = strsplit(trim(line), ';');
        if(v.size() == 2)
          adductlst.push_back(ADDUCT(v[1], v[0]));
        else
          continue;
      }
      faddlst.close();
    }
    else std::cout << "Unable to open adduct list" << std::endl;

    // Now for each feature search each adduct by running the standard lcmsannotate query.
    std::cout << "[" << std::endl;
    for(size_t j = 0; j < featlst.size(); j++){
      std::vector<std::string> jsonlst;
      for(size_t i = 0; i < adductlst.size(); i++){
        //inpstr = "mass: 347.2219 error: 25ppm add: 1.0079; tr: 9.05 error: 5% init: 5 final: 95 tg: 14 flow: 0.3 vm: 0.3099 vd: 0.375";
        inpstr = format("ms: %s error: %sppm add: %f; tr: %s emperror: %s%% prederror: %s%%  init: %s final: %s tg: %s flow: %s vm: %s vd: %s", featlst[j].mass.c_str(), argv[7], adductlst[i].ms, featlst[j].tr.c_str(), argv[9], argv[10], argv[11], argv[12], argv[13], argv[14], argv[15], argv[16]);
        //std::cout << "Searching for: " << inpstr << "\n" << std::endl;
        r = lcmsann->find(inpstr);
        if(r.size() > 0){
          for(size_t k = 0; k < r.size(); k++){
            jsonlst.push_back(Annotation2JSON(adductlst[i].name, r[k]));
          }
        }
        r.clear();
      }

      if(jsonlst.size() > 0){
        std::cout << "{" << std::endl;
        std::cout << "key:" << "\"" << featlst[j].origname << "\"," << std::endl;
        std::cout << "mass:" << "\"" << featlst[j].mass << "\"," << std::endl;
        std::cout << "tr:" << "\"" << featlst[j].tr << "\"," << std::endl;
        std::cout << "values: [" << std::endl;
        for(size_t i = 0; i < jsonlst.size()-1; i++){
          std::cout << jsonlst[i] << std::endl;
          std::cout << "," << std::endl;
        }
        std::cout << jsonlst[jsonlst.size()-1] << std::endl;
        std::cout << "]" << std::endl;

        if(j == featlst.size()-1){
          std::cout << "}" << std::endl;
        }
        else{
          std::cout << "}," << std::endl;
        }
      }
    }
    delete lcmsann;
    std::cout << "]" << std::endl;
  }
  else{
    std::cout << format("\nUsage: %s <host> <username> <password> <db name> <table name>  <feature list> <mass error> <adduct mass list> <tr empirical error> <tr predicted error> <gradient start> <gradient stop> <time gradient> <flow rate> <dead volume> <dwell volume>\n", argv[0]) << std::endl;
    std::cout << format("\nexample") << std::endl;
    std::cout << format("%s localhost User1 0000 stddb steroids list_of_feature.txt  5 list_of_adducts.txt 2% 5 95 14 0.3 0.983 0.375\n", argv[0]) << std::endl;
    std::cout << argv[0] << " was written by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>\n" << std::endl;
  }
  return 0;
}
