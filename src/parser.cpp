#include "parser.h"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include "lcmsannotate.h"
#include "stroperation.h"

/* Parsing features
 * Supported:
 *  - Progenesis (tr_mz)
 *  - Mass hunter (mz@tr)
 *  - Compound discoverer (????)
 *  - XCMS  (mz/tr)
 *  - Simple space/tab separed tr m/z format (tr mz or mz tr tr\tmz or mz\ttr)
 */
void FeatureRead(std::string finput, std::vector<FEATURE> &featlst){
  std::string origline;
  std::string line;
  std::ifstream f_featlst(finput);
  if(f_featlst.is_open()){
    while(getline(f_featlst, origline)){
      origline = trim(origline);
	  line = trim(origline);
      strreplace(line, "min", "");
      strreplace(line, "m/z", "");
      strreplace(line, "n", "");
      /*removeCharsFromString(v[0], "min/z");
      removeCharsFromString(v[1], "min/z");*/

      char delim;
      enum featuretype{
        progenesis = 0,
        masshunter,
        xcms,
        space,
        semicolon
      };

      featuretype ftype = progenesis;

      if(line.find("_") != std::string::npos){ // Progenesis
        delim = '_';
        ftype = progenesis;
      }
      else if(line.find("@") != std::string::npos){ // Mass hunter
        delim = '@';
        ftype = masshunter;
      }
      else if(line.find("/") != std::string::npos){ // XCMS
        delim = '/';
        ftype = xcms;
      }
      else if(line.find(" ") != std::string::npos){ // simple space
        delim = ' ';
        ftype = space;
      }
      else if(line.find(";") != std::string::npos){ // semicolon
        delim = ' ';
        ftype = semicolon;
      }
      else{
        std::cerr << ">> [DynMetId ERROR] -  m/z tr list not supported. Please contact Giuseppe Marco Randazzo <gmrandazzo@gmail.com> <<" << std::endl;
        return;
      }

      std::vector<std::string> v = strsplit(line, delim);
      if(v.size() == 2){
        std::string tr;
        std::string mass;
        if(ftype == featuretype::progenesis){
          tr = v[0];
          mass = v[1];
        }
        else if(ftype == featuretype::masshunter){
          tr = v[1];
          mass = v[0];
        }
        else if(ftype == featuretype::xcms){
          std::ostringstream strs;
          strs << stod_(v[1])/60.f;
          tr = strs.str();
          mass = v[0];
        }
        else if(ftype == featuretype::space){
          tr = v[0];
          mass = v[1];
        }
        else{ //(ftype == featuretype::semicolon)
          tr = v[0];
          mass = v[1];
        }
        featlst.push_back(FEATURE(mass, tr, origline));
      }
      else{
        //SKIPP
        continue;
      }

      /*if(v.size() == 2){
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
      }*/
    }
    f_featlst.close();
  }
  else std::cerr << "File not found!" << std::endl;
}

void AdductRead(std::string finput, std::vector<ADDUCT> &adductlst){
  std::string line;
  std::ifstream faddlst(finput);
  if(faddlst.is_open()){
    while(getline(faddlst, line)){
      std::vector<std::string> v = strsplit(trim(line), ';');
      if(v.size() == 3)
        adductlst.push_back(ADDUCT(v[0], v[1], v[2]));
      else
        continue;
    }
    faddlst.close();
  }
  else std::cerr << ">> Unable to open adduct list <<" << std::endl;
  adductlst.push_back(ADDUCT("NEUTRAL", "0.0", "1.0"));
}
