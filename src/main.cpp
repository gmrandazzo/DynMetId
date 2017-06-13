/* main.cpp (DynMetId)
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

#include "version.h"
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

//Data structure for validation
struct VOBJ{
  VOBJ(std::string name_, std::string trexp_, std::string trpred_): name(name_), trexp(trexp_), trpred(trpred_){}
  std::string name;
  std::string trexp;
  std::string trpred;
};

void PrintRes(std::string adductname, std::vector<std::string> r)
{
  //std::cout << "___________FOUND___________ " << std::endl;
  for(size_t i = 0; i < r.size(); i++)
    std::cout << "Adduct: " << adductname << ";" << r[i] << std::endl;
  //std::cout << "___________________________ " << std::endl;
}

/* convert the identification result string into a json object */
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
  if(argc >= 17 && argc <= 18){
    // Identify metabolites
    LCMSAnnotate *lcmsann = new LCMSAnnotate;

    //Variable definitions
    std::vector<ADDUCT> adductlst;
    std::vector<FEATURE> featlst;
    std::string inpstr;
    std::vector<std::string> r;
    std::string line;

    std::string mass_parameters = format("%s", argv[8]);
    std::string chromid_parameters = format("emperror: %s prederror: %s", argv[8], argv[9], argv[10]);

    std::string chromatographic_parameters;
    chromatographic_parameters = format("init: %s final: %s tg: %s flow: %s vm: %s vd: %s", argv[11], argv[12], argv[13], argv[14], argv[15], argv[16]);

    //MySQL parameters needed to init the database.
    lcmsann->init(argv[1], argv[2], argv[3], argv[4], argv[5]);

    if(argc == 18){
      // chromatographic_parameters needed!
      // Correct retention time shift
      lcmsann->setRTLinearCorrection(argv[17], chromatographic_parameters);
    }

    // Metabolomics parameters are treated here...
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
    else std::cout << ">> Unable to open file m/z tr list <<" << std::endl;

    std::ifstream faddlst(argv[7]);
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
    else std::cout << ">> Unable to open adduct list <<" << std::endl;

    // Now for each feature search each adduct by running the standard lcmsannotate query.
    std::cout << "[" << std::endl;
    for(size_t j = 0; j < featlst.size(); j++){
      std::vector<std::string> jsonlst;
      for(size_t i = 0; i < adductlst.size(); i++){
        // Example "mass: 347.2219 error: 25ppm add: 1.0079; tr: 9.05 error: 5% init: 5 final: 95 tg: 14 flow: 0.3 vm: 0.3099 vd: 0.375";
        inpstr = format("ms: %s error: %s add: %f; tr: %s %s %s",
                        featlst[j].mass.c_str(),
                        mass_parameters.c_str(),
                        adductlst[i].ms,
                        featlst[j].tr.c_str(),
                        chromid_parameters.c_str(),
                        chromatographic_parameters.c_str());

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
    std::cout << "]" << std::endl;
    delete lcmsann;
  }
  else if(argc >= 11 && argc <= 13){
    //Export the whole database in a CSV metascope file
    LCMSAnnotate *lcmsann = new LCMSAnnotate;

    //Variable definitions
    std::vector<std::string> r;

    std::string chromatographic_parameters;
    chromatographic_parameters = format("init: %s final: %s tg: %s flow: %s vm: %s vd: %s", argv[6], argv[7], argv[8], argv[9], argv[10], argv[11]);

    //MySQL parameters to init the database.
    lcmsann->init(argv[1], argv[2], argv[3], argv[4], argv[5]);

    if(argc == 13){
      // chromatographic_parameters needed!
      // Correct retention time shift
      lcmsann->setRTLinearCorrection(argv[12], chromatographic_parameters);
    }

    r = lcmsann->db2MetaScope(chromatographic_parameters);

    for(size_t i = 0; i < r.size(); i++){
      std::cout << r[i] << std::endl;
    }

    r.clear();
    delete lcmsann;
  }
  else if(strcmp(argv[6], "validate") == 0 && argc >= 14){
    //Validate the retention time prediction
    LCMSAnnotate *lcmsann = new LCMSAnnotate;

    //Variable definitions
    std::string chromatographic_parameters;
    chromatographic_parameters = format("init: %s final: %s tg: %s flow: %s vm: %s vd: %s", argv[8], argv[9], argv[10], argv[11], argv[12], argv[13]);

    //MySQL parameters to init the database.
    lcmsann->init(argv[1], argv[2], argv[3], argv[4], argv[5]);

    if(argc == 15){
      // chromatographic_parameters needed!
      // Correct retention time shift
      lcmsann->setRTLinearCorrection(argv[14], chromatographic_parameters);
    }

    std::vector<VOBJ> vobj;

    std::ifstream trlst(argv[7]);
    std::string line;

    if(trlst.is_open()){
      while(getline(trlst, line)){
        std::vector<std::string> v = strsplit(trim(line), ';');
        if(v.size() == 2){
          std::stringstream ss;
          ss << "name " << v[0] << "; tr " << v[1] << " emperror: 100.0 prederror: 100.0 " << chromatographic_parameters;
          std::vector<std::string> res = lcmsann->find(ss.str());
          for(size_t i = 0; i < res.size(); i++){
            std::vector<std::string> a = strsplit(res[i], ';');
            std::string name = trim(strsplit(a[0], ':').back());
            if(name.compare(v[0]) == 0){
              vobj.push_back(VOBJ(v[0], v[1], strsplit(a[3], ':')[1]));
              break;
            }
            else{
              continue;
            }
          }

          /*if(res.size() == 1){
            std::vector<std::string> a = strsplit(res[0], ';');
            std::vector<std::string> b = strsplit(a[3], ':');
            x.push_back(stod_(v[1]));
            y.push_back(stod_(b[1]));
          }
          else{
            continue;
          }*/
        }
      }
      trlst.close();
    }
    else{
      std::cout << "Unable to open compound retention time list" << std::endl;
      return 0;
    }

    //JSON output
    std::cout << "[";
    for(size_t i = 0; i < vobj.size()-1; i++){
      std::cout << "{";
      std::cout << "\"name\": "  << "\"" << vobj[i].name << "\",";
      std::cout << "\"trexp\": "  << "\"" << vobj[i].trexp << "\",";
      std::cout << "\"trpred\": " << "\"" << vobj[i].trpred << "\"";
      std::cout << "},";
    }
    size_t last = vobj.size()-1;
    std::cout << "{";
    std::cout << "\"name\": "  << "\"" << vobj[last].name << "\",";
    std::cout << "\"trexp\": "  << "\"" << vobj[last].trexp << "\",";
    std::cout << "\"trpred\": " << "\"" << vobj[last].trpred << "\"";
    std::cout << "}";
    std::cout << "]";
    delete lcmsann;
  }
  else{
    std::cout << format("DynMetId %d.%d.%d [Dynamic Metabolite Identification tool]", dynmetid_major, dynmetid_minor, dynmetid_patch) << std::endl;
    std::cout << format("Written by: Giuseppe Marco Randazzo <gmrandazzo@gmail.com>") << std::endl;
    std::cout << format("Software distributed under GPLv3 license\n", argv[0]) << std::endl;
    std::cout << format("Usage:") << std::endl;
    std::cout << format("Do annotation             :       %s [mysql parameters] [metabolomics parameters] [identification parameters] [chromatographic parameters] [optional parameters]", argv[0]) << std::endl;
    std::cout << format("Export the whole database :       %s [mysql parameters] [chromatographic parameters] [optional parameters]", argv[0]) << std::endl;
    std::cout << format("Validate tr prediction    :       %s [mysql parameters] validate <compounds and retention time list> [chromatographic parameters] [optional parameters]", argv[0]) << std::endl;
    std::cout << format("       mysql parameters           : <host> <username> <password> <db name> <table name>") << std::endl;
    std::cout << format("       identification parameters  : <ppm error> <tr empirical error> <tr predicted error>") << std::endl;
    std::cout << format("       metabolomics parameters    : <feature list> <adduct list>    i.e: \"9.12_306.2841m/z ...\" and \"1.0072764649920167;M+H ...\" ") << std::endl;
    std::cout << format("       chromatographic parameters : <gradient start> <gradient stop> <time gradient> <flow rate> <dead volume> <dwell volume> <correction shift list>") << std::endl;
    std::cout << format("       optional parameters        : <Correction retention time shift list>\n") << std::endl;
    //std::cout << format("\nUsage: %s <host> <username> <password> <db name> <table name>  <feature list> <mass error> <adduct mass list> <tr empirical error> <tr predicted error> <gradient start> <gradient stop> <time gradient> <flow rate> <dead volume> <dwell volume>\n", argv[0]) << std::endl;
    std::cout << format("\nExample:\n") << std::endl;
    std::cout << format("%s localhost User1 0000 stddb steroids list_of_feature.txt list_of_adducts.txt 5.0 2 4 5 95 14 0.3 0.983 0.375", argv[0]) << std::endl;
    std::cout << format("    where:") << std::endl;
    std::cout << format("    mysql parameters                   = localhost User1 0000 stddb steroid") << std::endl;
    std::cout << format("    metabolomics parameters            = list_of_feature.txt list_of_adducts.txt") << std::endl;
    std::cout << format("    identification parameters          = 5 2 4") << std::endl;
    std::cout << format("    chromatographic parameters         = 5 95 14 0.3 0.983 0.375") << std::endl;
    std::cout << std::endl;
  }
  return 0;
}
