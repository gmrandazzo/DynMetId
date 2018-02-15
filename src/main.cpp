/* main.cpp (DynMetId)
*
* Copyright (C) <2017>  Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
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
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <cctype>
#include <cstdarg>

#include "version.h"
#include "parser.h"
#include "lcmsannotate.h"
#include "stroperation.h"

// Compare features
bool featurecmp(const FEATURE& a, const FEATURE& b)
{
  return stod_(a.tr) < stod_(b.tr);
}

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
void Annotation2JSON(std::string adductname, std::string rstr, std::string& jsonobj, std::vector<std::string>& json_annotated_names){
  std::vector<std::string> v = strsplit(rstr, ';');
  jsonobj.append("{");
  jsonobj.append((std::string)"\"adduct\": " + (std::string)"\"" + adductname + (std::string)"\",");
  for(size_t i = 0; i < v.size()-1; i++){

    if(lower(v[i]).find("link") != std::string::npos){
      jsonobj.append((std::string)"\"link\": " +(std::string) "\"" + trim(purgestring(v[i], "link:")) + (std::string)"\"," );
    }
    else{
      std::vector<std::string> a = strsplit(v[i], ':');
      if(lower(a[0]).find("name") != std::string::npos){
        json_annotated_names.push_back(trim(a[1]));
      }
      jsonobj.append((std::string)"\""+trim(a[0]) + (std::string)"\": " + (std::string)"\"" + trim(a[1]) + (std::string)"\",");
    }
  }

  if(lower(v[v.size()-1]).find("link") != std::string::npos){
    jsonobj.append("\"link\": " + (std::string)"\"" + trim(purgestring(v[v.size()-1], (std::string)"link:")) + "\"" );
  }
  else{
    std::vector<std::string> a = strsplit(v[v.size()-1], ':');
    if(lower(a[0]).find("name") != std::string::npos){
      json_annotated_names.push_back(trim(a[1]));
    }
    jsonobj.append((std::string)"\""+trim(a[0]) + (std::string)"\": " + (std::string)"\"" + trim(a[1]) + (std::string)"\"");
  }
  jsonobj.append("}");
}

int main(int argc, char **argv)
{
  bool verbose_in_cerr = true;
  if(verbose_in_cerr) std::cerr << "Running DynMetId\n";
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
    std::string chromid_parameters = format("emperror: %s prederror: %s", argv[9], argv[10]);

    std::string chromatographic_parameters;
    chromatographic_parameters = format("init: %s final: %s tg: %s flow: %s vm: %s vd: %s", argv[11], argv[12], argv[13], argv[14], argv[15], argv[16]);

    //MySQL parameters needed to init the database.
    lcmsann->init(argv[1], argv[2], argv[3], argv[4], argv[5]);

    if(argc == 18){
      // chromatographic_parameters needed!
      // Correct retention time shift
      lcmsann->setRTLinearCorrection(argv[17], chromatographic_parameters);
    }

    /* Read features */
    FeatureRead(argv[6], featlst);
	if(verbose_in_cerr) std::cerr << "Loaded " << featlst.size() << " features\n";

    // Sort feature list
    std::sort(std::begin(featlst), std::end(featlst), featurecmp);

    /*Read adducts*/
    AdductRead(argv[7], adductlst);
	if(verbose_in_cerr) std::cerr << "Loaded " << adductlst.size() << " adducts\n";

    /* Now for each feature search each adduct by running the standard lcmsannotate query. */
    std::vector<std::string> annotated_results;
    for(size_t j = 0; j < featlst.size(); j++){
      std::vector<std::string> jsonlst;
      std::string annotated_names;
      std::vector<std::string> json_annotated_names;

      size_t nb_level2p = 0;
      size_t nb_level2 = 0;
      size_t nb_unknown = 0;

      for(size_t i = 0; i < adductlst.size(); i++){
        /* Query Example "mass: 347.2219 error: 25ppm add: 1.0079 1; tr: 9.05 error: 5% init: 5 final: 95 tg: 14 flow: 0.3 vm: 0.3099 vd: 0.375"; */
        inpstr = format("ms: %s %s error: %s add: %.8f %.2f; tr: %s %s %s",
                        featlst[j].mass.c_str(),
						featlst[j].origname.c_str(),
                        mass_parameters.c_str(),
                        adductlst[i].ms,
                        adductlst[i].mult,
                        featlst[j].tr.c_str(),
                        chromid_parameters.c_str(),
                        chromatographic_parameters.c_str());
        r = lcmsann->find(inpstr);
        if(r.size() > 0){
          for(size_t k = 0; k < r.size(); k++){
            if(lower(r[k]).find("experimental") != std::string::npos){
              nb_level2p += 1;
            }
            else if(lower(r[k]).find("predicted") != std::string::npos){
              nb_level2 += 1;
            }
            else{
              nb_unknown += 1;
            }
            std::string jsonobj;
            Annotation2JSON(adductlst[i].name, r[k], jsonobj, json_annotated_names);
            jsonlst.push_back(jsonobj);
          }
        }
        r.clear();
      }
      
      if(jsonlst.size() > 0 && json_annotated_names.size() > 0){
        std::string annotation;
        annotation += "{";
        annotation += "\"feature\": \"" + featlst[j].origname + "\",";
        annotation += "\"mass\": \"" + featlst[j].mass + "\",";
        annotation += "\"tr\": \"" + featlst[j].tr + "\",";
        annotation += "\"nb_level2p\": \"" + intost(nb_level2p) + "\",";
        annotation += "\"nb_level2\": \"" + intost(nb_level2) + "\",";
        annotation += "\"nb_unknown\": \"" + intost(nb_unknown) + "\",";
        std::string ann_names;
        for(size_t i = 0; i < json_annotated_names.size()-1; i++)
          ann_names.append((std::string)"\""+json_annotated_names[i]+(std::string)"\",");
        ann_names.append((std::string)"\""+json_annotated_names[json_annotated_names.size()-1]+(std::string)"\"");
        annotation += "\"annotated_names\": [" + ann_names + "],";
        annotation += "\"annotations\": [";
        for(size_t i = 0; i < jsonlst.size()-1; i++){
          annotation += jsonlst[i]+",";
        }
        annotation += jsonlst[jsonlst.size()-1]+"]}";
        annotated_results.push_back(annotation);
      }
    }

    if(verbose_in_cerr) std::cerr << "Found " << annotated_results.size() << " annotations\n";
    if(annotated_results.size()>0){
      std::cout << "[";
      for(size_t i = 0; i < annotated_results.size()-1; i++){
        std::cout << annotated_results[i]+",";
      }
      std::cout << annotated_results[annotated_results.size()-1];
      std::cout << "]";
    }
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
      stf::cerr << "Unable to open compound retention time list" << std::endl;
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
