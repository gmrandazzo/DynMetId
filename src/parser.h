/* parser.h (DynMetId)
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

#ifndef PARSER_H
#define PARSER_H

#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

/* Datastructure for adducts
 * The adduct is defined with a neutral mass function of the monoisotopic atom present in the adducts
 * and a multiplier which is function of the number of charges and number of monoisotopic atom present in the adduct.
 * See LCMSAnnotate::find for more explaination of the variable "mult"
 */
struct ADDUCT{
  ADDUCT(std::string name_, std::string ms_, std::string mult_) : name(name_), ms(atof(ms_.c_str())), mult(atof(mult_.c_str())){}
  std::string name;
  double ms;
  double mult;
};

/* Data structure for feature */
struct FEATURE{
  FEATURE(std::string mass_, std::string tr_, std::string origname_): mass(mass_), tr(tr_), origname(origname_){}
  std::string mass;
  std::string tr;
  std::string origname;
};

void FeatureRead(std::string finput, std::vector<FEATURE> &featlst);
void AdductRead(std::string finput, std::vector<ADDUCT> &adductlst);

#endif
