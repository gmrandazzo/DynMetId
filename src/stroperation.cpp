/* stroperation.cpp (MySQLDynMetId)
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

#include "stroperation.h"
#include <algorithm>
#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <locale>


/* Match a character from a list and remove this from a string */
 void removeCharsFromString(std::string &str, const char charsToRemove[]){
  for (unsigned int i = 0; i < strlen(charsToRemove); ++i){
   str.erase(std::remove(str.begin(), str.end(), charsToRemove[i]), str.end());
  }
}

/* Match a string and remove this from an origin string */
std::string purgestring(std::string s, std::string rm)
{
  std::string t = s;
  std::string::size_type i = t.find(rm);
   while(i != std::string::npos) {
     t.erase(i, rm.length());
     i = t.find(rm, i);
   }
  return t;
}

/* Split a string according a delimiter and get back a list of string*/
std::vector<std::string> strsplit(const std::string &s, char delim){
  std::vector<std::string> list;
  std::stringstream ss(s);
  std::string item;
  while(getline(ss, item, delim)){
    if(!item.empty())
      list.push_back(item);
  }
  return list;
}

/* Trim a string */
std::string trim(const std::string &s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && isspace(*it))
        it++;

    std::string::const_reverse_iterator rit = s.rbegin();
    while (rit.base() != it && isspace(*rit))
        rit++;

    return std::string(it, rit.base());
}


/* Check if the string is a number */
bool isNumber(std::string str){
  double d;
  std::istringstream is(str);
  is >> d;
  return !is.fail() && is.eof();
}

/* Convert a string to a double number */
double stod_(std::string str){
  std::replace(str.begin(), str.end(), ',', '.');
  return atof(str.c_str());
}

/* Convert a string to lower without modifing the input */
std::string lower(std::string str){
  std::string res = str;
  std::transform(res.begin(), res.end(), res.begin(), ::tolower);
  return res;
}
