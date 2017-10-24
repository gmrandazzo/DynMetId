/* stroperation.cpp (DynMetId)
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

#include "stroperation.h"
#include <algorithm>
#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <locale>
#include <cstdarg>


/* Match a substring (search) in a string (subject) and replace with another substring (replace)*/
void strreplace(std::string& subject, const std::string& search, const std::string& replace) {
    size_t pos = 0;
    while((pos = subject.find(search, pos)) != std::string::npos) {
         subject.replace(pos, search.length(), replace);
         pos += replace.length();
    }
}


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

/* Format various char * arguments and convert to a string
 * N.B.: Safe and convenient but not exactly efficient implementation
 */
std::string format(const char* fmt, ...){
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

/* Convert int number to string */
std::string intost(int i)
{
  std::stringstream ss;
  ss << i;
  return ss.str();
}

/* Convert a string to lower without modifing the input */
std::string lower(std::string str){
  std::string res = str;
  std::transform(res.begin(), res.end(), res.begin(), ::tolower);
  return res;
}
