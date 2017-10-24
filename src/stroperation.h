/* stroperation.h (DynMetId)
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

#ifndef STROPERATION_H
#define STROPERATION_H

#include <algorithm>
#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
/*
#define lower(x) std::transform (x.begin(), x.end(), x.begin(), ::tolower)
#define upper(x) std::transform (x.begin(), x.end(), x.begin(), ::toupper)
#define toucfirst(x) std::transform (x.begin(), x.begin()+1, x.begin(),  ::toupper); std::transform (x.begin()+1, x.end(),   x.begin()+1,::tolower)
*/

/* Match a substring (search) in a string (subject) and replace with another substring (replace)*/
void strreplace(std::string& subject, const std::string& search, const std::string& replace);

/* Match a character from a list and remove this from a string */
void removeCharsFromString(std::string &str, const char charsToRemove[]);

/* Match a string and remove this from an origin string */
std::string purgestring(std::string s, std::string rm);

/* Convert a floating point number (double or float) into a string according
 * a certain precision
 */
template <typename T> std::string FloatToString(T Number, int precision){
  std::ostringstream ss;
  ss << std::fixed << std::setprecision(precision) << Number;
  //ss << Number;
  return ss.str();
}

/* Convert a number (int, double or float) into a string */
template <typename T> T StringToNumber(const std::string &Text){
  std::replace(Text.begin(), Text.end(), ',', '.');
  std::istringstream ss(Text);
  T result;
  return ss >> result ? result : 0;
}

/* Split a string according a delimiter and get back a list of string*/
std::vector<std::string> strsplit(const std::string &s, char delim);

/* Format various char * arguments and get back a string */
std::string format(const char* fmt, ...);

/* Trim a string */
std::string trim(const std::string &s);

/* Check if the string is a number */
bool isNumber(std::string str);

/* Convert a string to a double number */
double stod_(std::string str);

/* Convert int number to string */
std::string intost(int i);

/* Convert a string to lower without modifing the input */
std::string lower(std::string str);
#endif
