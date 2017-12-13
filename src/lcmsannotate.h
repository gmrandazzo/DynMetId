/* lcmsannotate.h (DynMetId)
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

#ifndef LCMSANNOTATE_H
#define LCMSANNOTATE_H

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <limits>
#include <cmath>

/*
 * Dynamic Metabolite LCMSAnnotate
 *
 * This is software will annotate for
 * - m/z with molecular adducts
 * - retention time
 * TODO:
 * - isotopic pattern LCMSAnnotate
 * - ms/ms spectra LCMSAnnotate
 */

class LCMSAnnotate
{
public:
  LCMSAnnotate(){ nrow = ncol = 0; rtslope = 1.; rtintercept = 0.; }
  ~LCMSAnnotate(){ for(size_t i = 0; i < dbtable.size(); i++) dbtable[i].clear(); dbtable.clear(); };
  /* init the database structure. Charge in memory the mysql table */
  int init(std::string dbhost, std::string user, std::string password, std::string dbname, std::string dbtabname);

  /* find examples:
   * find(MS 233.2844 within error 5ppm at tR 10.25 within error 5%);
   * or
   * find(MS 233.2844 within error 5ppm);
   * or
   * find(tR 10.25 within error 5%);
   * or
   * find(Name Glucuronic acid);
   */
  std::vector<std::string> find(std::string qline);

  /* Export the whole database into a metascope file */
  std::vector<std::string> db2MetaScope(std::string chromparams);

  //std::vector<std::string> dbChromatogram(std::string chromparams);

  /*In case of retention time shift realign compounds according a simple linear regression. */
  void setRTLinearCorrection(std::string rttunfile, std::string qline);
  void setRTLinearCorrection(double rtslope_, double rtintercept_) { rtslope = rtslope_; rtintercept = rtintercept_; }

  /* Check if the database is empty */
  bool isEmpty(){ if(dbtable.size() == 0) return true; else return false; }

  /* Clean the database */
  void clear();

  /* Search in the header the id of name in order to make the search */
  int getdbid(std::string name);
private:
  // private methods
  template< typename T > typename std::vector<T>::iterator
    insert_sorted(std::vector<T> & vec, T const& item){
      return vec.insert(std::upper_bound(vec.begin(), vec.end(), item), item);
  }

  //std::vector<int> find_all_keys_id(std::vector<key_value> collection, std::string key);

  std::vector<std::string> parseqline(std::string qline);
  /* Calculate the retention time accordin LSS parameters and chromatographic conditions*/
  double rtpred(double logkw, double s, double vm, double vd, double flux,
    double init_B, double final_B, double tg);
  /*Calculate the mass error in dalton*/
  double DaltonError(double mass, double ppm);
  double PPMError(double mass, double theor_mass);
  void NameSearch(int idName, std::string name, std::vector<int> *found);
  void MSSearch(int idMS, double ms, double add, double mult, double mserror, bool is_neutral, std::vector<int> *found);
  void RTSearch(int idLogKw, int idS, int idFlag,
                double tr, double emp_trerr, double pred_trerr, double vm, double vd,
                double flow, double init_B, double final_B, double tg,
                std::vector<int> *found);
  /*reduce the floating point precision to a defined "precision" */
  double pround(double x, int precision);
  /* private data*/
  std::vector<std::vector<std::string>> dbtable; // database stored as string matrix
  std::vector<std::string> header;
  double rtslope, rtintercept;
  size_t nrow, ncol;
};

#endif
