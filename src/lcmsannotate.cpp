/* lcmsannotate.cpp (DynMetId)
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

#include "lcmsannotate.h"
#include "stroperation.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>



#include <mysql.h>

using namespace std;

/* Method to import from a mysql table into a local string matrix */
int LCMSAnnotate::init(std::string dbhost, std::string user, std::string password, std::string dbname, std::string dbtabname)
{
  MYSQL *connect;
  connect = mysql_init(NULL);

  if(!connect){
    std::cout << ">> MySQL initialization failed! <<" << std::endl;
    exit(1);
  }
  /*
  printf("MySQL Connection Info: %s \n", mysql_get_host_info(connect));
  printf("MySQL Client Info: %s \n", mysql_get_client_info());
  printf("MySQL Server Info: %s \n", mysql_get_server_info(connect));
  */
  connect = mysql_real_connect(connect, dbhost.c_str(), user.c_str(), password.c_str(), dbname.c_str(), 0, NULL, 0);

  if(!connect){
    std::cout << ">> MySQL connection failed! <<" << std::endl;
    exit(1);
  }


  MYSQL_RES *res_set;
  MYSQL_ROW row;
  /* Store the column names in memory */
  std::string query = "SHOW COLUMNS FROM "+dbtabname;
  mysql_query(connect, query.c_str());
  res_set = mysql_store_result(connect);
  while(((row = mysql_fetch_row(res_set)) !=NULL)){
    header.push_back(row[0]);
  }
  ncol = header.size();
  mysql_free_result(res_set);

  /* Store the mysql table in a string table */
  query = "SELECT * FROM "+dbtabname;
  mysql_query (connect, query.c_str());
  res_set = mysql_store_result(connect);
  //size_t numrows = mysql_num_rows(res_set);

  while(((row = mysql_fetch_row(res_set)) !=NULL)){
    dbtable.push_back(std::vector<std::string>());
    for(size_t i = 0; i < ncol; i++){
      dbtable.back().push_back(row[i]);
    }
  }
  mysql_free_result(res_set);

  nrow = dbtable.size();
  mysql_close (connect);
  return 0;
}

/*Clean the matrix */
void LCMSAnnotate::clear()
{
  for(size_t i = 0; i < dbtable.size(); i++){
    dbtable[i].clear();
  }
  dbtable.clear();
}

/*Select the column id using a name from header list */
int LCMSAnnotate::getdbid(std::string name){
  //std::cout << name << " " << lower_(name) << " " << name << std::endl;
  //std::transform(name.begin(), name.end(), name.begin(), ::tolower);
  std::vector<std::string>::iterator it = std::find(header.begin(), header.end(), lower(name));
  if (it == header.end()){
    // name not in vector
    return -1;
  }
  else{
    return std::distance(header.begin(), it);
  }
}

std::vector<std::string> LCMSAnnotate::parseqline(std::string qline)
{
  std::vector<std::string> r;
  std::vector<std::string> a = strsplit(qline, ';');
  for(size_t i = 0; i < a.size(); i++){
    std::vector<std::string> b = strsplit(a[i], ' ');
    for(size_t j = 0; j < b.size(); j++){
      removeCharsFromString(b[j], (char*)":%");
      r.push_back(b[j]);
    }
  }
  return r;
}

/* find examples:
 * find(MS: 233.2844 error: 5ppm; tR: 10.25 error: 5%);
 * or
 * find(MS: 233.2844 error: 5ppm;);
 * or
 * find(tR: 10.25 error: 5%);
 * or
 * find(Name: Glucuronic acid);
 */

struct mapres { // map results data structure...
  mapres(double ms_error_, double tr_error_, std::string row_) :
      ms_error(ms_error_),
      tr_error(tr_error_),
      row(row_){}
  double ms_error;
  double tr_error;
  std::string row;

  bool operator<(const mapres& other) const
  {
    if(ms_error < other.ms_error){
      if(tr_error < other.tr_error){
        return true;
      }
      else{
        return false;
      }
    }
    else if(fabs(ms_error - other.ms_error) < 1e-4){
      if(tr_error < other.tr_error){
        return true;
      }
      else{
        return false;
      }
    }
    else{
      return false;
    }
  }
};

std::vector<std::string> LCMSAnnotate::find(std::string qline)
{
  //std::cout << ">>>>> Search: <<<<<\n" << qline << std::endl;
  std::vector<std::string> q = parseqline(qline);
  std::vector<int> found; // here we put the db id!!
  bool refine = false;
  /* DEBUG QUERY STRING
  std::cout << "Parsing out: " << std::endl;
  for(int i = 0; i < (int)q.size(); i++){
    std::cout << q[i] << std::endl;
  }
  std::cout << "---------------" << std::endl;
  */

  // Global variables
  double ms = 0.f;
  double add = 0.f;
  double tr = 0.f;
  double pred_trerr = 0.f;
  double emp_trerr = 0.f;
  double vm = 0.f;
  double vd = 0.f;
  double flow = 0.f;
  double init_B = 0.f;
  double final_B = 0.f;
  double tg = 0.f;

  /* select the headers ID
   * N.B.: these name must be fixed in mysql database!
   */
  int idName = getdbid("name");
  int idMS = getdbid("ms");
  int idFlag = getdbid("flag");
  int idLogKw = getdbid("logkw");
  int idS = getdbid("s");

  //std::cout << "ID PARAMETERS " << idName << " " << idMS << " " << idFlag << " " << idLogKw << " " << idS << std::endl;

  for(size_t i = 0; i < q.size(); i+=2){
    if(q[i].compare("name") == 0){
      for(size_t j = 0; j < dbtable.size(); j++){
        if(dbtable[j][idName].compare(q[i+1]) == 0){
          found.push_back(j);
        }
        else{
          continue;
        }
      }
      refine = true;
    }
    else if(q[i].compare("ms") == 0){
      if(q[i+1].find("m/z") != std::string::npos)
        purgestring(q[i+1], "m/z");
      else if(q[i+1].find("n") != std::string::npos)
        purgestring(q[i+1], "n");

      //std::cout << "MS converted: " << q[i+1] << std::endl;
      ms = stod_(q[i+1]);
      double ppm = 0.f;
      size_t sz = 0;
      if(i+5 <= q.size())
        sz = i+5;
      else
        sz = q.size();

      for(size_t j = 0; j < sz; j++){
        if(q[j].compare("error") == 0){
          ppm = stod_(purgestring(q[j+1], "ppm"));
        }
        else if(q[j].compare("add") == 0){
          add = stod_(q[j+1]);
        }
        else
          continue;
      }

      double mserror = DaltonError(ms, ppm);

      //std::cout << "Search for... " << FloatToString(ms, 4) << " " << ppm << " " << mserror << " "<<  add << std::endl;
      if(found.size() == 0 && refine == false){ // search starting from MS then refine...
        for(size_t j = 0; j < dbtable.size(); j++){
          //std::cout << stod_(dbtable[j][idMS]) << " " << add << " " << stod_(dbtable[j][idMS]) + add << " " << ms << std::endl;
          if(std::fabs((stod_(dbtable[j][idMS])+add) - ms) <= mserror){
            //std::cout << stod_(dbtable[j][idMS]) << " " << add << " " << stod_(dbtable[j][idMS]) + add << " " << ms << " " << mserror << std::endl;
            found.push_back(j);
          }
          else{
            continue;
          }
        }
      }
      else{ // refine search by MS
        size_t j = 0;
        while(j < found.size()){
          if(std::fabs(stod_(dbtable[found[j]][idMS]) - ms) <= mserror){
            continue;
          }
          else{
            found.erase(found.begin()+j);
            j = 0;
          }
        }
      }
      refine = true;
    }
    else if(q[i].compare("tr") == 0){
      tr = stod_(q[i+1]);

      /*if(i+1+15 > q.size()){
        std::cout << "Error in method definition!" << std::endl;
      }*/

      for(size_t j = i+1; j < q.size(); j++){
        if(q[j].compare("emperror") == 0){
          emp_trerr = stod_(purgestring(q[j+1], "%"));
        }
        if(q[j].compare("prederror") == 0){
          pred_trerr = stod_(purgestring(q[j+1], "%"));
        }
        else if(q[j].compare("init") == 0){
          init_B = stod_(q[j+1]);
          if(init_B > 1)
            init_B /= 100.f;
        }
        else if(q[j].compare("final") == 0){
          final_B = stod_(q[j+1]);
          if(final_B > 1)
            final_B /= 100.f;
        }
        else if(q[j].compare("tg") == 0){
          tg = stod_(q[j+1]);
        }
        else if(q[j].compare("flow") == 0){
          //std::cout << q[j] << " " << q[j+1] << std::endl;
          flow = stod_(q[j+1]);
        }
        else if(q[j].compare("vm") == 0){
          vm = stod_(q[j+1]);
        }
        else if(q[j].compare("vd") == 0){
          vd = stod_(q[j+1]);
        }
        else{
          continue;
        }
      }

      //std::cout << "tr: " << tr << " empirical err: " << emp_trerr  << " predicted error: " << pred_trerr << " vm: " << vm << " vd: " << vd << " flow: " << flow << " init: " << init_B << " final: " << final_B << " tg: " << tg << std::endl;
      if(found.size() == 0 && refine == false){ // search starting from tR
        for(size_t j = 0; j < dbtable.size(); j++){
          if(tr > -1 && emp_trerr > -1 && pred_trerr > -1){
            double logkw = stod_(dbtable[j][idLogKw]);
            double s = stod_(dbtable[j][idS]);
            double tr_pred = rtpred(logkw, s, vm, vd, flow, init_B, final_B, tg);
            //Check if is it a predicted or experimental logkw and s
            double perr;
            if(lower(dbtable[j][idFlag]).compare("experimental") == 0){
              perr = emp_trerr;
            }
            else{
              perr = pred_trerr;
            }

            if(std::fabs((tr - tr_pred)/tr)*100.f <= perr){
              found.push_back(j);
            }
            else{
              continue;
            }
          }
          else{
            // no retention error comparisson needed
            continue;
          }
        }
      }
      else{ //refine search
        size_t j = 0;
        while(j < found.size()){
          if(tr > -1 && emp_trerr > -1 && pred_trerr > -1){
            double logkw = stod_(dbtable[found[j]][idLogKw]);
            double s = stod_(dbtable[found[j]][idS]);
            double tr_pred = rtpred(logkw, s, vm, vd, flow, init_B, final_B, tg);
            //Check if is it a predicted or experimental logkw and s
            double perr;
            if(lower(dbtable[found[j]][idFlag]).compare("experimental") == 0){
              perr = emp_trerr;
            }
            else{
              perr = pred_trerr;
            }
            if(std::fabs((tr - tr_pred)/tr)*100.f <= perr){
              j++;
            }
            else{
              found.erase(found.begin()+j);
              j = 0;
            }
          }
          else{
            // no retention error comparisson needed
            //j++;
          }
        }
      }
      refine = true;
    }
    else{
      continue;
    }
  }

  // Give output
  // Name;mass_identified;error_ms;tr;error_tr %;
  // rank output from the small errors in retention time and ppm
  std::vector<std::string> res;
  if(found.size() == 0){
    //std::cout << "No object found!" << std::endl;
    return res;
  }
  else{
    //std::vector<std::string> lines = getqline(found);
    std::vector<mapres> tmprow;
    for(size_t i = 0; i < found.size(); i++){
      //std::vector<std::string> v = strsplit(lines[i], ';');
      //std::string row;
      std::stringstream row;
      double ms_error = 0.f;
      double tr_error = 0.f;
      double mass_calculated = 0.f;
      //std::cout << "ms calculation & rt calculation" << std::endl;
      //std::cout << tr << " " << vm << " " << vd << std::endl;
      if(init_B > 0 && final_B > 0 && tg > 0 && vm > 0 && vd > 0){ // rt calculation
        double logkw = stod_(dbtable[found[i]][idLogKw]);
        double s = stod_(dbtable[found[i]][idS]);
        //std::cout << logkw << " " << s << std::endl;
        double tr_pred = rtpred(logkw, s, vm, vd, flow, init_B, final_B, tg);
        //std::cout << tr_pred << std::endl;
        if(ms > 0){
          // ppm error = (mass_observed - mass_calcuated) / mass_calcualted * 1e6
          mass_calculated = stod_(dbtable[found[i]][idMS])+add;
          //std::cout << ms << " " << v[idMS] << " " << add << " " << mass_calcuated << std::endl;
          ms_error = PPMError(ms, mass_calculated);
        }

        if(tr > 0){
          tr_error = (fabs(tr_pred-tr)/tr)*100;
        }

        row << "name: " << dbtable[found[i]][idName] << ";"; // name
        row << "mass: " << FloatToString(mass_calculated, 4) << ";"; // MS
        row << "mass_error: " << FloatToString(ms_error, 2) << ";"; // MSERROR
        row << "tr: " << FloatToString(tr_pred, 2) << ";"; // Retention time predicted
        row << "tr_error: " << FloatToString(tr_error, 2) << ";"; // Retention time error

        // add other info!
        for(size_t j = 0; j < header.size()-1; j++){
          if(j != (size_t)idName && j != (size_t)idMS && j != (size_t)idLogKw && j != (size_t)idS){
            row << header[j] << ": " << dbtable[found[i]][j] << ";";
          }
          else{
            continue;
          }
        }

        int j = header.size()-1;
        if(j != idName && j != idMS && j != idLogKw && j != idS){
          row << header.back() << ": "+dbtable[found[i]].back();
        }
      }
      else{
        row << "name: " << dbtable[found[i]][idName] << ";"; // name
        row << "mass: " << FloatToString(mass_calculated, 4) << ";"; // MS

        mass_calculated = stod_(dbtable[found[i]][idMS])+add;
        ms_error = PPMError(ms, mass_calculated);

        row << "mass_error: " << FloatToString(ms_error, 3) << ";"; // mass error

        // add other info!
        for(size_t j = 0; j < header.size()-1; j++){
          if(j != (size_t)idName && j != (size_t)idMS && j != (size_t)idLogKw && j != (size_t)idS){
            row << header[j] << ": " << dbtable[found[i]][j] << ";";
          }
          else{
            continue;
          }
        }

        int j = header.size()-1;
        if(j != idName && j != idMS && j != idLogKw && j != idS){
          row << header.back() << ": " << dbtable[found[i]].back();
        }
      }
      tmprow.push_back(mapres(ms_error, tr_error, row.str()));
    }

    /*sort the output*/
    std::sort(tmprow.begin(), tmprow.end());
    for(size_t i = 0; i < tmprow.size(); i++)
      res.push_back(tmprow[i].row);
    return res;
  }
}

void LCMSAnnotate::setRTLinearCorrection(std::string rttunfile, std::string qline)
{
  rtslope = 1.f;
  rtintercept = 0.f;
  std::vector<double> x;
  std::vector<double> y;
  std::string line;
  std::ifstream frttun(rttunfile.c_str());
  if(frttun.is_open()){
    while(getline(frttun, line)){
      std::vector<std::string> v = strsplit(line, ';');
      if(v.size() == 2){
        std::stringstream ss;
        ss << "name " << v[0] << "; tr " << v[1] << " emperror: 100.0 prederror: 100.0 " << qline;

        /* qline is of type: Name: HMDB00253; tR: -1 error: -1 init: 5 final: 95 tg: 14 flow: 0.3 vm: 0.3099 vd: 0.375 */
        std::vector<std::string> res = find(ss.str());


        if(res.size() == 1){
          std::cout << res[0] << std::endl;
          std::vector<std::string> a = strsplit(res[0], ';');
          std::vector<std::string> b = strsplit(a[3], ':');
          x.push_back(stod_(v[1]));
          y.push_back(stod_(b[1]));
        }
        else{
          continue;
        }
      }
    }

    std::vector<std::string> q = parseqline(qline);

    /*Calculate rtslope and rtintercept by linear regression.
     * y = ax + b
     *  b = [Sum y_i(x_i-X)] / [Sum x_i(x_i-X)] were X is the x average;
     *  a = Y-Xb were Y is the y average;
     */

    // Calculating averages
    double x_avg = 0.f; double y_avg = 0.f;
    for(size_t i = 0; i < x.size(); i++){
      x_avg += x[i];
      y_avg += y[i];
    }
    x_avg /= (double)x.size();
    y_avg /= (double)y.size();

    double n = 0.f; double d = 0.f;
    for(size_t i = 0; i < x.size(); i++){
      n += y[i]*(x[i]-x_avg);
      d += x[i]*(x[i]-x_avg);
    }

    // Calculation tR corrective parameters.
    rtslope = n/d;
    rtintercept = y_avg - x_avg*rtslope;
    frttun.close();
  }
}

double LCMSAnnotate::rtpred(double logkw, double s,
  double vm, double vd, double flow, double init_B, double final_B, double tg)
{
  double t0 = vm/flow;
  double td = vd/flow;
  double DeltaFi = final_B - init_B;
  double b = (t0 * DeltaFi * s) / tg;
  double k0 = pow(10, (logkw - s*(init_B)));
  double trpred = ((t0/b) * log10(2.3*k0*b))+ t0 + td;
  /* y = ax + b where in this case rtslope
  * and rtintercept are estimated using already identified
  * compounds in a run
  * N.B.: if no setRTTuning is run, rtslope = 1 and rtintercept = 0
  */
  return (trpred - rtintercept)/rtslope;
}

std::vector<std::string> LCMSAnnotate::db2MetaScope(std::string chromparams)
{
  std::vector<std::string> q = parseqline(chromparams);
  // Global variables
  double vm = 0.f;
  double vd = 0.f;
  double flow = 0.f;
  double init_B = 0.f;
  double final_B = 0.f;
  double tg = 0.f;
  std::vector<std::string> r;

  for(size_t i = 0; i < q.size(); i++){
    if(q[i].compare("init") == 0){
      init_B = stod_(q[i+1]);
      if(init_B > 1)
        init_B /= 100.f;
    }
    else if(q[i].compare("final") == 0){
      final_B = stod_(q[i+1]);
      if(final_B > 1)
        final_B /= 100.f;
    }
    else if(q[i].compare("tg") == 0){
      tg = stod_(q[i+1]);
    }
    else if(q[i].compare("flow") == 0){
      flow = stod_(q[i+1]);
    }
    else if(q[i].compare("vm") == 0){
      vm = stod_(q[i+1]);
    }
    else if(q[i].compare("vd") == 0){
      vd = stod_(q[i+1]);
    }
    else{
      continue;
    }
  }

  int idName = getdbid("name");
  int idMS = getdbid("ms");
  int idFormula = getdbid("formula");
  int idFlag = getdbid("flag");
  int idLogKw = getdbid("logkw");
  int idS = getdbid("s");

  if(init_B > 0 && final_B > 0 && tg > 0 && vm > 0 && vd > 0){ // rt calculation
    r.push_back(format("Name;Compound ID;Neutral Mass;Retention time (min);Formula;Flag"));
    for(size_t i = 0; i < dbtable.size(); i++){
      double logkw = stod_(dbtable[i][idLogKw]);
      double s = stod_(dbtable[i][idS]);
      //std::cout << logkw << " " << s << std::endl;
      double tr_pred = rtpred(logkw, s, vm, vd, flow, init_B, final_B, tg);
      r.push_back(format("%s;%s;%s;%.2f;%s;%s", dbtable[i][idName].c_str(), dbtable[i][idName].c_str(), dbtable[i][idMS].c_str(), tr_pred, dbtable[i][idFormula].c_str(), dbtable[i][idFlag].c_str()));
    }
    return r;
  }
  else{
    return r;
  }
}

// Private methods
double LCMSAnnotate::DaltonError(double mass, double ppm)
{
  // error(Da) = (theoretical mass) / ((1/errorPPM)*1,000,000)
  return mass/((1/ppm) * 1000000);
}

double LCMSAnnotate::PPMError(double mass, double theor_mass)
{
  return ((mass - theor_mass)/theor_mass)*1e6;
}

double LCMSAnnotate::pround(double x, int precision)
{
  /*stringstream stream;
  stream << fixed << setprecision(precision) << x;
  stream.str();*/
  /*if (x == 0.f)
    return x;
  int ex = floor(log10(fabs(x))) - precision + 1;
  double div = pow(10, ex);
  return floor(x / div + 0.5) * div;*/
  int dst;
  double tmp;
  double result;
  tmp = x * pow(10, precision);
  if(tmp < 0) {//negative double
     dst = (int)(tmp - 0.5);
  }else {
     dst = (int)(tmp + 0.5);
  }
  result = (double)((double)dst * pow(10, -precision));
  return result;

}
