#include "particle.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

void Particle::setName(std::string& name) {
  name_ = name;
}

std::string Particle::getName() const {
  return name_;
}

void Particle::setPartId(int pid) {
  pid_ = pid;
}

int Particle::getDegen() const {
  return degen_;
}

void Particle::setDegen(int degen) {
  degen_ = degen;
}

int Particle::getPartId() const {
  return pid_;
}

void Particle::setMass(double mass) {
  mass_ = mass;
}

double Particle::getMass() const {
  return mass_;
}

void Particle::setStrangeness(int strangeness) {
  strangeness_ = strangeness;
}

int Particle::getStrangeness() const {
  return strangeness_;
}

void Particle::setBottomness(int bottomness) {
  bottomness_ = bottomness;
}

int Particle::getBottomness() const {
  return bottomness_;
}

void Particle::setCharmness(int charmness) {
  charmness_ = charmness;
}

int Particle::getCharmness() const {
  return charmness_;
}

void Particle::setBaryonNumber(int baryonNumber) {
  baryonNumber_ = baryonNumber;
}

int Particle::getBaryonNumber() const {
  return baryonNumber_;
}

void Particle::setCharge(int charge) {
  charge_ = charge;
}

int Particle::getCharge() const {
  return charge_;
}

void Particle::setIsospin3(double isospin3) {
  isospin3_ = isospin3;
}

double Particle::getIsospin3() const {
  return isospin3_;
}

void Particle::read_fastreso_components(const std::string& folderpath, const std::vector<double>& temparr, const std::vector<double>& pbararr) {
  
  size_t Ntemp = temparr.size() ;
  size_t Npbar = pbararr.size() ;
  
  temp.resize(Ntemp) ;
  pbar.resize(Npbar) ;
  
  temp = temparr;
  pbar = pbararr;
  
  std::sort(temp.begin(), temp.end());
  
  mintemp = temp[0];
  maxtemp = temp[Ntemp-1];
  minpstar = pbar[0];
  maxpstar = pbar[Npbar - 1];
  
  f1.resize(Ntemp, std::vector<double>(Npbar, 0.0));
  f2.resize(Ntemp, std::vector<double>(Npbar, 0.0));
  
  std::ostringstream filename;
  
  for (size_t i = 0; i < Ntemp; ++i) {    
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(4) << temp[i];
    std::string formattedtemperature = ss.str();
    
    filename.str("");

    filename << folderpath + "PDGid_" + std::to_string(pid_) + "_total_T" + formattedtemperature + "_Fj.out";
   
    std::ifstream filedata(filename.str());
    
    // Skip the first three lines
    for (int line = 0; line < 3; ++line) {
      std::string dummy_line;
      getline(filedata, dummy_line);
    }
      
    double dummymass = 0, dummypbar = 0;
    
    for (size_t line = 0; line < Npbar; ++line) {
      std::string line_data;
      
      getline(filedata, line_data);
      std::istringstream iss(line_data);
      int j = line-3;
      iss >> dummypbar >> dummymass >> f1[i][j] >> f2[i][j];
    }
    filedata.close();
  }
  
  size_t nx = Ntemp; // number of files
  size_t ny = Npbar; // number of lines in file
  
  // Allocate memory for arrays
  double zaf1[Ntemp * Npbar];
  double zaf2[Ntemp * Npbar];
      
  // GSL Interpolation setup
  splinef1 = gsl_spline2d_alloc(T, nx, ny);
  splinef2 = gsl_spline2d_alloc(T, nx, ny);
  xacc = gsl_interp_accel_alloc();
  yacc = gsl_interp_accel_alloc();
  
  // Set GSL Interpolation data
  for (size_t i = 0; i < nx; ++i) {
    for (size_t j = 0; j < ny; ++j) {
      gsl_spline2d_set(splinef1, zaf1, i, j, f1[i][j]);
      gsl_spline2d_set(splinef2, zaf2, i, j, f2[i][j]);
    }
  }
  
  gsl_spline2d_init(splinef1, temp.data(), pbar.data(), zaf1, nx, ny);
  gsl_spline2d_init(splinef2, temp.data(), pbar.data(), zaf2, nx, ny);
}

void Particle::clear(){
  gsl_spline2d_free(splinef1);
  gsl_spline2d_free(splinef2);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
}

double Particle::interpolated_f1(double temp_value, double pbar_value) const {
  if (temp_value < mintemp || temp_value > maxtemp || pbar_value < minpstar || pbar_value > maxpstar) {
    return 0.0;
  }
  
  return gsl_spline2d_eval(splinef1, temp_value, pbar_value, xacc, yacc);
}

double Particle::interpolated_f2(double temp_value, double pbar_value) const {
  if (temp_value < mintemp || temp_value > maxtemp || pbar_value < minpstar || pbar_value > maxpstar) {
    return 0.0;
  }
  
  return gsl_spline2d_eval(splinef2, temp_value, pbar_value, xacc, yacc);
}
