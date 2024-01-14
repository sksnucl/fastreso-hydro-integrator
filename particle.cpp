#include "particle.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

Particle::Particle(const std::string& name, double mass, int strangeness, int charmness, int baryonNumber, int charge, int isospin3, int spin)
  : name_(name), mass_(mass), strangeness_(strangeness), charmness_(charmness), baryonNumber_(baryonNumber), charge_(charge), isospin3_(isospin3), spin_(spin){
  
  degen_ = (spin_ + 1)*(isospin3_ + 1);
  temp = nullptr;
  pbar = nullptr;
  
}

Particle::~Particle() {
  delete temp;
  delete pbar;
  delete f1;
  delete f2;
  gsl_spline2d_free(splinef1);
  gsl_spline2d_free(splinef2);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
  delete[] zaf1;
  delete[] zaf2;
}

std::string Particle::getName() const {
  return name_;
}

double Particle::getMass() const {
  return mass_;
}

int Particle::getStrangeness() const {
  return strangeness_;
}

int Particle::getCharmness() const {
  return charmness_;
}

int Particle::getBaryonNumber() const {
  return baryonNumber_;
}

int Particle::getCharge() const {
  return charge_;
}

int Particle::getSpin() const {
  return spin_;
}

int Particle::getIsospin3() const {
  return isospin3_;
}

void Particle::readf1f2(size_t Ntemp,size_t Npbar) {
  
  temp = new std::vector<double>(Ntemp);
  pbar = new std::vector<double>(Npbar);
  f1 = new std::vector<std::vector<double> >(Ntemp, std::vector<double>(Npbar, 0.0));
  f2 = new std::vector<std::vector<double> >(Ntemp, std::vector<double>(Npbar, 0.0));
  
  //std::vector<double>(Ntemp) templist;
  std::ostringstream filename[Ntemp];
  
  for (size_t i = 0; i < Ntemp; ++i) {
    //templist[i] = mintemp + (maxtemp - mintemp) * static_cast<double>(i) / (Ntemp - 1);
    //(*temp)[i] = templist[i];
    
    (*temp)[i] = mintemp + (maxtemp - mintemp) * static_cast<double>(i) / (Ntemp - 1);
    
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(4) << (*temp)[i];
    std::string formattedtemperature = ss.str();
    
    filename[i] << "../fi_PDG2016+/" + name_ + "_total_T" + formattedtemperature + "_Fj.out";
    
    std::ifstream filedata(filename[i].str());
    
    // Skip the first three lines
    for (int line = 0; line < 3; ++line) {
      std::string dummy_line;
      getline(filedata, dummy_line);
    }
      
    double dummymass = 0, dummypbar = 0;
    
    for (size_t j = 0; j < Npbar; ++j) {
      if (i == 0) {
	filedata >> (*pbar)[j] >> dummymass >> (*f1)[i][j] >> (*f2)[i][j];
      } else {
	filedata >> dummypbar >> dummymass >> (*f1)[i][j] >> (*f2)[i][j];
      }
    }
    
    filedata.close();
  }
  
  nx = Ntemp; // number of files
  ny = Npbar; // number of lines in file
  
  // Allocate memory for arrays
  zaf1 = new double[nx * ny];
  zaf2 = new double[nx * ny];
  
  // GSL Interpolation setup
  splinef1 = gsl_spline2d_alloc(T, nx, ny);
  splinef2 = gsl_spline2d_alloc(T, nx, ny);
  xacc = gsl_interp_accel_alloc();
  yacc = gsl_interp_accel_alloc();
  
  // Set GSL Interpolation data
  for (size_t i = 0; i < nx; ++i) {
    for (size_t j = 0; j < ny; ++j) {
      gsl_spline2d_set(splinef1, zaf1, i, j, (*f1)[i][j]);
      gsl_spline2d_set(splinef2, zaf2, i, j, (*f2)[i][j]);
    }
  }
  
  gsl_spline2d_init(splinef1, temp->data(), pbar->data(), zaf1, nx, ny);
  gsl_spline2d_init(splinef2, temp->data(), pbar->data(), zaf2, nx, ny);
  
  maxpstar = (*pbar)[ny - 1];
  minpstar = (*pbar)[0];
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
