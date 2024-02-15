#include "particle.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>

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

void Particle::read_fastreso_components(const size_t& decayflag, const std::string& folderpath, const std::vector<double>& temparr, const std::vector<double>& pbararr) {
  
  size_t Ntemp = temparr.size() ;
  size_t Npbar = pbararr.size() ;
  
  f1.resize(Ntemp, std::vector<double>(Npbar, 0.0));
  f2.resize(Ntemp, std::vector<double>(Npbar, 0.0));
  
  std::ostringstream filename;
  
  for (size_t i = 0; i < Ntemp; ++i) {    
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(4) << temparr[i];
    std::string formattedtemperature = ss.str();
    
    filename.str("");

    if(decayflag == 0){
      filename << folderpath + "PDGid_" + std::to_string(pid_) + "_thermal_T" + formattedtemperature + "_Fj.out";
    }else{
      filename << folderpath + "PDGid_" + std::to_string(pid_) + "_total_T" + formattedtemperature + "_Fj.out";
    }
   
    std::ifstream filedata(filename.str());
    
    // Skip the first three lines of datafiles
    for (int line = 0; line < 3; ++line) {
      std::string dummy_line;
      getline(filedata, dummy_line);
    }
      
    double dummymass = 0, dummypbar = 0;
    
    for (size_t line = 0; line < Npbar; ++line) {
      std::string line_data;
      
      getline(filedata, line_data);
      std::istringstream iss(line_data);
      iss >> dummypbar >> dummymass >> f1[i][line] >> f2[i][line];
    }
    filedata.close();
  }
  
}
