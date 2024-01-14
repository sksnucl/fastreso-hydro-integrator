#include "freezeout.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

Freezeout::Freezeout(const std::string& filename) : filename_(filename){
  readFile(filename_);
}

void Freezeout::readFile(const std::string& datafile) {
  
  std::ifstream surfdata(datafile); //the actual data
  
  if (!surfdata.is_open()) {
    std::cerr << "Error opening file: " << datafile << std::endl;
    exit(1);
  }
  
  surface_element fo_surf;
  
  mysurface.clear();
  
  while (surfdata >> fo_surf.xmu.tau >> fo_surf.xmu.x >> fo_surf.xmu.y >> fo_surf.xmu.eta) {
    surfdata >> fo_surf.dSigma.tau >> fo_surf.dSigma.x >> fo_surf.dSigma.y >> fo_surf.dSigma.eta;
    surfdata >> fo_surf.umu.tau >> fo_surf.umu.x >> fo_surf.umu.y >> fo_surf.umu.eta;
    surfdata >> fo_surf.temperature;
    
    surfdata >> fo_surf.muB >> fo_surf.muS >> fo_surf.muQ >> fo_surf.energy >> fo_surf.entropy;
    
    surfdata >> fo_surf.bulk;
    
    surfdata >> fo_surf.shear.element[0] >> fo_surf.shear.element[1] >> fo_surf.shear.element[2]
	     >> fo_surf.shear.element[3] >> fo_surf.shear.element[4] >> fo_surf.shear.element[5]
	     >> fo_surf.shear.element[6] >> fo_surf.shear.element[7] >> fo_surf.shear.element[8] 
	     >> fo_surf.shear.element[9];
    
    mysurface.push_back(fo_surf);
  }
  
  surfdata.close();
  
  NumberofCells = mysurface.size();
  //std::cout << "Number of freezeout cells = " << NumberofCells << std::endl;
}

double Freezeout::EdNd3p(const Particle& particle, double pT, double y, double phi, int eps) {
  
  double sum = 0.0;
  double pdsig, mT, pu, eqdist, mu, temp;
  double pmu[4], umu[4], dSigma[4];
  surface_element surf;
  
  mT = sqrt(particle.getMass()*particle.getMass()+pT*pT);
  
  for (int i = 0; i < NumberofCells; ++i) {
    
    surf = mysurface[i];
    
    mu = particle.getBaryonNumber()*surf.muB + particle.getStrangeness()*surf.muS + particle.getCharge()*surf.muQ;
    temp = surf.temperature;
    
    pmu[0] = mT*cosh(y-surf.xmu.eta);
    pmu[1] = pT*cos(phi);
    pmu[2] = pT*sin(phi);
    pmu[3] = mT*sinh(y-surf.xmu.eta)/surf.xmu.tau;
    
    dSigma[0] = surf.dSigma.tau;
    dSigma[1] = surf.dSigma.x;
    dSigma[2] = surf.dSigma.y;
    dSigma[3] = surf.dSigma.eta;
    
    umu[0] = surf.umu.tau;
    umu[1] = surf.umu.x;
    umu[2] = surf.umu.y;
    umu[3] = surf.umu.eta;
    
    pdsig = pmu[0]*dSigma[0] + pmu[1]*dSigma[1] + pmu[2]*dSigma[2] + pmu[3]*dSigma[3];  // GeV fm^3
    pu = pmu[0]*umu[0] - pmu[1]*umu[1] - pmu[2]*umu[2] - pow(surf.xmu.tau,2)*pmu[3]*umu[3];
    
    eqdist = 1.0/(exp((pu - mu)/temp) + eps);
    
    sum = sum + pdsig*eqdist ;
  }
  
  return sum;
}
