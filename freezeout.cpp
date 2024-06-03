#include "freezeout.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

Freezeout::Freezeout(const std::string& filename) : filename_(filename){
  //readFile(filename_); 
  readFileCart(filename_);
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

//  double c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16;
//  
//  while (surfdata >> fo_surf.xmu.tau >> fo_surf.xmu.x >> fo_surf.xmu.y >> fo_surf.xmu.eta) {
//    surfdata >> fo_surf.dSigma.tau >> fo_surf.dSigma.x >> fo_surf.dSigma.y >> fo_surf.dSigma.eta;
//    surfdata >> fo_surf.umu.tau >> fo_surf.umu.x >> fo_surf.umu.y >> fo_surf.umu.eta;
//    surfdata >> fo_surf.temperature;
//    
//    surfdata >> fo_surf.muB >> fo_surf.muQ >> fo_surf.muS >> c1 >> c2;
//    
//    surfdata >> c3;
//    
//    surfdata >> c4 >> c5 >> c6 >> c7 >> c8 >> c9 >> c10 >> c11 >> c12 >> c13 >> c14 >> c15 >> c16;
//    
//    mysurface.push_back(fo_surf);
//  }
  
  surfdata.close();
  
  NumberofCells = mysurface.size();
  //std::cout << "Number of freezeout cells = " << NumberofCells << std::endl;
}

void Freezeout::readFileCart(const std::string& datafile) {
  
  std::ifstream surfdata(datafile); //the actual data
  
  if (!surfdata.is_open()) {
    std::cerr << "Error opening file: " << datafile << std::endl;
    exit(1);
  }
  
  surface_element fo_surf;
  double dummy,vx,vy,vz;
  mysurface.clear();
  
  while (surfdata >> fo_surf.xmu.tau >> fo_surf.xmu.x >> fo_surf.xmu.y >> fo_surf.xmu.eta) {
    surfdata >> fo_surf.dSigma.tau >> fo_surf.dSigma.x >> fo_surf.dSigma.y >> fo_surf.dSigma.eta;
    surfdata >> dummy >> dummy >> vx >> vy >> vz ;
    
    fo_surf.umu.tau = 1.0/sqrt(1.0-vx*vx-vy*vy-vz*vz);
    fo_surf.umu.x = vx/sqrt(1.0-vx*vx-vy*vy-vz*vz);
    fo_surf.umu.y = vy/sqrt(1.0-vx*vx-vy*vy-vz*vz);
    fo_surf.umu.eta = vz/sqrt(1.0-vx*vx-vy*vy-vz*vz);
    
    surfdata >> dummy >> fo_surf.muB >> fo_surf.temperature >> dummy;

//    surfdata >> fo_surf.shear.element[0] >> fo_surf.shear.element[1] >> fo_surf.shear.element[2]
//	     >> fo_surf.shear.element[3] >> fo_surf.shear.element[4] >> fo_surf.shear.element[5]
//	     >> fo_surf.shear.element[6] >> fo_surf.shear.element[7] >> fo_surf.shear.element[8] 
//	     >> fo_surf.shear.element[9];
//	         
//    surfdata >> fo_surf.bulk;
    
    surfdata >> fo_surf.dbeta[0][0] >> fo_surf.dbeta[0][1] >> fo_surf.dbeta[0][2] >> fo_surf.dbeta[0][3]
	     >> fo_surf.dbeta[1][0] >> fo_surf.dbeta[1][1] >> fo_surf.dbeta[1][2] >> fo_surf.dbeta[1][3]
	     >> fo_surf.dbeta[2][0] >> fo_surf.dbeta[2][1] >> fo_surf.dbeta[2][2] >> fo_surf.dbeta[2][3]
	     >> fo_surf.dbeta[3][0] >> fo_surf.dbeta[3][1] >> fo_surf.dbeta[3][2] >> fo_surf.dbeta[3][3];
    
    mysurface.push_back(fo_surf);
  }

//  double c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16;
//  
//  while (surfdata >> fo_surf.xmu.tau >> fo_surf.xmu.x >> fo_surf.xmu.y >> fo_surf.xmu.eta) {
//    surfdata >> fo_surf.dSigma.tau >> fo_surf.dSigma.x >> fo_surf.dSigma.y >> fo_surf.dSigma.eta;
//    surfdata >> fo_surf.umu.tau >> fo_surf.umu.x >> fo_surf.umu.y >> fo_surf.umu.eta;
//    surfdata >> fo_surf.temperature;
//    
//    surfdata >> fo_surf.muB >> fo_surf.muQ >> fo_surf.muS >> c1 >> c2;
//    
//    surfdata >> c3;
//    
//    surfdata >> c4 >> c5 >> c6 >> c7 >> c8 >> c9 >> c10 >> c11 >> c12 >> c13 >> c14 >> c15 >> c16;
//    
//    mysurface.push_back(fo_surf);
//  }
  
  surfdata.close();
  
  NumberofCells = mysurface.size();
  //std::cout << "Number of freezeout cells = " << NumberofCells << std::endl;
}
