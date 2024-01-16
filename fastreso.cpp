#include "fastreso.h"
#include <fstream>
#include <sstream>
#include <iostream>

FastReso::FastReso(const Particle& particle, const Freezeout& freezeout)
  : particle_(particle), freezeout_(freezeout) {
  
  pTarr.resize(NpT);
  for (size_t i = 0; i < NpT; ++i) {
    pTarr[i] = pT_min + (pT_max - pT_min) * i / (NpT - 1);
  }
  
  phiarr.resize(Nphi);
  for (size_t j = 0; j < Nphi; ++j) {
    phiarr[j] = phi_min + (phi_max - phi_min) * j / (Nphi - 1);
  }
  
  yarr.resize(Ny);
  for (size_t k = 0; k < Ny; ++k) {
    yarr[k] = y_min + (y_max - y_min) * k / (Ny - 1);
  }
  
  EdNd3p.resize(NpT, std::vector<std::vector<double>>(Nphi, std::vector<double>(Ny, 0.0)));
  dNpTdpTdy.resize(NpT, std::vector<double>(Ny, 0.0));
  v1.resize(NpT, std::vector<double>(Ny, 0.0));
  v2.resize(NpT, std::vector<double>(Ny, 0.0));
  v3.resize(NpT, std::vector<double>(Ny, 0.0));
  v4.resize(NpT, std::vector<double>(Ny, 0.0));
  dNdy.resize(Ny);
}

void FastReso::calc_EdNd3p() {
  
  double gdsig, pu, eqdist, mu, temp, pbar;
  double pmu[4], umu[4], dSigma[4], pf1[4], pf2[4], gmu[4];
  
  for (size_t pT_idx = 0; pT_idx < pTarr.size(); ++pT_idx) {
    double pT = pTarr[pT_idx]; 
    for (size_t phi_idx = 0; phi_idx < phiarr.size(); ++phi_idx) {
      double phi = phiarr[phi_idx];
      for (size_t y_idx = 0; y_idx < yarr.size(); ++y_idx) {
	double y = yarr[y_idx];
	double EdNd3p_ = 0.0 ;
	for (size_t i = 0; i < freezeout_.mysurface.size(); ++i) {
	  const surface_element& surf = freezeout_.mysurface[i];
	  
	  double mT = sqrt(particle_.getMass()*particle_.getMass()+pT*pT);
	  mu = particle_.getBaryonNumber()*surf.muB + particle_.getStrangeness()*surf.muS + particle_.getCharge()*surf.muQ;
	  temp = surf.temperature;
	  
	  pmu[0] = mT*cosh(y-surf.xmu.eta);
	  pmu[1] = pT*cos(phi);
	  pmu[2] = pT*sin(phi);
	  pmu[3] = mT*sinh(y-surf.xmu.eta)/surf.xmu.tau;
	  
	  umu[0] = surf.umu.tau;
	  umu[1] = surf.umu.x;
	  umu[2] = surf.umu.y;
	  umu[3] = surf.umu.eta;
          
	  pu = pmu[0]*umu[0] - pmu[1]*umu[1] - pmu[2]*umu[2] - pow(surf.xmu.tau,2)*pmu[3]*umu[3];
	  pbar = sqrt(pu*pu - particle_.getMass()*particle_.getMass());

	  double f1_ = particle_.interpolated_f1(temp, pbar);
	  double f2_ = particle_.interpolated_f2(temp, pbar);

	  f1_ = f1_/pbar;
          f2_ = f2_/pbar;
		
	  pf1[0] = pmu[0] - pu*umu[0];
	  pf1[1] = pmu[1] - pu*umu[1];
	  pf1[2] = pmu[2] - pu*umu[2];
	  pf1[3] = pmu[3] - pu*umu[3];
          
	  pf2[0] = pu*umu[0];
	  pf2[1] = pu*umu[1];
	  pf2[2] = pu*umu[2];
	  pf2[3] = pu*umu[3];       
	            
	  gmu[0] = f1_*pf1[0] + f2_*pf2[0] ;
	  gmu[1] = f1_*pf1[1] + f2_*pf2[1] ;
	  gmu[2] = f1_*pf1[2] + f2_*pf2[2] ;
	  gmu[3] = f1_*pf1[3] + f2_*pf2[3] ;
	  
	  dSigma[0] = surf.dSigma.tau;
	  dSigma[1] = surf.dSigma.x;
	  dSigma[2] = surf.dSigma.y;
	  dSigma[3] = surf.dSigma.eta;
          
	  gdsig = gmu[0]*dSigma[0] + gmu[1]*dSigma[1] + gmu[2]*dSigma[2] + gmu[3]*dSigma[3];  // GeV fm^3
          
	  //degen = (particle_.spin + 1)*(particle_.isospin3 + 1);
          
	  // Whether to multiply degen or not?
	  
	  EdNd3p_ = EdNd3p_ + gdsig;
	}   
        
	EdNd3p[pT_idx][phi_idx][y_idx] = EdNd3p_;     
      }
    }
  }
}

void FastReso::calc_dNpTdpTdy() {
  // Integrate over phi
  for (size_t pT_idx = 0; pT_idx < pTarr.size(); ++pT_idx) {
    for (size_t y_idx = 0; y_idx < yarr.size(); ++y_idx) {
      double result = 0.0;
      for (size_t phi_idx = 0; phi_idx < phiarr.size() - 1; ++phi_idx) {
	double delta_phi = phiarr[phi_idx + 1] - phiarr[phi_idx];
	result += (EdNd3p[pT_idx][phi_idx][y_idx] + EdNd3p[pT_idx][phi_idx + 1][y_idx]) * delta_phi / 2.0;
      }
      dNpTdpTdy[pT_idx][y_idx] = result;
    }
  }
}

void FastReso::calc_dNdy() {
  // Integrate over phi and pT
  for (size_t y_idx = 0; y_idx < yarr.size(); ++y_idx) {
    double result = 0.0;
    for (size_t pT_idx = 0; pT_idx < pTarr.size() - 1; ++pT_idx) {
      double mid_pT = (pTarr[pT_idx + 1] + pTarr[pT_idx]) / 2.0;
      double delta_pT = pTarr[pT_idx + 1] - pTarr[pT_idx];
      result += mid_pT * (dNpTdpTdy[pT_idx][y_idx] + dNpTdpTdy[pT_idx + 1][y_idx]) * delta_pT / 2.0;
    }
    dNdy[y_idx] = result;
  }
}

void FastReso::calc_vn() {
  
  for (size_t pT_idx = 0; pT_idx < pTarr.size(); ++pT_idx) {
    for (size_t y_idx = 0; y_idx < yarr.size(); ++y_idx) {
      double num1 = 0.0, num2 = 0.0, num3 = 0.0, num4 = 0.0;
      double den = 0.0;
      for (size_t phi_idx = 0; phi_idx < phiarr.size() - 1; ++phi_idx) {
	double delta_phi = phiarr[phi_idx + 1] - phiarr[phi_idx];
	double mid_phi = (phiarr[phi_idx + 1] + phiarr[phi_idx]) / 2.0;
        
	num1 += (std::cos(phiarr[phi_idx]) * EdNd3p[pT_idx][phi_idx][y_idx]
                 + std::cos(phiarr[phi_idx + 1]) * EdNd3p[pT_idx][phi_idx + 1][y_idx]) * delta_phi/ 2.0;
	num2 += (std::cos(2 * phiarr[phi_idx]) * EdNd3p[pT_idx][phi_idx][y_idx] 
                 + std::cos(2 * phiarr[phi_idx + 1]) * EdNd3p[pT_idx][phi_idx + 1][y_idx]) * delta_phi/ 2.0;
	num3 += (std::cos(3 * phiarr[phi_idx]) * EdNd3p[pT_idx][phi_idx][y_idx]
                 + std::cos(3 * phiarr[phi_idx + 1]) * EdNd3p[pT_idx][phi_idx + 1][y_idx]) * delta_phi/ 2.0;
	num4 += (std::cos(4 * phiarr[phi_idx]) * EdNd3p[pT_idx][phi_idx][y_idx] 
                 + std::cos(4 * phiarr[phi_idx + 1]) * EdNd3p[pT_idx][phi_idx + 1][y_idx]) * delta_phi/ 2.0;
	den += (EdNd3p[pT_idx][phi_idx][y_idx] + EdNd3p[pT_idx][phi_idx + 1][y_idx]) * delta_phi/ 2.0;
      }
      if (den != 0.0) {
	v1[pT_idx][y_idx] = num1/den;
	v2[pT_idx][y_idx] = num2/den;
	v3[pT_idx][y_idx] = num3/den;
	v4[pT_idx][y_idx] = num4/den;
      }
    }
  }
}

void FastReso::calc_observables() {
  calc_EdNd3p() ;
  calc_dNpTdpTdy();
  calc_dNdy();
  calc_vn();
}

void FastReso::output() {
  std::ostringstream filenameStream;
  filenameStream << particle_.getName() << "_dNpTdpTdy_vn.dat";
  std::string filename = filenameStream.str();
  
  std::ofstream outputfile(filename);
  
  if (!outputfile.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }
  
  outputfile << "#pT\t\t y\t\t dNpTdpTdy\t\t v1\t\t v2\t\t v3\t\t v4\n";
  
  for (size_t pT_idx = 0; pT_idx < pTarr.size(); ++pT_idx) {
    for (size_t y_idx = 0; y_idx < yarr.size(); ++y_idx) {
      outputfile << pTarr[pT_idx] << "\t" << yarr[y_idx] << "\t" 
                 << dNpTdpTdy[pT_idx][y_idx] << "\t" << v1[pT_idx][y_idx] 
                 << "\t" << v2[pT_idx][y_idx] << "\t" << v3[pT_idx][y_idx] 
                 << "\t" << v4[pT_idx][y_idx] << "\n";
    }
  }
  
  outputfile.close();
  
  filenameStream.str(""); // Clear the stream
  filenameStream << particle_.getName() << "_dNdy.dat";
  filename = filenameStream.str();
  
  std::ofstream outputfile2(filename);
  
  if (!outputfile2.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }
  
  outputfile2 << "#y\t\t dNdy\n";
  
  for (size_t y_idx = 0; y_idx < yarr.size(); ++y_idx) {
    outputfile2 << yarr[y_idx] << "\t" << dNdy[y_idx] << "\n";
  }
  
  outputfile2.close();
}
