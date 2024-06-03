#include "fastreso.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <omp.h>

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
  
  deltapT = (pT_max - pT_min) / (NpT - 1);
  deltaphi = (phi_max - phi_min) / (Nphi - 1);
  deltay = (y_max - y_min) / (Ny - 1);
  
  EdNd3p.resize(NpT, std::vector<std::vector<double>>(Nphi, std::vector<double>(Ny, 0.0)));
  dNpTdpTdy.resize(NpT, std::vector<double>(Ny, 0.0));
  v1.resize(NpT, std::vector<double>(Ny, 0.0));
  v2.resize(NpT, std::vector<double>(Ny, 0.0));
  v3.resize(NpT, std::vector<double>(Ny, 0.0));
  v4.resize(NpT, std::vector<double>(Ny, 0.0));
  dNdy.resize(Ny);
}

//**************************************************************************************************************

double FastReso::calc_EdNd3p(const double& pT, const double& y, const double& phi, gsl_spline2d *splinef1, gsl_spline2d *splinef2, gsl_interp_accel *xacc, gsl_interp_accel *yacc, const double& mint, const double& maxt, const double& minp, const double& maxp) {
  
  double EdNd3p_ = 0.0 ;
  #pragma omp parallel for reduction(+:EdNd3p_)
  
  for (size_t i = 0; i < freezeout_.mysurface.size(); ++i) {
    const surface_element& surf = freezeout_.mysurface[i];
    double particlemass =  particle_.getMass();
    double mT = sqrt(particlemass*particlemass + pT*pT);
    double temp = surf.temperature;
    
    double pmu[4], umu[4], pu, pbar, pf1[4], pf2[4], gmu[4], dSigma[4], gdsig;
    
    pmu[0] = mT*std::cosh(y-surf.xmu.eta);
    pmu[1] = pT*std::cos(phi);
    pmu[2] = pT*std::sin(phi);
    pmu[3] = mT*std::sinh(y-surf.xmu.eta)/surf.xmu.tau;
    
    umu[0] = surf.umu.tau;
    umu[1] = surf.umu.x;
    umu[2] = surf.umu.y;
    umu[3] = surf.umu.eta;
    
    pu = pmu[0]*umu[0] - pmu[1]*umu[1] - pmu[2]*umu[2] - pow(surf.xmu.tau,2)*pmu[3]*umu[3];
    
    if ((pu*pu - particlemass*particlemass) < 0){
      std::cout << "Error: pbar is NaN" << std::endl;  
      exit(0);
    }
    
    pbar = std::sqrt(pu*pu - particlemass*particlemass);
    
    double f1_, f2_;
    
    if (temp > mint && temp < maxt && pbar > minp && pbar < maxp) {
      f1_ = gsl_spline2d_eval(splinef1, temp, pbar, xacc, yacc);
      f2_ = gsl_spline2d_eval(splinef2, temp, pbar, xacc, yacc);
    }else{
      f1_ = 0.0;
      f2_ = 0.0;
    }
    
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
    
    EdNd3p_ = EdNd3p_ + particle_.getDegen()*gdsig/(pow(2*M_PI,3)*pow(0.197,3)); //(pow(2*M_PI,3)); 
  }   
  
  return EdNd3p_;     // GeV^[-2]
}

//**************************************************************************************************************

double FastReso::calc_EdNd3p_cart(const double& pT, const double& y, const double& phi, gsl_spline2d *splinef1, gsl_spline2d *splinef2, gsl_interp_accel *xacc, gsl_interp_accel *yacc, const double& mint, const double& maxt, const double& minp, const double& maxp) {
  
  double EdNd3p_ = 0.0 ;
  #pragma omp parallel for reduction(+:EdNd3p_)
  
  for (size_t i = 0; i < freezeout_.mysurface.size(); ++i) {
    const surface_element& surf = freezeout_.mysurface[i];
    double particlemass =  particle_.getMass();
    double mT = sqrt(particlemass*particlemass + pT*pT);
    double temp = surf.temperature;
    
    double pmu[4], umu[4], pu, pbar, pf1[4], pf2[4], gmu[4], dSigma[4], gdsig;
    
    pmu[0] = mT*std::cosh(y);
    pmu[1] = pT*std::cos(phi);
    pmu[2] = pT*std::sin(phi);
    pmu[3] = mT*std::sinh(y);
    
    // umu is cartesian
    umu[0] = surf.umu.tau;
    umu[1] = surf.umu.x;
    umu[2] = surf.umu.y;
    umu[3] = surf.umu.eta;
    
    pu = pmu[0]*umu[0] - pmu[1]*umu[1] - pmu[2]*umu[2] - pmu[3]*umu[3];
    
    if ((pu*pu - particlemass*particlemass) < 0){
      std::cout << "Error: pbar is NaN" << std::endl;  
      exit(0);
    }
    
    pbar = std::sqrt(pu*pu - particlemass*particlemass);
    
    double f1_, f2_;
    
    if (temp > mint && temp < maxt && pbar > minp && pbar < maxp) {
      f1_ = gsl_spline2d_eval(splinef1, temp, pbar, xacc, yacc);
      f2_ = gsl_spline2d_eval(splinef2, temp, pbar, xacc, yacc);
    }else{
      f1_ = 0.0;
      f2_ = 0.0;
    }

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
    
    //dSigma is cartesian
    dSigma[0] = surf.dSigma.tau;
    dSigma[1] = surf.dSigma.x;
    dSigma[2] = surf.dSigma.y;
    dSigma[3] = surf.dSigma.eta;
    
    gdsig = gmu[0]*dSigma[0] + gmu[1]*dSigma[1] + gmu[2]*dSigma[2] + gmu[3]*dSigma[3];  // GeV fm^3
    
    EdNd3p_ = EdNd3p_ + particle_.getDegen()*gdsig/(pow(2*M_PI,3)*pow(0.197,3)); 
  }   
  
  return EdNd3p_;     // GeV^[-2]
}

//**************************************************************************************************************

void FastReso::calc_observables(gsl_spline2d *splinef1, gsl_spline2d *splinef2, gsl_interp_accel *xacc, gsl_interp_accel *yacc, const double& mint, const double& maxt, const double& minp, const double& maxp) {

  double m2 = std::pow(particle_.getMass(),2);

  for (size_t pT_idx = 0; pT_idx < pTarr.size(); ++pT_idx) {
    double pT = pTarr[pT_idx]; 
    std::cout << pT << std::endl;
    for (size_t y_idx = 0; y_idx < yarr.size(); ++y_idx) {
      double mT = std::sqrt(m2 + pT*pT);
      double pz = pT*std::sinh(yarr[y_idx]);
      double y = std::asinh(pz/mT);
      double c2hy = std::pow(std::cosh(y), 2);
      double jac = std::sqrt(1-m2/(mT*mT*c2hy));
      // Integrate over phi
      for (size_t phi_idx = 0; phi_idx < phiarr.size(); ++phi_idx) {
        double phi = phiarr[phi_idx];
        //double EdNd3p_ = calc_EdNd3p(pT, y, phi, splinef1, splinef2, xacc, yacc, mint, maxt, minp, maxp);
        double EdNd3p_ = calc_EdNd3p_cart(pT, y, phi, splinef1, splinef2, xacc, yacc, mint, maxt, minp, maxp);
        EdNd3p[pT_idx][phi_idx][y_idx] = jac*EdNd3p_;
      }
    }
  }
}

//**************************************************************************************************************

void FastReso::output(const size_t& decayflag) {
  std::ostringstream filenameStream;
  // Compute only pT dependent quantities at rapidity 0
  if(decayflag==0){
    filenameStream << particle_.getName() << "_dNpTdpTdy_vn_thermal.dat";
  }else{
    filenameStream << particle_.getName() << "_dNpTdpTdy_vn_total.dat";
  }
  
  std::string filename = filenameStream.str();
  std::ofstream outputfile(filename);
  
  if (!outputfile.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }
  
  outputfile << "#pT \t\t dNpTdpTdy \t\t v1\t\t v2\t\t v3\t\t v4\n";
  size_t y_idx = (Ny-1)/2;
  
  for (size_t pT_idx = 0; pT_idx < pTarr.size(); ++pT_idx) {
    double v1cos = 0.0, v2cos = 0.0, v3cos = 0.0, v4cos = 0.0;
    double v1sin = 0.0, v2sin = 0.0, v3sin = 0.0, v4sin = 0.0;
    double den = 0.0;
    for (size_t phi_idx = 0; phi_idx < phiarr.size() - 1 ; ++phi_idx) {
      double delta_phi = phiarr[phi_idx + 1] - phiarr[phi_idx];
      
      v1cos += (std::cos(phiarr[phi_idx]) * EdNd3p[pT_idx][phi_idx][y_idx]
               + std::cos(phiarr[phi_idx + 1]) * EdNd3p[pT_idx][phi_idx + 1][y_idx]) * delta_phi/ 2.0;
      v2cos += (std::cos(2 * phiarr[phi_idx]) * EdNd3p[pT_idx][phi_idx][y_idx] 
               + std::cos(2 * phiarr[phi_idx + 1]) * EdNd3p[pT_idx][phi_idx + 1][y_idx]) * delta_phi/ 2.0;
      v3cos += (std::cos(3 * phiarr[phi_idx]) * EdNd3p[pT_idx][phi_idx][y_idx]
               + std::cos(3 * phiarr[phi_idx + 1]) * EdNd3p[pT_idx][phi_idx + 1][y_idx]) * delta_phi/ 2.0;
      v4cos += (std::cos(4 * phiarr[phi_idx]) * EdNd3p[pT_idx][phi_idx][y_idx] 
               + std::cos(4 * phiarr[phi_idx + 1]) * EdNd3p[pT_idx][phi_idx + 1][y_idx]) * delta_phi/ 2.0;
      v1sin += (std::sin(phiarr[phi_idx]) * EdNd3p[pT_idx][phi_idx][y_idx]
               + std::sin(phiarr[phi_idx + 1]) * EdNd3p[pT_idx][phi_idx + 1][y_idx]) * delta_phi/ 2.0;
      v2sin += (std::sin(2 * phiarr[phi_idx]) * EdNd3p[pT_idx][phi_idx][y_idx] 
               + std::sin(2 * phiarr[phi_idx + 1]) * EdNd3p[pT_idx][phi_idx + 1][y_idx]) * delta_phi/ 2.0;
      v3sin += (std::sin(3 * phiarr[phi_idx]) * EdNd3p[pT_idx][phi_idx][y_idx]
               + std::sin(3 * phiarr[phi_idx + 1]) * EdNd3p[pT_idx][phi_idx + 1][y_idx]) * delta_phi/ 2.0;
      v4sin += (std::sin(4 * phiarr[phi_idx]) * EdNd3p[pT_idx][phi_idx][y_idx] 
               + std::sin(4 * phiarr[phi_idx + 1]) * EdNd3p[pT_idx][phi_idx + 1][y_idx]) * delta_phi/ 2.0;
      den += (EdNd3p[pT_idx][phi_idx][y_idx] + EdNd3p[pT_idx][phi_idx + 1][y_idx]) * delta_phi / 2.0;
    }
    if (den != 0.0) {
      outputfile << pTarr[pT_idx] << "\t" << den << "\t" << v1cos/den << "\t" << v1sin/den
                 << "\t" << v2cos/den << "\t" << v2sin/den << "\t" << v3cos/den
                 << "\t" << v3sin/den << "\t" << v4cos/den << "\t" << v4sin/den << "\n";
    }
  }
  
  outputfile.close();
  
  // Compute only y dependent quantities by integrating over pT and phi
  filenameStream.str(""); // Clear the stream
  if(decayflag==0){
    filenameStream << particle_.getName() << "_dNdeta_vn_thermal.dat";
  }else{
    filenameStream << particle_.getName() << "_dNdeta_vn_total.dat";
  }

  filename = filenameStream.str();
  std::ofstream outputfile2(filename);
  
  if (!outputfile2.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }
  
  outputfile2 << "#eta\t\t dNdeta\t\t v1\t\t v2\t\t v3\t\t v4\n";
  
  for (size_t y_idx = 0; y_idx < yarr.size(); ++y_idx) {
    double v1cos = 0.0, v2cos = 0.0, v3cos = 0.0, v4cos = 0.0;
    double v1sin = 0.0, v2sin = 0.0, v3sin = 0.0, v4sin = 0.0;
    double den = 0.0;
    for (size_t pT_idx = 0; pT_idx < pTarr.size() - 1; ++pT_idx) {
      for (size_t phi_idx = 0; phi_idx < phiarr.size() - 1 ; ++phi_idx) {
        double delta_phi = phiarr[phi_idx + 1] - phiarr[phi_idx];
        double delta_pT = pTarr[pT_idx + 1] - pTarr[pT_idx];
      
        v1cos += (std::cos(phiarr[phi_idx]) * pTarr[pT_idx] * EdNd3p[pT_idx][phi_idx][y_idx]
                 + std::cos(phiarr[phi_idx + 1]) * pTarr[pT_idx] * EdNd3p[pT_idx][phi_idx + 1][y_idx]
                 + std::cos(phiarr[phi_idx]) * pTarr[pT_idx + 1] * EdNd3p[pT_idx + 1][phi_idx][y_idx]
                 + std::cos(phiarr[phi_idx + 1]) * pTarr[pT_idx + 1] * EdNd3p[pT_idx + 1][phi_idx + 1][y_idx]) * delta_pT * delta_phi/ 4.0;
        v2cos += (std::cos(2 * phiarr[phi_idx]) * pTarr[pT_idx] * EdNd3p[pT_idx][phi_idx][y_idx] 
                 + std::cos(2 * phiarr[phi_idx + 1]) * pTarr[pT_idx] * EdNd3p[pT_idx][phi_idx + 1][y_idx]
                 + std::cos(2 * phiarr[phi_idx]) * pTarr[pT_idx + 1] * EdNd3p[pT_idx + 1][phi_idx][y_idx] 
                 + std::cos(2 * phiarr[phi_idx + 1]) * pTarr[pT_idx + 1] * EdNd3p[pT_idx + 1][phi_idx + 1][y_idx])* delta_pT * delta_phi/ 4.0;
        v3cos += (std::cos(3 * phiarr[phi_idx]) * pTarr[pT_idx] * EdNd3p[pT_idx][phi_idx][y_idx]
                 + std::cos(3 * phiarr[phi_idx + 1]) * pTarr[pT_idx] * EdNd3p[pT_idx][phi_idx + 1][y_idx]
                 + std::cos(3 * phiarr[phi_idx]) * pTarr[pT_idx + 1] * EdNd3p[pT_idx + 1][phi_idx][y_idx]
                 + std::cos(3 * phiarr[phi_idx + 1]) * pTarr[pT_idx + 1] * EdNd3p[pT_idx + 1][phi_idx + 1][y_idx])* delta_pT * delta_phi/ 4.0;
        v4cos += (std::cos(4 * phiarr[phi_idx]) * pTarr[pT_idx] * EdNd3p[pT_idx][phi_idx][y_idx] 
                 + std::cos(4 * phiarr[phi_idx + 1]) * pTarr[pT_idx] * EdNd3p[pT_idx][phi_idx + 1][y_idx]
                 + std::cos(4 * phiarr[phi_idx]) * pTarr[pT_idx + 1] * EdNd3p[pT_idx + 1][phi_idx][y_idx] 
                 + std::cos(4 * phiarr[phi_idx + 1]) * pTarr[pT_idx + 1] * EdNd3p[pT_idx + 1][phi_idx + 1][y_idx])* delta_pT * delta_phi/ 4.0;
        v1sin += (std::sin(phiarr[phi_idx]) * pTarr[pT_idx] * EdNd3p[pT_idx][phi_idx][y_idx]
                 + std::sin(phiarr[phi_idx + 1]) * pTarr[pT_idx] * EdNd3p[pT_idx][phi_idx + 1][y_idx]
                 + std::sin(phiarr[phi_idx]) * pTarr[pT_idx + 1] * EdNd3p[pT_idx + 1][phi_idx][y_idx]
                 + std::sin(phiarr[phi_idx + 1]) * pTarr[pT_idx + 1] * EdNd3p[pT_idx + 1][phi_idx + 1][y_idx])* delta_pT * delta_phi/ 4.0;
        v2sin += (std::sin(2 * phiarr[phi_idx]) * pTarr[pT_idx] * EdNd3p[pT_idx][phi_idx][y_idx] 
                 + std::sin(2 * phiarr[phi_idx + 1]) * pTarr[pT_idx] * EdNd3p[pT_idx][phi_idx + 1][y_idx]
                 + std::sin(2 * phiarr[phi_idx]) * pTarr[pT_idx + 1] * EdNd3p[pT_idx + 1][phi_idx][y_idx] 
                 + std::sin(2 * phiarr[phi_idx + 1]) * pTarr[pT_idx + 1] * EdNd3p[pT_idx + 1][phi_idx + 1][y_idx])* delta_pT * delta_phi/ 4.0;
        v3sin += (std::sin(3 * phiarr[phi_idx]) * pTarr[pT_idx] * EdNd3p[pT_idx][phi_idx][y_idx]
                 + std::sin(3 * phiarr[phi_idx + 1]) * pTarr[pT_idx] * EdNd3p[pT_idx][phi_idx + 1][y_idx]
                 + std::sin(3 * phiarr[phi_idx]) * pTarr[pT_idx + 1] * EdNd3p[pT_idx + 1][phi_idx][y_idx]
                 + std::sin(3 * phiarr[phi_idx + 1]) * pTarr[pT_idx + 1] * EdNd3p[pT_idx + 1][phi_idx + 1][y_idx])* delta_pT * delta_phi/ 4.0;
        v4sin += (std::sin(4 * phiarr[phi_idx]) * pTarr[pT_idx] * EdNd3p[pT_idx][phi_idx][y_idx] 
                 + std::sin(4 * phiarr[phi_idx + 1]) * pTarr[pT_idx] * EdNd3p[pT_idx][phi_idx + 1][y_idx]
                 + std::sin(4 * phiarr[phi_idx]) * pTarr[pT_idx + 1] * EdNd3p[pT_idx + 1][phi_idx][y_idx] 
                 + std::sin(4 * phiarr[phi_idx + 1]) * pTarr[pT_idx + 1] * EdNd3p[pT_idx + 1][phi_idx + 1][y_idx])* delta_pT * delta_phi/ 4.0;
        den += (EdNd3p[pT_idx][phi_idx][y_idx] * pTarr[pT_idx] + EdNd3p[pT_idx][phi_idx + 1][y_idx] * pTarr[pT_idx]
                 + EdNd3p[pT_idx + 1][phi_idx][y_idx] * pTarr[pT_idx + 1] 
                 + EdNd3p[pT_idx + 1][phi_idx + 1][y_idx] * pTarr[pT_idx + 1])* delta_pT * delta_phi/ 4.0;
      }
    }
    if (den != 0.0) {
      outputfile2 << yarr[y_idx] << "\t" << den << "\t" << v1cos/den << "\t" << v1sin/den
                  << "\t" << v2cos/den << "\t" << v2sin/den << "\t" << v3cos/den
                  << "\t" << v3sin/den << "\t" << v4cos/den << "\t" << v4sin/den << "\n";
    }
  }  
  outputfile2.close();  
}
