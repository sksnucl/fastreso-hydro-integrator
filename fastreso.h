#ifndef FASTRESO_H
#define FASTRESO_H

#include <vector>
#include <cmath>
#include "particle.h"
#include "freezeout.h"

class FastReso {
public:
  FastReso(const Particle& particle, const Freezeout& freezeout);
  //~FastReso();
  
  void calc_observables();
  void output();
  
private:
  const Particle& particle_;
  const Freezeout& freezeout_;
  
  // Define the range for pT, phi, and y
  double pT_min = 0.0, pT_max = 3.0;
  double phi_min = 0.0, phi_max = 2 * M_PI;
  double y_min = -1.0, y_max = 1.0;
  
  // Specify the number of points in each dimension
  const size_t NpT = 20;  /* Number of points in pT */
  const size_t Nphi = 20; /* Number of points in phi */
  const size_t Ny = 20; /* Number of points in y */
  
  double deltapT, deltaphi, deltay;
  
  std::vector<double> pTarr;
  std::vector<double> phiarr;
  std::vector<double> yarr;
  
  std::vector<std::vector<std::vector<double> > > EdNd3p;
  std::vector<std::vector<double> > dNpTdpTdy;
  std::vector<std::vector<double> > v1;
  std::vector<std::vector<double> > v2;
  std::vector<std::vector<double> > v3;
  std::vector<std::vector<double> > v4;
  std::vector<double> dNdy;
  
  double calc_EdNd3p(const double& pT, const double& y, const double& phi);
  void calc_dNpTdpTdy();
  void calc_dNdy();
  void calc_vn();
};

#endif
