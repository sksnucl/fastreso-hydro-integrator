#ifndef FASTRESO_H
#define FASTRESO_H

#include <vector>
#include <cmath>
#include "particle.h"
#include "freezeout.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

class FastReso {
public:
  FastReso(const Particle& particle, const Freezeout& freezeout);
  //~FastReso();
  
  void calc_observables(gsl_spline2d *splinef1, gsl_spline2d *splinef2, gsl_interp_accel *xacc, gsl_interp_accel *yacc, const double& mint, const double& maxt, const double& minp, const double& maxp);
  void output(const size_t& decayflag);
  
private:
  const Particle& particle_;
  const Freezeout& freezeout_;
  
  // Define the range for pT, phi, and y
  double pT_min = 0.1, pT_max = 3.1;
  double phi_min = 0.0, phi_max = 2 * M_PI;
  double y_min = -5.0, y_max = 5.0;
  
  // Specify the number of points in each dimension
  const size_t NpT = 21;  /* Number of points in pT */
  const size_t Nphi = 21; /* Number of points in phi */
  const size_t Ny = 51; /* Number of points in y */
  
  double deltapT, deltaphi, deltay;
  
  std::vector<double> pTarr;
  std::vector<double> phiarr;
  std::vector<double> yarr;
  
  std::vector<std::vector<std::vector<double> > > EdNd3p;
/*  std::vector<std::vector<std::vector<double> > > EdNd3p(NpT, std::vector<std::vector<double>>(Nphi, std::vector<double>(Ny, 0.0)));*/
  std::vector<std::vector<double> > dNpTdpTdy;
  std::vector<std::vector<double> > v1;
  std::vector<std::vector<double> > v2;
  std::vector<std::vector<double> > v3;
  std::vector<std::vector<double> > v4;
  std::vector<double> dNdy;
  
  double calc_EdNd3p(const double& pT, const double& y, const double& phi, gsl_spline2d *splinef1, gsl_spline2d *splinef2, gsl_interp_accel *xacc, gsl_interp_accel *yacc, const double& mint, const double& maxt, const double& minp, const double& maxp);
  double calc_EdNd3p_cart(const double& pT, const double& y, const double& phi, gsl_spline2d *splinef1, gsl_spline2d *splinef2, gsl_interp_accel *xacc, gsl_interp_accel *yacc, const double& mint, const double& maxt, const double& minp, const double& maxp);

};

#endif
