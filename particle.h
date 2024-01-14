#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <string>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

class Particle {
public:
  Particle(const std::string& name, double mass, int strangeness, int charmness, int baryonNumber, int charge, int isospin3, int spin);
  ~Particle();
  
  std::string getName() const;
  double getMass() const;
  int getStrangeness() const;
  int getCharmness() const;
  int getBaryonNumber() const;
  int getCharge() const;
  int getSpin() const;
  int getIsospin3() const;
  void readf1f2(size_t Ntemp,size_t Npbar); 
  
  double interpolated_f1(double temp_value, double pbar_value) const;
  double interpolated_f2(double temp_value, double pbar_value) const;
  
private:
  std::string name_;
  double mass_;
  int strangeness_;
  int charmness_;
  int baryonNumber_;
  int charge_;
  int isospin3_;
  int spin_;
  int degen_;
  
  std::vector<double>* temp;   // Pointer to dynamically allocated vector
  std::vector<double>* pbar;   // Pointer to dynamically allocated vector
  std::vector<std::vector<double> >* f1;
  std::vector<std::vector<double> >* f2;
  
  const gsl_interp2d_type *T = gsl_interp2d_bilinear;
  double *zaf1;
  double *zaf2;
  gsl_spline2d *splinef1 ;
  gsl_spline2d *splinef2 ;
  gsl_interp_accel *xacc ;
  gsl_interp_accel *yacc ;
  size_t nx;
  size_t ny;
  double maxpstar;
  double minpstar;
  const double mintemp = 0.11;
  const double maxtemp = 0.16;
};

#endif
