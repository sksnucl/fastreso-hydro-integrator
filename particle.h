#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <string>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

class Particle {
public:
  //Particle();
  //~Particle();
  
  std::string getName() const;
  double getMass() const;
  int getStrangeness() const;
  int getCharmness() const;
  int getBottomness() const;
  int getBaryonNumber() const;
  int getCharge() const;
  int getPartId() const;
  int getDegen() const;
  double getIsospin3() const;
  
  void setName(std::string& name);
  void setMass(double mass);
  void setStrangeness(int strangeness);
  void setCharmness(int charmness);
  void setBottomness(int bottomness);
  void setBaryonNumber(int baryonNumber);
  void setCharge(int charge);
  void setPartId(int pid);
  void setDegen(int degen);
  void setIsospin3(double isospin);
  
  void read_fastreso_components(const std::string& folderpath, const std::vector<double>& temparr, const std::vector<double>& pbararr);
  void clear();
  
  double interpolated_f1(double temp_value, double pbar_value) const;
  double interpolated_f2(double temp_value, double pbar_value) const;
  
private:
  std::string name_;
  int pid_;
  double mass_;
  int strangeness_;
  int charmness_;
  int baryonNumber_;
  int charge_;
  double isospin3_;
  int bottomness_;
  int degen_;
  int partid_;
  
  std::vector<double> temp;
  std::vector<double> pbar;
  std::vector<std::vector<double> > f1;
  std::vector<std::vector<double> > f2;
  
  const gsl_interp2d_type *T = gsl_interp2d_bilinear;
  gsl_spline2d *splinef1 ;
  gsl_spline2d *splinef2 ;
  gsl_interp_accel *xacc ;
  gsl_interp_accel *yacc ;
  double maxpstar;
  double minpstar;
  double mintemp;
  double maxtemp;
};

#endif
