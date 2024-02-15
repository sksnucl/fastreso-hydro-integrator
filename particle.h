#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <string>

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
  
  void read_fastreso_components(const size_t& decayflag, const std::string& folderpath, const std::vector<double>& temparr, const std::vector<double>& pbararr);
  
  std::vector<std::vector<double> > f1;
  std::vector<std::vector<double> > f2;
  
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
  
};

#endif
