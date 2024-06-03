#ifndef FREEZEOUT_H
#define FREEZEOUT_H

#include <vector>
#include <string>
#include <stdexcept>
#include <array>
#include "particle.h"

struct FourVector {
  double tau, x, y, eta ;
  
  double getComponent(int i) const {
    switch (i) {
    case 0: return tau;
    case 1: return x;
    case 2: return y;
    case 3: return eta;
    default: 
      throw std::out_of_range("Wrong vector index");
    }
  }
};

struct SymmetricTensor {
  std::array<double, 10> element;
  
  double getComponent(int i, int j) {
    if (i < 0 || i > 3 || j < 0 || j > 3) {
      throw std::out_of_range("Invalid indices in symmmetric tensor");
    }
    
    if (i > j) {
      std::swap(i, j);
    }
    
    return element[i * (7 - i) / 2 + j];
  }
};

struct surface_element{
  FourVector xmu;
  FourVector dSigma;
  FourVector umu;
  double temperature;
  double energy;
  double entropy;
  double muB;
  double muS;
  double muQ;
  double bulk;
  SymmetricTensor shear;
  double dbeta[4][4];
};

class Freezeout {
public:
  Freezeout(const std::string& filename);
    
  std::vector<surface_element> mysurface;
  size_t NumberofCells;
    
private:
  std::string filename_;
  void readFile(const std::string& dataFile);
  void readFileCart(const std::string& datafile);
};

#endif
