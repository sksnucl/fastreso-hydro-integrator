#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <vector>
#include <string>
#include <algorithm>
#include "freezeout.h"
#include "fastreso.h"

void getArrays(const size_t& decayflag, const std::string& folderPath, std::vector<Particle>& particles,
                   std::vector<double>& templist, std::vector<double>& pbarlist) {

  bool firstfile = false;

  std::vector<int> pdglist;
  std::vector<std::string> filenames;
  std::string findstr;
  
  if(decayflag == 0){
    findstr = "thermal";
  }else{
    findstr = "total";
  }
  
  for (const auto& entry : std::filesystem::directory_iterator(folderPath)) {
    if (entry.is_regular_file()) {
      std::string filename = entry.path().filename().string();
      
      //if(filename.find("total") != std::string::npos) {
      if(filename.find(findstr) != std::string::npos) {
        filenames.push_back(filename);
	
	// Read the PDG id of particle from the filename and store it in pdglist
	size_t startpos = filename.find_first_of("0123456789");
	std::string pdgsubstr = filename.substr(startpos);
	int pdgid = std::stoi(pdgsubstr);
        
	auto ip = std::find(pdglist.begin(), pdglist.end(), pdgid);
        
	if(ip == pdglist.end()){
	  pdglist.push_back(pdgid) ;
	  std::ifstream filedata(folderPath+filename);
	  
	  std::string line;
	  std::getline(filedata, line);
	  std::getline(filedata, line);
	  
	  Particle tmp;
	  std::istringstream iss(line);
	  
	  std::string name;
	  int pid, deg, Qb, Qs, Qc, Qbot, charge, ndecays;
	  double mass, isospin, width, nyield; 
	  std::string dummychar;
	  
	  iss >> dummychar >> pid >> name >> mass >> width >> deg >> Qb >> Qs >> Qc >> Qbot >> isospin >> charge >> ndecays >> nyield ;
	  
	  tmp.setName(name);
	  tmp.setPartId(pid);
	  tmp.setMass(mass);
	  tmp.setDegen(deg);
	  tmp.setBaryonNumber(Qb);
	  tmp.setStrangeness(Qs);
	  tmp.setCharmness(Qc);
	  tmp.setBottomness(Qbot);
	  tmp.setIsospin3(isospin);
	  tmp.setCharge(charge);
	  
	  particles.push_back(tmp);
	  filedata.close();
	  
	}
	
	// Read the temperature from the filename and store it in templist
	startpos = filename.find("T0.");
	std::string tempsubstr = filename.substr(startpos+1); // Extract the substring starting from one place after T (which is 0)
	double var = std::stod(tempsubstr); // Convert the substring to double
	auto it = std::find(templist.begin(), templist.end(), var);
        
	if(it == templist.end()){
	  templist.push_back(var) ;
	}
        
	// Read the first file to fill the pbarlist
	if(firstfile == false){
	  std::ifstream filedata(folderPath+filename);
	  std::string line;
          
	  //Skip three lines
	  std::getline(filedata, line);
	  std::getline(filedata, line);
	  std::getline(filedata, line);
          
	  while (std::getline(filedata, line)) {
	    std::istringstream iss(line);
	    double pbarValue;
	    if (iss >> pbarValue) {
	      pbarlist.push_back(pbarValue) ;
	    }
	  }
          
	  firstfile = true ;
	}
	
      }
    }
  }
  std::sort(templist.begin(), templist.end());
}

// Function to create and initialize spline object
gsl_spline2d* create_spline(const std::vector<double>& xa, const std::vector<double>& ya,
                            const std::vector<double>& za, const gsl_interp2d_type *T) {
    size_t nx = xa.size();
    size_t ny = ya.size();
    gsl_spline2d *spline = gsl_spline2d_alloc(T, nx, ny);
    gsl_spline2d_init(spline, xa.data(), ya.data(), za.data(), nx, ny);
    return spline;
}

// main program

int main(int argc, char *argv[]) {
  
  if (argc < 4) {
    //std::cout << "You did not provide the folder path." << std::endl;
    std::cerr << "Not all input arguments provided." << std::endl;
    //std::cerr << "Syntax: " << argv[0] << " <folder_path>" << std::endl;
    return 1; 
  }
  
  size_t decayflag = std::stoi(argv[1]); // Flag for thermal or total
  std::string folderPath = argv[2];
  std::string fofile = argv[3];  // Name+Path of freezeout file
  
  // define the freezeout object
  Freezeout freezeout(fofile);
  
  // define array to store particles, temperature and pbar
  std::vector<Particle> particles;
  std::vector<double> templist;
  std::vector<double> pbarlist;
  
  getArrays(decayflag, folderPath, particles, templist, pbarlist);
  
  std::cout << "Number of particles: " << particles.size() << std::endl;
  std::cout << "Number of temp: " << templist.size() << std::endl;
  std::cout << "Number of pbar: " << pbarlist.size() << std::endl;
  
  const gsl_interp2d_type *T = gsl_interp2d_bilinear;
  size_t Ntemp = templist.size();
  size_t Npbar = pbarlist.size();
  
  double minpstar = pbarlist[0];
  double maxpstar = pbarlist[pbarlist.size() - 1];
  double mintemp = templist[0];
  double maxtemp = templist[templist.size() - 1];
  
  std::vector<double> zaf1(Ntemp * Npbar, 0.0);
  std::vector<double> zaf2(Ntemp * Npbar, 0.0);

  gsl_spline2d *splinef1;
  gsl_spline2d *splinef2;
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();
  
  for (size_t i = 0; i < particles.size(); ++i) {
    std::cout << particles[i].getPartId() << " " << particles[i].getName() <<std::endl;
    particles[i].read_fastreso_components(decayflag, folderPath, templist, pbarlist);

    for (size_t k = 0; k < Ntemp; ++k) {
      for (size_t l = 0; l < Npbar; ++l) {
        zaf1[l*Ntemp + k] = particles[i].f1[k][l];
        zaf2[l*Ntemp + k] = particles[i].f2[k][l];
        //gsl_spline2d_set(splinef1, zaf1.data(), k, l, particles[i].f1[k][l]);
        //gsl_spline2d_set(splinef2, zaf2.data(), k, l, particles[i].f2[k][l]);
      }
    }

    splinef1 = create_spline(templist, pbarlist, zaf1, T);
    splinef2 = create_spline(templist, pbarlist, zaf2, T);
    FastReso fastreso(particles[i] , freezeout);
    fastreso.calc_observables(splinef1, splinef2, xacc, yacc, mintemp, maxtemp, minpstar, maxpstar);
    fastreso.output(decayflag);
    
    gsl_spline2d_free(splinef1);
    gsl_spline2d_free(splinef2);
  }

  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
  return 0;
}
