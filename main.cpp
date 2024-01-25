#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <vector>
#include <string>
#include <algorithm>
#include "freezeout.h"
#include "fastreso.h"

int main(int argc, char *argv[]) {
  
  if (argc < 3) {
    std::cout << "You did not provide the folder path." << std::endl;
    std::cerr << "Syntax: " << argv[0] << " <folder_path>" << std::endl;
    return 1; 
  }
  
  std::string folderPath = argv[1];
  std::string fofile = argv[2];  // Name+Path of freezeout file
  
  // define the freezeout object
  Freezeout freezeout(fofile);
  
  // define array to store particles
  std::vector<Particle> particles;
 
  //vector arrays for file reading
  std::vector<int> pdglist;
  std::vector<double> templist;
  std::vector<double> pbarlist;
  std::vector<std::string> filenames;
  
  bool firstfile = false;
  
  for (const auto& entry : std::filesystem::directory_iterator(folderPath)) {
    if (entry.is_regular_file()) {
      std::string filename = entry.path().filename().string();
      
      if(filename.find("total") != std::string::npos) {
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
  
  std::cout << "Number of particles: " << pdglist.size() << std::endl;
  std::cout << "Number of temp: " << templist.size() << std::endl;
  std::cout << "Number of pbar: " << pbarlist.size() << std::endl;
  
  for (size_t i = 0; i < particles.size(); ++i) {
    Particle tmp = particles[i];
    //std::cout << tmp.getPartId() << " " << tmp.getName() <<std::endl;
    tmp.read_fastreso_components(folderPath, templist, pbarlist);
    //std::cout << "components read" <<std::endl;
    FastReso fastreso(tmp , freezeout);
    //std::cout << "calculating" <<std::endl;
    fastreso.calc_observables();
    fastreso.output();
    tmp.clear();
  }

  return 0;
}
