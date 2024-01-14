#include <iostream>
#include "freezeout.h"
#include "fastreso.h"

int main(int argc, char **argv) {

    Freezeout freezeout("./freezeout.dat");
    
    Particle pion("pi+", 0.139, 0, 0, 0, 1, 2, 0);
    Particle proton("p", 1.0, 0, 0, 1, 1, 1, 1);
    Particle kaon("k+", 0.494, 1, 0, 0, 1, 1, 0);
    
    const int Ntemp = 10;
    const int Npbar = 10;
    
    pion.readf1f2(Ntemp, Npbar);
    proton.readf1f2(Ntemp, Npbar);
    kaon.readf1f2(Ntemp, Npbar);
    
    FastReso fastreso_pion(pion, freezeout);
    FastReso fastreso_proton(proton, freezeout);
    FastReso fastreso_kaon(kaon, freezeout);
    
    fastreso_pion.calc_observables();
    fastreso_proton.calc_observables();
    fastreso_kaon.calc_observables();
    
    fastreso_pion.output();
    fastreso_proton.output();
    fastreso_kaon.output();

 return 0;
}
