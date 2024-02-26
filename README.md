# fastreso-hydro-integrator
g++ -std=c++11 -o myprogram main.cpp particle.cpp freezeout.cpp fastreso.cpp -I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib -lgsl -lgslcblas -lm
