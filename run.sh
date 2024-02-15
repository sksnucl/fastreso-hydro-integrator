g++ -fopenmp -O3 -Wall -std=c++17 -o myprogram main.cpp particle.cpp freezeout.cpp fastreso.cpp -lgsl -lgslcblas -lm
time ./myprogram 1 ./fastresodata/ ./fodata/surfaceb5p0.dat
rm myprogram
