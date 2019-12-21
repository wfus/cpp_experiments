CXX=g++
CXXFLAGS=-march=native -std=c++11 -fopenmp -O3



percolation: sim_site_percolation.cpp
	$(CXX) $(CXXFLAGS) sim_site_percolation.cpp -o percolation 

all: percolation

clean:
	rm -rf *.o percolation



