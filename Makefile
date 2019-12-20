CXX=g++
CXXFLAGS=-march=native -std=c++11



percolation:
	$(CXX) $(CXXFLAGS) sim_site_percolation.cpp -o percolation 

all: percolation

clean:
	rm -rf *.o percolation



