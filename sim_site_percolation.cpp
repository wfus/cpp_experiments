/**
 * Quick simulation for 2D site percolation.
 */
#include <random>
#include <chrono>
#include <algorithm>
#include <iostream>
#include <vector>
#include <stack>
#include <queue>
#include <tuple>

/**
 * Simulate the total resistance of a random resistor network.
 * Consider N points chosen randomly and uniformly in a 2d box of side sqrt(N).
 * Put a resistor between each two points with resistance proportional to e^(Cr)
 * where r is the distance between the points.
 * Find the average resistance of the system between the left-most
 * and right_most points.
 */
double sim_random_resistor_percolation(double N) {
	return 0.0;
}

/**
 * Simulate 2D site percolation.
 * Consider a square with side 1. N rods of length La nd negligible
 * width are thrown into the square, with their center uniformly chosen in the
 * square and their angle uniformly chosen.
 */

struct Point {
	double x;
	double y;
} typedef Point;


/* Check if three points are in CCW order */
bool ccw(Point A, Point B, Point C) {
    return (C.y - A.y) * (B.x - A.x) > (B.y - A.y) * (C.x - A.x);
}

/* Check if line defined by endpoints AB intersects with CD */
bool intersect(Point A, Point B, Point C, Point D) {
    return (ccw(A,C,D) != ccw(B,C,D)) && (ccw(A,B,C) != ccw(A,B,D));
}


bool sim_rod_percolation(int N, double L) {
	// Generate the rod centers and angles for our simulation.
	std::random_device rnd_device;
	std::mt19937 mersenne_engine {rnd_device()};
	// Uniform distribution is considered inclusive
	std::uniform_real_distribution<double> uniform_dist {0.0, 1.0};
	std::uniform_real_distribution<double> uniform_angle {0.0, 2.0 * M_PI};

	auto loc_gen = [&uniform_dist, &mersenne_engine]() {
		return uniform_dist(mersenne_engine);
	};

	auto angle_gen = [&uniform_angle, &mersenne_engine]() {
		return uniform_angle(mersenne_engine);
	};


	std::vector<double> xlocs(N), ylocs(N), angles(N);
	std::vector<Point> start_points, end_points;
	
	generate(begin(xlocs), end(ylocs), loc_gen);
	generate(begin(ylocs), end(ylocs), loc_gen);
	generate(begin(angles), end(angles), angle_gen);

	// turn this into start_points and end-points
	for (int idx = 0; idx < xlocs.size(); idx++) {
		double xloc = xlocs[idx];
		double yloc = ylocs[idx];
		double angle = angles[idx];
	}
	return false;
}


/**
 * Simulate 2D site percolation, with prob as the probability of
 * spreading. Spreading only works if two cells are directly adjacent.
 * N is the side length of the 2D array.
 */
bool sim_2d_site_percolation(int N, double prob) {
	std::vector<std::vector<int>> lattice(N, std::vector<int>(N, 0));
	std::random_device rnd_device;
	std::mt19937 mersenne_engine {rnd_device()};
	// Uniform distribution is considered inclusive
	std::uniform_int_distribution<int> dist {0, N};
	std::uniform_real_distribution<double> uniform_dist {0, 1};

	auto rand_gen = [&dist, &mersenne_engine]() {
		return dist(mersenne_engine);
	};

	// Pregenerate all of the random numbers we need: one approach
	// const int NUM_SAMPLES = 10000;
	// std::vector<int> rand_numbers(NUM_SAMPLES);
	// generate(begin(rand_numbers), end(rand_numbers), rand_gen);

	std::queue<std::tuple<int, int>> curr_stack;

	// Generate the first tree to be "open"
	for (int i = 0; i < N; i++) {
		if (uniform_dist(mersenne_engine) < prob) {
			lattice[i][0] = 1;
			curr_stack.push(std::make_tuple(i, 0));
		}
	}

	// Now, we go through our stack and keep propagating. If we
	// hit something that was already a 1, then we don't add it
	// back to the stack. We will return True if we ever hit
	// the last row - otherwise we will return False
	
	while (!curr_stack.empty()) {
		auto indices = curr_stack.front();
		curr_stack.pop();
		int i = std::get<0>(indices);
		int j = std::get<1>(indices);

		if (j == (N - 1)) return true;

		// Otherwise, loop through all neighbors. See the
		// chance that it spreads and only add to queue if
		// we have not seen it before.
		for (auto di: {-1, 1}) {
			for (auto dj: {-1, 1}) {
				if (i + di < 0 || i + di >= N) continue;	
				if (j + dj < 0 || j + dj >= N) continue;	
				if (lattice[i + di][j + dj] == 1) continue;
					
				if (uniform_dist(mersenne_engine) < prob) {
					lattice[i + di][j + dj] = 1;
					// std::cout << "Adding tuple " << i + di << " " << j + dj
					// << std::endl;
					curr_stack.push(std::make_tuple(i + di, j + dj));
				}

			}	
		}


	}
	

	return false;
}



int main(int argc, char** argv) {


	// sim_rod_percolation(100, 0.1);


	const int NUM_EXPERIMENTS = 100;
	std::vector<double> probs({0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
			0.8, 0.9, 1.0});
		
#pragma omp parallel for
	for (int i = 0; i < probs.size(); i++) {
		double prob = probs[i];
		int ctr = 0;
		for (int step = 0; step < NUM_EXPERIMENTS; step++) {
			if (sim_2d_site_percolation(200, prob)) {
				ctr ++;
			}
		}
		std::cout << "Prob: " << prob << ": " << (ctr / (float) NUM_EXPERIMENTS)
			  << std::endl;
	}
	return 0;
}
