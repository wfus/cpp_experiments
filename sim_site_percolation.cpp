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
	std::uniform_int_distribution<int> uniform_dist {0, 1};

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


int main() {
	const int NUM_EXPERIMENTS = 100;
	std::vector<double> probs({0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0});
		
#pragma omp parallel for
	for (int i = 0; i < probs.size(); i++) {
		double prob = probs[i];
		int ctr = 0;
		for (int step = 0; step < NUM_EXPERIMENTS; step++) {
			if (sim_2d_site_percolation(2000, prob)) {
				ctr ++;
			}
		}
		std::cout << "Prob: " << prob << ": " << (ctr / (float) NUM_EXPERIMENTS)
			  << std::endl;
	}
	return 0;
}
