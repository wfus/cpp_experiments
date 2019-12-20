/**
 * Quick simulation for 2D site percolation.
 */
#include <random>
#include <algorithm>
#include <iostream>
#include <vector>


double sim_2d_site_percolation(int N) {
	random_device rnd_device;
	mt19937 mersenne_engine {rnd_device()};
	uniform_int_distribution<int> dist {0, N};

	auto rand_gen = [&dist, &mersenne_engine]() {
		return dist(mersenne_engine);
	};

	const int NUM_SAMPLES = 10000;

	std::vector<int> rand_numbers(NUM_SAMPLES);
	generate(begin(rand_numbers), end(rand_numbers), rand_gen);
	
	for (auto i : rand_numbers) {
		std::cout << i << " ";
	}
}


int main() {
	sim_2d_site_percolation(100);
	std::cout << "Results were: " << std::endl;
	return 0;
}
