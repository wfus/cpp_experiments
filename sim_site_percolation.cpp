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
#include <set>
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

	generate(begin(xlocs), end(xlocs), loc_gen);
	generate(begin(ylocs), end(ylocs), loc_gen);
	generate(begin(angles), end(angles), angle_gen);

	// turn this into start_points and end-points
	for (int idx = 0; idx < xlocs.size(); idx++) {
		double xloc = xlocs[idx];
		double yloc = ylocs[idx];
		double angle = angles[idx];
		double dx = L * std::cos(angle);
		double dy = L * std::sin(angle);

		Point start_point, end_point;
		start_point.x = xloc - dx;
		end_point.x = xloc + dx;
		start_point.y = yloc - dy;
		end_point.y = yloc + dy;

		start_points.push_back(start_point);
		end_points.push_back(end_point);
	}

	std::vector<std::vector<int>> connections(N, std::vector<int>({}));
	// Find out which rods are connected to each other.
	Point start_a, start_b, end_a, end_b;
	for (int a = 0; a < N; a++) {
		for (int b = 0; b < N; b++) {
			start_a = start_points[a];
			end_a = end_points[a];
			start_b = start_points[b];
			end_b = end_points[b];
			if (intersect(start_a, end_a, start_b, end_b)) {
				connections[a].push_back(b);
			}
		}
	}

	std::vector<int> source_idxs;
	std::set<int> sink_idxs;
	for (int i = 0; i < start_points.size(); i++) {
		start_a = start_points[i];
		end_a = end_points[i];
		if (start_a.x < 0 || end_a.x < 0) {
			source_idxs.push_back(i);	
		} else if (start_a.x > 1. || end_a.x > 1.) {
			sink_idxs.insert(i);	
		}
	}

	for (int start_idx: source_idxs) {
		std::queue<int> q;
		std::set<int> seen;
		q.push(start_idx);

		while (!q.empty()) {
			int curr = q.front();
			q.pop();
			if (seen.count(curr) < 1) {
				seen.insert(curr);
				if (sink_idxs.count(curr) > 0) {
					return true;
				}
				for (int neighbor: connections[curr]) {
					q.push(neighbor);
				}
			}
		}
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


template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}


int main(int argc, char** argv) {

	const int NUM_EXPERIMENTS = 100;

	std::vector<double> n_choices = linspace(1.0, 100.0, 100);
	std::vector<double> prob_choices = linspace(0.0, 1.0, 101);

	for (int i = 0; i < n_choices.size(); i++)  {
		for (int j = 0; j < prob_choices.size(); j++)  {
			double n = n_choices[i];
			double p = prob_choices[i];
			int ctr = 0;

			#pragma omp parallel for
			for (int i = 0; i < NUM_EXPERIMENTS; i++) {
				if (sim_rod_percolation(n, p)) {
					ctr++;
				}
			}

			// Output it as N,P,AVERAGE_SUCCESS.
			double avg_success = (double) ctr / NUM_EXPERIMENTS;
			std::cout << n << "," << 
				p << "," << 
				avg_success << std::endl;
		}
	}
	
	
	
	/*
	const int NUM_EXPERIMENTS = 100;
	std::vector<double> probs;
	for (double prob = 0.40; prob < 0.55; prob += 0.005) {
		probs.push_back(prob);
	}
		
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
	*/
	return 0;
}

