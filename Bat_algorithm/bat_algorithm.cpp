#include "bat_algorithm.h"
#include <time.h>

//Checked
void Bat_algorithm::find_best_bat() {
	int best_bat = 0, current_bat = 1;

	for (; current_bat < population; current_bat++)
		if (function_values[current_bat] < function_values[best_bat])
			best_bat = current_bat;

	int best_bat_index = best_bat * dimension;
	for (int i = 0; i < dimension; i++)
		best_bat_values.push_back(bat_positions[best_bat_index + i]);	

	min_value = function_values[best_bat];
}

//Checked
void Bat_algorithm::bats_init() {
	for (int i = 0; i < population; i++) {
		freq.push_back(freq_min);

		srand(time(NULL));
		for (int j = 0; j < dimension; j++) {
			double rnd = static_cast<double> (rand()) / static_cast<double> (RAND_MAX);
			bat_velocities.push_back(0);
			bat_positions.push_back(lower_bound + (upper_bound - lower_bound) * rnd);
		}

		function_values.push_back(function(bat_positions.begin() + i*dimension, bat_positions.begin() + (i + 1)*dimension));
	}

	find_best_bat();
}

//Checked
void Bat_algorithm::cmp_position2bound(vector<double> &bat_positions) {
	for (auto it = bat_positions.begin(); it != bat_positions.end(); ++it)
		if (*it < lower_bound)
			*it = lower_bound;
		else
			if (*it > upper_bound)
				*it = upper_bound;
}


void Bat_algorithm::move_bats() {
	double rnd;
	double func_new;
	vector<double> curr_positions;

	bats_init();

	srand(time(NULL));
	for (int it = 0; it < iterations; ++it) {
		for (int i = 0; i < population; ++i) {
			rnd = static_cast<double> (rand()) / static_cast<double> (RAND_MAX);
			freq[i] = freq_min + (freq_max - freq_min) * rnd;

			int curr_bat = i*dimension;
			for (int j = 0; j < dimension; ++j) {
				bat_velocities[curr_bat + j] += (bat_positions[curr_bat + j] - best_bat_values[j])*freq[i];
				curr_positions[curr_bat + j] = bat_positions[curr_bat + j] + bat_velocities[curr_bat + j];

				cmp_position2bound(curr_positions);
			}

			//TODO generate local solution 
			rnd = static_cast<double> (rand()) / static_cast<double> (RAND_MAX);

			if (rnd > pulse_rate) {
				//random decimal in range [-1; 1]
				double loud_coef = static_cast<double> (rand()) / static_cast<double> (RAND_MAX / 2) - 1;
				for (int i = 0; i < dimension; ++i)
					best_bat_values[i] += loud_coef * loudness;
			}

			//TODO generate new solution by flying randomly
			vector<double> bat(dimension);
			//std::copy(bat_positions.begin() + i*dimension, bat_positions.begin() + (i + 1)*dimension, bat.begin());

			func_new = function(bat_positions.begin() + i*dimension, bat_positions.begin() + (i + 1)*dimension);
			rnd = static_cast<double> (rand()) / static_cast<double> (RAND_MAX);

			if (rnd < loudness && func_new < function_values[i]) {
				for (int j = 0; j < dimension; j++) 
					bat_positions[curr_bat + j] = curr_positions[curr_bat + j];
				function_values[i] = func_new;
				loudness *= alpha;
				pulse_rate *= (1 - exp(-alpha*it));
			}

			//TODO check func_new solution for the best global solution
			if (func_new < min_value) {
				for (int j = 0; j < dimension; ++j)
					best_bat_values[j] = curr_positions[curr_bat + j];
				min_value = func_new;
			}	
		}
	}
}