#include "bat_algorithm.h"
#include <time.h>
#include <iostream>

using namespace std;

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
	double rnd = 0;
	for (int i = 0; i < population; i++) {
		freq.push_back(freq_min);

		for (int j = 0; j < dimension; j++) {
			rnd = static_cast<double> (rand()) / static_cast<double> (RAND_MAX);
			bat_velocities.push_back(0);
			bat_positions.push_back(lower_bound + (upper_bound - lower_bound) * rnd);
		}

		bat_loudness.push_back(loudness_max);
		bat_pulse_rate.push_back(pulse_rate_max);
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


void Bat_algorithm::check_curr_best(const double func_new, const vector<double> &curr_positions) {
	if (func_new < min_value) {
		for (int j = 0; j < dimension; ++j)
			best_bat_values[j] = curr_positions[j];
		min_value = func_new;
	}
}

//CHECKED
void Bat_algorithm::adjust_bat_parameters(int bat_num, vector<double> &curr_positions) {
	double rnd = static_cast<double> (rand()) / static_cast<double> (RAND_MAX);
	freq[bat_num] = freq_min + (freq_max - freq_min) * rnd;				//CHECK !!!
																		//with '-' solution seems a little bit better
	int curr_bat = bat_num*dimension;
	for (int j = 0; j < dimension; ++j) {
		bat_velocities[curr_bat + j] += (bat_positions[curr_bat + j] - best_bat_values[j])*freq[bat_num];
		curr_positions[j] = bat_positions[curr_bat + j] + bat_velocities[curr_bat + j];
																		//HERE TOO better
		cmp_position2bound(curr_positions);
	}
}


void Bat_algorithm::move_bats() {
	double rnd;
	double func_new;
	vector<double> curr_positions(dimension);
	int tmp = 0;

	bats_init();

	for (int it = 0; it < iterations; ++it) {
		for (int i = 0; i < population; ++i) {

			adjust_bat_parameters(i, curr_positions);

			//TODO generate local solution 
			rnd = static_cast<double> (rand()) / static_cast<double> (RAND_MAX);

			if (rnd > pulse_rate_max /*bat_pulse_rate[i]*/) {
				//random decimal in range [-1; 1]
				double loud_coef = static_cast<double> (rand()) / static_cast<double> (RAND_MAX / 2) - 1;
				for (int j = 0; j < dimension; ++j)
					curr_positions[j] = best_bat_values[j] + loud_coef * 0.001; //bat_loudness[i];				//Something is not ok too
			}

			//TODO generate new solution by flying randomly
			
			func_new = function(curr_positions.begin(), curr_positions.end());

			rnd = static_cast<double> (rand()) / static_cast<double> (RAND_MAX);// * loudness_max;

			int curr_bat = i * dimension;
			if (rnd < loudness_max /*bat_loudness[i]*/ && func_new <= function_values[i]) {
				for (int j = 0; j < dimension; j++) 
					bat_positions[curr_bat + j] = curr_positions[j];
				function_values[i] = func_new;
				bat_loudness[i] *= alpha;
				bat_pulse_rate[i] *= (1 - exp(-alpha*it));
			}

			//TODO check func_new solution for the best global solution
			check_curr_best(func_new, curr_positions);
		}
	}
}