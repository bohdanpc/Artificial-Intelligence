#define _USE_MATH_DEFINES
#include <cmath>
#include "bat_algorithm.h"
#include <time.h>
#include <iostream>

using namespace std;

double gaussian_random(double mue, double sigma) {
	double x1, x2, w, y;

	do {
		x1 = 2.0 * rand() / (RAND_MAX + 1) - 1.0;
		x2 = 2.0 * rand() / (RAND_MAX + 1) - 1.0;
		w = x1 * x1 + x2 * x2;
	} while (w >= 1.0);
	double llog = log(w);
	w = sqrt((-2.0 * llog) / w);
	y = x1 * w;

	return mue + sigma * y;
}

// tgammal prototype is in <math.h>
//lambda in range [0.3, 1.99], 
double mantegna_random(double lambda) {
	long double sigmaX, sigmaY = 1;
	double x, y;
	long double tg = tgammal(lambda + 1);
	sigmaX = tgammal(lambda + 1) * sin(M_PI * lambda / 2);
	double divider = tgammal((lambda) / 2) * lambda * pow(2., (lambda - 1) / 2);
	sigmaX /= divider;
	sigmaX = pow(sigmaX, (long double)1. / lambda);

	x = gaussian_random(0, sigmaX);
	y = fabs(gaussian_random(0, 1.));

	return x / pow(y, 1. / lambda);

}


Bat_algorithm::Bat_algorithm(int population, int iterations, int dimension, double freq_min, double freq_max, 
							 double lower_bound, double upper_bound, double loudness_max) {
	this->population = population;
	this->iterations = iterations;
	this->dimension = dimension;
	this->freq_min = freq_min;
	this->freq_max = freq_max;
	this->lower_bound = lower_bound;
	this->upper_bound = upper_bound;

	this->loudness_max = loudness_max;
	this->alpha = 0.99;
}

void Bat_algorithm::set_alpha(const double alpha) {
	this->alpha = alpha;
}

double Bat_algorithm::get_result() {
	return min_value;
}

double Bat_algorithm::get_dimension() {
	return dimension;
}

void Bat_algorithm::find_best_bat() {
	int best_bat_ind = 0, current_bat = 1;

	for (; current_bat < population; current_bat++)
		if (function_values[current_bat] < function_values[best_bat_ind])
			best_bat_ind = current_bat;

	int best_bat_beg = best_bat_ind * dimension;
	for (int j = 0; j < dimension; j++)
		best_bat[j] = bats_positions[best_bat_beg + j];	

	min_value = function_values[best_bat_ind];
}


void Bat_algorithm::bats_init() {
	double rnd, loudness;
	for (int i = 0; i < dimension; i++)
		best_bat.push_back(lower_bound);

	for (int curr_bat = 0; curr_bat < population; ++curr_bat) {
		rnd = static_cast<double> (rand()) / static_cast<double> (RAND_MAX);
		freq.push_back(freq_min + (freq_max - freq_min) * rnd);

		for (int j = 0; j < dimension; j++) {
			rnd = static_cast<double> (rand()) / static_cast<double> (RAND_MAX);
			bat_velocities.push_back(0);
			bats_positions.push_back(lower_bound + (upper_bound - lower_bound) * rnd);
		}

		//initial loundess value in range [loudness_max - 1, loudness_max]
		rnd = static_cast<double> (rand()) / static_cast<double> (RAND_MAX);
		loudness = loudness_max - rnd;
		
		bats_loudness.push_back(loudness);

		//initial pulse rate is close to 0
		rnd = static_cast<double> (rand()) / static_cast<double> (RAND_MAX) * 0.01;
		bats_pulse_rate.push_back(rnd);
		function_values.push_back(function(bats_positions.begin() + curr_bat*dimension, 
										   bats_positions.begin() + (curr_bat + 1)*dimension));
	}
	find_best_bat();
}


void Bat_algorithm::cmp_position2bound(const double_iterator begin, const double_iterator end) {
	for (vector<double>::iterator it = begin; it != end; ++it)
		if (*it < lower_bound)
			*it = lower_bound;
		else
			if (*it > upper_bound)
				*it = upper_bound;
}


void Bat_algorithm::adjust_bat_parameters(const int bat_num, vector<double> &curr_positions) {
	int curr_bat_beg = bat_num * dimension;

	for (int j = 0; j < dimension; ++j) {
		bat_velocities[curr_bat_beg + j] += (bats_positions[curr_bat_beg + j] - best_bat[j])*freq[bat_num];
		curr_positions[j] = bats_positions[curr_bat_beg + j] + bat_velocities[curr_bat_beg + j];		//Here better too
	}
	//cmp_position2bound(curr_positions);
}


double Bat_algorithm::get_average_loudness() {
	double res = 0;
	for (int i = 0; i < population; i++)
		res += bats_loudness[i];
	return res;
}


void Bat_algorithm::gen_local_sol_around_best(const int curr_bat, vector<double> &curr_bat_positions) {
	double rnd = static_cast<double> (rand()) / static_cast<double> (RAND_MAX);

	if (rnd > bats_pulse_rate[curr_bat]) {
		//random decimal in range [-1; 1]
		double loud_coef = static_cast<double> (rand()) / static_cast<double> (RAND_MAX / 2) - 1;
		for (int j = 0; j < dimension; ++j)
			curr_bat_positions[j] = best_bat[j] + loud_coef * get_average_loudness();	
	}
}


void Bat_algorithm::check_curr_to_best(const double func_new, const vector<double> &curr_positions) {
	if (func_new < min_value) {
		for (int j = 0; j < dimension; ++j)
			best_bat[j] = curr_positions[j];
		min_value = func_new;
	}
}


void Bat_algorithm::check_new_solution(const int curr_bat, const int iteration, vector<double> &curr_bat_positions) {
	int curr_bat_beg = curr_bat * dimension;
	double func_new = function(curr_bat_positions.begin(), curr_bat_positions.end());
	//random value in range [0, loudness_max]
	double rnd = static_cast<double> (rand()) / static_cast<double> (RAND_MAX) * loudness_max;

	if (rnd < bats_loudness[curr_bat] && func_new <= function_values[curr_bat]) {
		for (int j = 0; j < dimension; j++)
			bats_positions[curr_bat + j] = curr_bat_positions[j];
		function_values[curr_bat] = func_new;
		bats_loudness[curr_bat] *= alpha;
		bats_pulse_rate[curr_bat] *= (1 - exp(-alpha * iteration));
	}

	//check for the best global solution
	//check_curr_to_best(func_new, curr_bat_positions);
}


void Bat_algorithm::move_bats() {
	vector<double> curr_bat_positions(dimension);

	bats_init();

	for (int iteration = 0; iteration < iterations; ++iteration) {
		for (int curr_bat = 0; curr_bat < population; ++curr_bat) {
			int rand_bat = rand() % population;
			while (rand_bat == curr_bat)
				rand_bat = rand() % population;

			double f1 = freq_min + (freq_max - freq_min) * static_cast<double> (rand()) / static_cast<double> (RAND_MAX);
			double f2 = freq_min + (freq_max - freq_min) * static_cast<double> (rand()) / static_cast<double> (RAND_MAX);

			int curr_bat_beg = curr_bat * dimension;
			if (function_values[rand_bat] < function_values[curr_bat])
				for (int j = 0; j < dimension; j++)
					curr_bat_positions[j] = bats_positions[curr_bat_beg + j] +
					(best_bat[j] - bats_positions[curr_bat_beg + j]) * f1 +
					(bats_positions[rand_bat * dimension + j] - bats_positions[curr_bat_beg + j])*f2;
			else
				for (int j = 0; j < dimension; j++)
					curr_bat_positions[j] = bats_positions[curr_bat_beg + j] + (best_bat[j] - bats_positions[curr_bat_beg + j]) * f1;

			cmp_position2bound(curr_bat_positions.begin(), curr_bat_positions.end());
			//END ADJUSTING

			//BEGIN LOCAL SEARCH PART
			double rnd = static_cast<double> (rand()) / static_cast<double> (RAND_MAX);
			if (rnd > bats_pulse_rate[curr_bat]) {
				//random decimal in range [-1; 1]
				double loud_coef = static_cast<double> (rand()) / static_cast<double> (RAND_MAX / 2) - 1;
				double wi0 = (upper_bound - lower_bound) / 4.0;
				double wi1 = wi0 / 100.0;
				double wi = (wi0 - wi1) / (1 - iterations) * (iteration - iterations) + wi1;

				for (int j = 0; j < dimension; ++j)
					curr_bat_positions[j] = best_bat[j] + loud_coef * wi * get_average_loudness();
			}
			//END LOCAL SEARCH PART

			//CHECK NEW SOLUTION
			curr_bat_beg = curr_bat * dimension;
			double func_new = function(curr_bat_positions.begin(), curr_bat_positions.end());
			//random value in range [0, loudness_max]
			rnd = static_cast<double> (rand()) / static_cast<double> (RAND_MAX) * loudness_max;

			if (rnd < bats_loudness[curr_bat] && func_new <= function_values[curr_bat]) {
				for (int j = 0; j < dimension; j++)
					bats_positions[curr_bat + j] = curr_bat_positions[j];
				function_values[curr_bat] = func_new;
				bats_loudness[curr_bat] = loudness_max / (1.0 - iterations) * (iteration - iterations);
				bats_pulse_rate[curr_bat] = (0. - 1) / (1. - iterations) * (iteration - iterations) + 1;
			}

			if (function_values[curr_bat] < min_value) {
				for (int j = 0; j < dimension; j++)
					best_bat[j] = curr_bat_positions[j];
				min_value = func_new;
			}
			
		}
			//END CHECK NEW SOLUTION

		//find_best_bat();
	}
}