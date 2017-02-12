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
//input value in range [0.3, 1.99], 
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
							 double lower_bound, double upper_bound, double loudness_max, double loudness_min, 
							 double pulse_rate_max, double pulse_rate_min) {
	this->population = population;
	this->iterations = iterations;
	this->dimension = dimension;
	this->freq_min = freq_min;
	this->freq_max = freq_max;
	this->lower_bound = lower_bound;
	this->upper_bound = upper_bound;

	this->loudness_max = loudness_max;
	this->loudness_max = loudness_min;
	this->pulse_rate_max = pulse_rate_max;
	this->pulse_rate_min = pulse_rate_min;

	mantegna_scale = (lambda_max - lambda_min) / (loudness_max - loudness_min);

	bat_best.resize(dimension);
	fitness_func.resize(population);
	bat_position.resize(population * dimension);

	bat_pulse_rate.resize(population, pulse_rate_min);
	bat_loudness.resize(population, loudness_max);
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
		if (fitness_func[current_bat] < fitness_func[best_bat_ind])
			best_bat_ind = current_bat;

	int best_bat_beg = best_bat_ind * dimension;
	for (int j = 0; j < dimension; j++)
		bat_best[j] = bat_position[best_bat_beg + j];	

	min_value = fitness_func[best_bat_ind];
}


void Bat_algorithm::bats_init() {
	double rnd;

	for (int curr_bat = 0; curr_bat < population; ++curr_bat) {
		rnd = static_cast<double> (rand()) / RAND_MAX;

		int curr_bat_idx = curr_bat * dimension;
		for (int j = 0; j < dimension; j++) {
			rnd = static_cast<double> (rand()) / RAND_MAX;
			bat_position[curr_bat_idx + j] = lower_bound + (upper_bound - lower_bound) * rnd;
		}

		fitness_func[curr_bat] = function(bat_position.begin() + curr_bat*dimension, 
										  bat_position.begin() + (curr_bat + 1)*dimension);
	}
	find_best_bat();
}

/*
* Description:
*	return random bat excluding the current one
*/
int Bat_algorithm::select_rand_bat(const int curr_bat) {
	int rand_bat = rand() % population;
	while (rand_bat == curr_bat)
		rand_bat = rand() % population;
	return rand_bat;
}

void Bat_algorithm::cmp_position2bound(const double_iterator begin, const double_iterator end) {
	for (vector<double>::iterator it = begin; it != end; ++it)
		if (*it < lower_bound)
			*it = lower_bound;
		else
			if (*it > upper_bound)
				*it = upper_bound;
}


void Bat_algorithm::update_bat_position(const int bat_num, vector<double> &one_bat_curr_pos) {
	int rand_bat = select_rand_bat(bat_num);
	int curr_bat_idx = bat_num * dimension;
	
	double f1 = freq_min + (freq_max - freq_min) * static_cast<double> (rand()) / RAND_MAX;
	double f2 = freq_min + (freq_max - freq_min) * static_cast<double> (rand()) / RAND_MAX;
	
	if (fitness_func[rand_bat] < fitness_func[bat_num])
		for (int j = 0; j < dimension; j++)
			one_bat_curr_pos[j] = bat_position[curr_bat_idx + j] +
			(bat_best[j] - bat_position[curr_bat_idx + j]) * f1 +
			(bat_position[rand_bat * dimension + j] - bat_position[curr_bat_idx + j])*f2;
	else
		for (int j = 0; j < dimension; j++)
			one_bat_curr_pos[j] = bat_position[curr_bat_idx + j] +
								 (bat_best[j] - bat_position[curr_bat_idx + j]) * f1;

	cmp_position2bound(one_bat_curr_pos.begin(), one_bat_curr_pos.end());
}


void Bat_algorithm::check_curr2best(const double func_new, const vector<double> &one_bat_pos) {
	if (func_new < min_value) {
		for (int j = 0; j < dimension; ++j)
			bat_best[j] = one_bat_pos[j];
		min_value = func_new;
	}
}


void Bat_algorithm::check_new_solution(const int bat_num, const int iter, vector<double> &one_bat_pos) {
	double rnd = static_cast<double> (rand()) / static_cast<double> (RAND_MAX);
	double func_new = function(one_bat_pos.begin(), one_bat_pos.end());

	//accept new solution of the bat if current sol is better than previous one
	if (rnd < bat_loudness[bat_num] && func_new <= fitness_func[bat_num]) {
		for (int j = 0; j < dimension; j++)
			bat_position[bat_num * dimension + j] = one_bat_pos[j];
		fitness_func[bat_num] = func_new;
		bat_loudness[bat_num] = (loudness_max / (1.0 - iterations)) * (iter - iterations);
		bat_pulse_rate[bat_num] = ((0. - 1) / (1. - iterations)) * (iter - iterations) + 1;
	
	}
	//compare current solution with the global best solution
	check_curr2best(func_new, one_bat_pos);
}


double Bat_algorithm::get_average_loudness() {
	double res = 0;
	for (int i = 0; i < population; i++)
		res += bat_loudness[i];
	return res;
}


/*
* Description:
*	lambda_min for mantegna_random is 0.3, lambda_max is 1.99;
*	mantegna_scale is evaluated in ctor
*   depending on loudness min and max values
*/
void Bat_algorithm::local_search(const int bat_num, vector<double> &one_bat_pos) {
	double rnd = static_cast<double> (rand()) / RAND_MAX;

	if (rnd > bat_pulse_rate[bat_num]) {
		double lambda = lambda_min + mantegna_scale * (bat_loudness[bat_num] - loudness_min);
		for (int j = 0; j < dimension; ++j)
			one_bat_pos[j] = bat_best[j] + mantegna_random(lambda) * get_average_loudness();
	}
}


void Bat_algorithm::move_bats() {
	vector<double> one_bat_pos(dimension);

	bats_init();

	for (int iter = 0; iter < iterations; ++iter) {
		for (int curr_bat = 0; curr_bat < population; ++curr_bat) {
			update_bat_position(curr_bat, one_bat_pos);
			local_search(curr_bat, one_bat_pos);
			check_new_solution(curr_bat, iter, one_bat_pos);
		}
	}
}