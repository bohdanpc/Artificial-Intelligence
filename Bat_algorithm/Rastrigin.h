#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include "bat_algorithm.h"

class Rastrigin : public Bat_algorithm {
public:
	Rastrigin(int population, int iterations, int dimension, double freq_min, double freq_max,
		double lower_bound, double upper_bound, double loudness_max, double loudness_min,
		double pulse_rate_max, double pulse_rate_min) :
		Bat_algorithm(population, iterations, dimension, freq_min, freq_max, lower_bound, upper_bound,
			loudness_max, loudness_min, pulse_rate_max, pulse_rate_min) {};

	double function(const double_iterator beg, const double_iterator end) {
		double sum = 0;
		for (double_iterator it = beg; it != end; ++it) {
			sum += pow(*it, 2) - 10 * cos(2 * M_PI * *it);
		}
		return 10 * get_dimension() + sum;
	}
};