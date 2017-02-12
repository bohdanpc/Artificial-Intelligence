#pragma once
#include <math.h>
#include "bat_algorithm.h"

class Sphere : public Bat_algorithm {
public:
	Sphere() {};
	Sphere(int population, int iterations, int dimension, double freq_min, double freq_max,
		double lower_bound, double upper_bound, double loudness_max, double loudness_min,
		double pulse_rate_max, double pulse_rate_min) :
		Bat_algorithm(population, iterations, dimension, freq_min, freq_max, lower_bound, upper_bound,
			loudness_max, loudness_min, pulse_rate_max, pulse_rate_min) {}

	double function(const double_iterator beg, const double_iterator end) {
		double sum = 0;
		for (vector<double>::iterator it = beg; it != end; ++it) 
			sum += pow(*it, 2);

		return sum;
	}
};