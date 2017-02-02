#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include "bat_algorithm.h"

class Ackley : public Bat_algorithm {
public:
	Ackley(int population, int iterations, int dimension, double freq_min, double freq_max,
		double lower_bound, double upper_bound, double loudness_max = 2) :
		Bat_algorithm(population, iterations, dimension, freq_min, freq_max, lower_bound, upper_bound, loudness_max) 
	{
		a = 20;
		b = 0.2;
		c = 2 * M_PI;
	};

	double function(const double_iterator beg, const double_iterator end) {
		double sum1 = 0, sum2 = 0, left_op, right_op;

		for (vector<double>::iterator it = beg; it != end; it++) {
			sum1 += pow(*it, 2);
			sum2 += cos(c * *it);
		}
		left_op = -a * exp(-b * sqrt(sum1 / get_dimension()));
		right_op = -exp(sum2 / get_dimension()) + a + exp(1.0);

		return (left_op + right_op);
	}

private:
	double a, b, c;
};