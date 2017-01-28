#pragma once
#include "bat_algorithm.h"
class Sphere : public Bat_algorithm {
public:
	double function(const double_iterator beg, const double_iterator end) {
		double sum = 0;

		for (auto it = beg; it != end; ++it)
			sum += pow(*it, 2);

		return sum;
	}

	Sphere() : Bat_algorithm(20, 1000, 5, 0, 2, -10, 10, 0.5, 0.5) {}
};