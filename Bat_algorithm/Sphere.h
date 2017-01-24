#pragma once
#include "bat_algorithm.h"
class Sphere : public Bat_algorithm {
public:
	double function(const double_iterator beg, const double_iterator end) {
		double sum = 0;
		double tmp;
		for (auto it = beg; it != end; ++it) 
			sum += pow(*it, 2);

		return sum;
	}

	Sphere() : Bat_algorithm(10, 20, 50, 1, 100, -100, 100, 0.8, 0.8) {}
};