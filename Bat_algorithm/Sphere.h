#pragma once
#include "bat_algorithm.h"
class Sphere : public Bat_algorithm {
public:
	double function(double_iterator start_pos, double_iterator end_pos) {
		double sum = 0;
		for (auto it = start_pos; it != end_pos; ++it)
			sum += pow(*it, 2);

		return sum;
	}

	Sphere() : Bat_algorithm(50, 50, 0, 1, -100, 100, 0.8, 0.8) {}
	

};