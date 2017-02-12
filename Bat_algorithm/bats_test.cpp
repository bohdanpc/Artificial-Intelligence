#include <iostream>
#include "time.h"
#include "Sphere.h"
#include "Ackley.h"
#include "Griewank.h"
#include "Rastrigin.h"
#include "Rosenbrock.h"

using std::cout;
using std::endl;

int main() {
	srand((unsigned int) time(NULL));
	//population,  iterations,  dimension,    freq_min,    freq_max,
	//lower_bound, upper_bound, loudness_max, loudness_min, 
	//pulse_rate_max, pulse_rate_min

	Sphere sphere(30, 500, 50, 0, 2, -100, 100, 1, 0.5, 0.7, 0.1);
	Ackley ackley(30, 1000, 10, 0, 1, -32.768, 32.768, 0.5, 0.1, 0.7, 0.1);
	Griewank griewank(50, 1000, 2, 0, 1, -600, 600, 0.9, 0.6, 0.7, 0.1);	//50
	Rastrigin rastrigin(50, 1000, 2, 0, 1, -5.12, 5.12, 0.9, 0.6, 0.7, 0.1);	//30
	Rosenbrock rosenbrock(50, 1000, 2, 0, 1, -5, 10, 0.9, 0.6, 0.7, 0.1);	//30

	sphere.move_bats();
	cout << "Sphere function, minimum: " << sphere.get_result() << endl;

	ackley.move_bats();
	cout << "Ackley function, minimum: " << ackley.get_result() << endl;

	griewank.move_bats();
	cout << "Griewank function, minimum: " << griewank.get_result() << endl;

	rastrigin.move_bats();
	cout << "Rastrigin function, minimum: " << rastrigin.get_result() << endl;

	rosenbrock.move_bats();
	cout << "Rosenbrock function, minimum: " << rosenbrock.get_result() << endl;

	return 0;
}