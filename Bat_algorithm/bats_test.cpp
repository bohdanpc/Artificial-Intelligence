#include <iostream>
#include "Sphere.h"
#include "Ackley.h"
#include "Griewank.h"
#include "Rastrigin.h"
#include "Rosenbrock.h"
#include "time.h"

using namespace std;

int main() {
	srand(time(NULL));
	Sphere sphere(50, 2000, 50, 0, 2, -100, 100, 0.9, 0.6, 0.7, 0.1);
	Ackley ackley(100, 2000, 10, 0, 1, -32.768, 32.768, 0.5, 0.1, 0.7, 0.1);
	Griewank griewank(50, 1000, 50, 0, 1, -600, 600, 0.9, 0.6, 0.7, 0.1);
	Rastrigin rastrigin(50, 1000, 30, 0, 1, -5.12, 5.12, 0.9, 0.6, 0.7, 0.1);
	Rosenbrock rosenbrock(50, 1000, 30, 0, 1, -5, 10, 0.9, 0.6, 0.7, 0.1);

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