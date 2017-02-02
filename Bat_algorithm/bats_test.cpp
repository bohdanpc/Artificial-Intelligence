#include <iostream>
#include "Sphere.h"
#include "Ackley.h"
#include "Griewank.h"
#include "Rastrigin.h"
#include "Rosenbrock.h"

using namespace std;

int main() {
	Sphere sphere(30, 100, 50, 0, 1, -100, 100, 2);
	Ackley ackley(30, 100, 10, 0, 1, -32.768, 32.768, 2);
	Griewank griewank(50, 100, 50, 0, 1, -600, 600, 2);
	Rastrigin rastrigin(50, 100, 30, 0, 1, -5.12, 5.12, 2);
	Rosenbrock rosenbrock(50, 100, 30, 0, 1, -5, 10, 2);

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