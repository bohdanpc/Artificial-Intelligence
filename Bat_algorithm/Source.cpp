#include <iostream>
#include <time.h>
#include "bat_algorithm.h"
#include "Sphere.h"

using namespace std;

int main() {
	Sphere *sp1;
	srand(time(NULL));
	for (int i = 0; i < 10; i++) {
		sp1 = new Sphere;
		sp1->move_bats();
		cout << sp1->get_result() << endl;
	}
	return 0;
}