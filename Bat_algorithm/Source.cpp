#include <iostream>
#include <time.h>
#include "bat_algorithm.h"
#include "Sphere.h"
using namespace std;

int main() {
	Sphere sp1;

	sp1.move_bats();
	cout << sp1.get_result() << endl;
	return 0;
}