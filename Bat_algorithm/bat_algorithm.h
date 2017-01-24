#pragma once
#include <vector>

using namespace std;

class Bat_algorithm {
public:
	typedef vector<double>::iterator double_iterator;

	Bat_algorithm() {};

	Bat_algorithm(int population, int iterations, int dimension, double freq_min, double freq_max, 
		double lower_bound, double upper_bound, double loudness = 0.9, double pulse_rate = 0.9) {
		this->population = population;
		this->iterations = iterations;
		this->dimension = dimension;
		this->freq_min = freq_min;
		this->freq_max = freq_max;
		this->lower_bound = lower_bound;
		this->upper_bound = upper_bound;

		loudness = 0.9;
		pulse_rate = 0.9;
		alpha = 0.9;
	}

	virtual double function(const double_iterator begin, const double_iterator end) = 0;
	void move_bats();
	
	double get_result() {return min_value;}

protected:
	void bats_init();
	void find_best_bat();
	void cmp_position2bound(vector<double> &bats);

	vector<double> freq;
	vector<double> function_values;
	vector<double> best_bat_values;

	//2-dimensional arrays: [amount of bats] [amount of dimensions]
	vector<double> bat_velocities;
	vector<double> bat_positions;

	int dimension;
	int population;
	int iterations;
	double freq_min;
	double freq_max;
	double lower_bound;
	double upper_bound;

	double min_value;

	double loudness;	//[0; 1]
	double pulse_rate;	//[0; 1]
	double alpha;		//[0; 1]  coefficients used to change loudness/pulse_rate
};