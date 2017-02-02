#pragma once
#include <vector>

using std::vector;

class Bat_algorithm {
public:
	typedef vector<double>::iterator double_iterator;
	Bat_algorithm() {};
	Bat_algorithm(int population, int iterations, int dimension, double freq_min, double freq_max,
			      double lower_bound, double upper_bound, double loudness_max = 2);

	virtual double function(const double_iterator begin, const double_iterator end) = 0;
	void move_bats();
	double get_result();
	double get_dimension();
	void set_alpha(const double alpha);

private:
	void bats_init();
	void find_best_bat();
	//check if the current function value is better than global minimum value
	void check_curr_to_best(const double func_new, const vector<double> &curr_positions);
	void cmp_position2bound(vector<double> &bats);

	//update velocities and generate new possible solution
	void adjust_bat_parameters(const int bat_num, vector<double> &parameters);
	double get_average_loudness();
	void gen_local_sol_around_best(const int curr_bat, vector<double> &curr_bat_positions);
	void check_new_solution(const int curr_bat, const int iteration, vector<double> &curr_bat_positions);

	vector<double> freq;
	vector<double> function_values;
	vector<double> best_bat;
	vector<double> bats_loudness;
	vector<double> bats_pulse_rate;

	//2-dimensional arrays: [amount of bats] [amount of dimensions]
	vector<double> bat_velocities;
	vector<double> bats_positions;

	int dimension;
	int population;
	int iterations;
	double freq_min;
	double freq_max;
	double lower_bound;
	double upper_bound;
	double min_value;

	double loudness_max;	
	double alpha;		//[0; 1]  coefficients used to change loudness/pulse_rate
};