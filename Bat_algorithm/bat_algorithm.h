#pragma once
#include <vector>

using std::vector;

class Bat_algorithm {
public:
	typedef vector<double>::iterator double_iterator;
	Bat_algorithm() {};
	Bat_algorithm(int population, int iterations, int dimension, double freq_min, double freq_max,
			      double lower_bound, double upper_bound, double loudness_max, double loudness_min,
				  double pulse_rate_max, double pulse_rate_min);

	virtual double function(const double_iterator begin, const double_iterator end) = 0;
	void move_bats();
	double get_result();
	double get_dimension();

private:
	int Bat_algorithm::select_rand_bat(const int curr_bat);

	void bats_init();
	void find_best_bat();
	//check if the current function value is better than global minimum value
	void check_curr_to_best(const double func_new, const vector<double> &curr_positions);
	void cmp_position2bound(const double_iterator begin, const double_iterator end);

	void local_search(const int bat_num, vector<double> &parameters);
	//update velocities and generate new possible solution
	void update_bat_position(const int bat_num, vector<double> &parameters);
	double get_average_loudness();
	void gen_local_sol_around_best(const int curr_bat, vector<double> &curr_bat_positions);
	void check_new_solution(const int curr_bat, const int iteration, vector<double> &curr_bat_positions);

	vector<double> bat_freq;
	vector<double> fitness_func;
	vector<double> bat_best;
	vector<double> bat_loudness;
	vector<double> bat_pulse_rate;

	//2-dimensional array: [amount of bats] [amount of dimensions]
	vector<double> bat_position;

	int dimension;
	int population;
	int iterations;
	double freq_min;
	double freq_max;
	double lower_bound;
	double upper_bound;
	double min_value;

	double loudness_max;
	double loudness_min;
	double pulse_rate_max;
	double pulse_rate_min;

	const double lambda_min = 0.3;
	const double lambda_max = 1.99;
	double mantegna_scale;
};