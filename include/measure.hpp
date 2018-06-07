#pragma once

#include <vector>
#include <numeric>


using namespace std;

vector<double> calculate_wasted_times(vector<vector<double>>& xs, vector<double>& vmaxs, vector<double>& ts, double race_distance);

double calculate_total_wasted_time(vector<vector<double>>& xs, vector<double>& vmaxs, vector<double>& ts, double race_distance);

double calculate_average_wasted_time(vector<vector<double>>& xs, vector<double>& vmaxs, vector<double>& ts, double race_distance);

vector<double> get_finish_times(vector<vector<double>>& xs, vector<double>& vmaxs, vector<double>& ts, double race_distance);
