#include "measure.hpp"
#include <stdexcept>
#include <string>
#include <cassert>
vector<double> calculate_wasted_times(vector<vector<double>>& xs, vector<double>& vmaxs, vector<double>& ts, double race_distance)
{
    auto num_runners = vmaxs.size();
    auto timesteps = xs.size();
    auto time_startline = vector<double>(num_runners, -1);
    auto time_finishline = vector<double>(num_runners, -1);
    auto time_taken_from_start_to_finish = vector<double>(num_runners, 0.);
    auto wasted_times = vector<double>(num_runners, 0.);

    for(int r=0; r<num_runners;r++)
    {
        for(int i=0; i<timesteps-1;i++)
        {
            if(xs[i][r]<0 && xs[i+1][r] >= 0)
                time_startline[r] = (ts[i]+ts[i+1])/2.; // could calculate a better proxy
            if(xs[i][r] < race_distance && xs[i+1][r] >= race_distance)
            {
                time_finishline[r] = (ts[i]+ts[i+1])/2.; // same issue
                break;
            }
        }
        if(time_startline[r] < 0)
            throw logic_error("Runner " + to_string(r) + " never crossed the start-line.");
        else if(time_finishline[r] > time_startline[r])
            time_taken_from_start_to_finish[r] = time_finishline[r]-time_startline[r];
        else if(time_finishline[r] < 0)
        {
            // estimate how much longer it would take the runner to cross assuming there is no more congestion
            double distance_left = race_distance - xs[timesteps-1][r];
            assert(distance_left>0);
            time_taken_from_start_to_finish[r] = ts[timesteps-1] - time_startline[r] + distance_left/vmaxs[r];
        }
        else
            throw logic_error("Unexpected case");
        wasted_times[r] = time_taken_from_start_to_finish[r] - race_distance/vmaxs[r];
    }
    return wasted_times;
}

double calculate_total_wasted_time(vector<vector<double>>& xs, vector<double>& vmaxs, vector<double>& ts, double race_distance)
{
    auto wasted_times = calculate_wasted_times(xs, vmaxs, ts, race_distance);
    return accumulate(begin(wasted_times), end(wasted_times), 0.);
}

double calculate_average_wasted_time(vector<vector<double>>& xs, vector<double>& vmaxs, vector<double>& ts, double race_distance)
{
    auto total_wasted_time = calculate_total_wasted_time(xs, vmaxs, ts, race_distance);
    return total_wasted_time/vmaxs.size();
}
