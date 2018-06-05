#include "particlemodel.hpp"
#include <random>
#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;
using namespace std;

template <typename T>
void sort_indexes(const vector<T> &v, vector<int> &I) {

  // initialize original index locations
  iota(I.begin(), I.end(), 0);

  // sort indexes based on comparing values in v
  sort(I.begin(), I.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
}

class DensitySpeedFunction
{
    private:
        double delta_x = 10;
        double width = 7.5;
        double v_min = 1.0;
        double rho_2 = 2.;
        double rho_1 = .5;
        vector<double> rhos;
        vector<double> rhos_sorted;
        vector<double> x_sorted;
        vector<int> x_order;
        vector<double> vmaxs;
    public:
        DensitySpeedFunction(vector<double>& _vmaxs): rhos(_vmaxs.size(), 0), x_sorted(_vmaxs.size(), 0), x_order(_vmaxs.size(), 0), vmaxs(_vmaxs), rhos_sorted(_vmaxs.size()) {}
        void operator() (const vector<double> &xs, vector<double> &dxdt, const double t)
        {
            sort_indexes(xs, x_order);
            for(int i=0; i<xs.size(); i++)
            {
                x_sorted[i] = xs[x_order[i]];
            }
           
            auto it = x_sorted.begin();
            for(int i=0; i<xs.size(); i++)
            {
                auto up = upper_bound(it, x_sorted.end(), x_sorted[i] + delta_x);
                //rhos[x_order[i]] = distance(it, up)/(delta_x * width);
                double act_width = abs(x_sorted[i]-5000) < 10 ? width * 0.1 : width;
                rhos_sorted[i] = (up-it)/(delta_x * act_width);
                it++;
            }
            for(int i=0; i<xs.size(); i++)
            {
                rhos[x_order[i]] = rhos_sorted[i];
            }
            for(int i=0; i<xs.size(); i++)
                dxdt[i] = min(max(vmaxs[i] + ((v_min-vmaxs[i])/(rho_2-rho_1)) * (rhos[i] - rho_1), v_min), vmaxs[i]);
        }
};

pair<vector<double>, vector<vector<double>>> solve_particle_model(int num_runners)
{
    
    //int num_runners = 10;
    double w_0 = 7.5; // width of street
    double start_box_density = 4.;
    double start_box_length = num_runners/(w_0 * start_box_density);
    auto xs = vector<double>(num_runners);

    auto mersenne_engine = mt19937(1);
    //auto unif_dist = uniform_real_distribution<double>(-start_box_length, 0.);

    //generate(begin(xs), end(xs), [&mersenne_engine, &unif_dist](){return unif_dist(mersenne_engine);});
    for(int i=0; i<xs.size(); i++)
    {
        xs[i] = -start_box_density + start_box_density * ((double)i)/xs.size();
    }

    double mu = 5.;
    double sd = 0.5;

    auto normal_dist = normal_distribution<double>(mu, sd);
    auto vmaxs = vector<double>(num_runners);
    generate(begin(vmaxs), end(vmaxs), [&mersenne_engine, &normal_dist](){return normal_dist(mersenne_engine);});

    auto rho = vector<double>(num_runners);
    auto rhs_class = DensitySpeedFunction(vmaxs);

    //function<void(const vector<double>&, vector<double>&, const double)> rhs = [&](const vector<double> &x, vector<double> &dxdt, const double t)
    //{
    //    //simple_speed_function(dxdt, vmaxs, x);
    //    rhs_class.dxdt(dxdt, x, t);
    //};
    double tmax = 30*60;

    auto res_x = vector<vector<double>>();
    auto res_t = vector<double>();
    function<void(const vector<double>&, double)> observer = [&res_x, &res_t](const vector<double>& x, double t)
    {
        res_x.push_back(x);
        res_t.push_back(t);
    };
    integrate(rhs_class, xs, 0., tmax, 1., observer); 
    return make_pair(res_t, res_x);
};