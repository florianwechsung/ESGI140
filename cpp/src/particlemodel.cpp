#include "particlemodel.hpp"
#include <algorithm>
#include <boost/numeric/odeint.hpp>
#include <functional>
#include <iterator>
#include <numeric>
#include <random>

using namespace boost::numeric::odeint;
using namespace std;



class DensitySpeedFunction {
   public:
    int rho_bin_shift;
    shared_ptr<vector<double>> rhos;
    shared_ptr<vector<double>> rho_binned;
    shared_ptr<vector<double>> velocity;
    vector<double> vmaxs;
    vector<double> width_binned;
    vector<double> tstart;
    function<double(double)> width;
    SpeedSettings settings;
    DensitySpeedFunction(vector<double> &_vmaxs, vector<double>& xs,
                         vector<double> &_tstart, function<double(double)> _width, double tmax, SpeedSettings& speedSettings)
        : rhos(make_shared<vector<double>>(_vmaxs.size(), 0)),
          vmaxs(_vmaxs),
          velocity(make_shared<vector<double>>(_vmaxs.size(), 0)),
          width(_width),
          tstart(_tstart),
          settings(speedSettings)

    {
        double furthest_distance_possible =
            *max_element(begin(vmaxs), end(vmaxs)) * tmax;
        int start_shift = -floor(xs[xs.size() - 1]);
        rho_binned = make_shared<vector<double>>(
            start_shift + int(furthest_distance_possible) + 100, 0.);
        rho_bin_shift = start_shift;
        width_binned = vector<double>(rho_binned->size(), 0.);
        for(int i=0; i<rho_binned->size(); i++)
            width_binned[i] = width(i - rho_bin_shift);
    }
    void operator()(const vector<double> &xs, vector<double> &dxdt,
                    const double t) {

        for (int i = 0; i < rho_binned->size(); i++) (*rho_binned)[i] = 0;
        double xx;
        int xx_floor;
        for (int i = 0; i < xs.size(); i++) {
            xx = xs[i];
            xx_floor = (int)(xx);
            (*rho_binned)[xx_floor + rho_bin_shift] += (xx_floor + 1 - xx);
            (*rho_binned)[xx_floor + 1 + rho_bin_shift] += (xx - xx_floor);
        }


        int foresight = 10;
        for (int i = 0; i < xs.size(); i++) {
            (*rhos)[i] = 0;
            int start = int(xs[i] + 1);
            for (int j = start; j < start + foresight; j++) {
                (*rhos)[i] += (*rho_binned)[j + rho_bin_shift] / width_binned[j + rho_bin_shift];
            }
            (*rhos)[i] *= 1. / foresight;
        }
        auto rho_1 = settings.rho_1;
        auto rho_2 = settings.rho_2;
        auto v_min = settings.v_min;
        for (int i = 0; i < xs.size(); i++) {
            if (t < tstart[i]) continue;
            dxdt[i] =
                min(max(vmaxs[i] + ((v_min - vmaxs[i]) / (rho_2 - rho_1)) *
                                       ((*rhos)[i] - rho_1),
                        v_min),
                    vmaxs[i]);
            (*velocity)[i] = dxdt[i];
        }
    }
};

tuple<vector<double>, vector<vector<double>>, vector<vector<double>>,
      vector<vector<double>>>
solve_particle_model_simple(int num_runners, function<double(double)> width) {
    double mu = 5.;
    double sd = 0.5;
    auto mersenne_engine = mt19937(1);
    auto normal_dist = normal_distribution<double>(mu, sd);
    auto vmaxs = vector<double>(num_runners);
    generate(begin(vmaxs), end(vmaxs), [&mersenne_engine, &normal_dist]() {
        return normal_dist(mersenne_engine);
    });
    auto tstarts = vector<double>(num_runners, 0.0);
    auto settings = SpeedSettings(0.1, 0.5, 2.0);
    return solve_particle_model(vmaxs, tstarts, width, settings);
}

tuple<vector<double>, vector<vector<double>>, vector<vector<double>>,
      vector<vector<double>>>
solve_particle_model(vector<double> &vmaxs, vector<double> &tstarts,
                     function<double(double)>& width, SpeedSettings& speedSettings, double tmax, double dt, double rho_start) {
    int num_runners = vmaxs.size();
    double start_box_length = num_runners / (width(0.) * rho_start);
    auto xs = vector<double>(num_runners);

    for (int i = 0; i < xs.size(); i++)
        xs[i] = (-(i + 1)) / (rho_start * width(0.));

    auto rho = vector<double>(num_runners);
    auto rhs_class = DensitySpeedFunction(
        vmaxs, xs, tstarts, width, tmax, speedSettings);

    auto res_x = vector<vector<double>>();
    auto res_rho = vector<vector<double>>();
    auto res_vel = vector<vector<double>>();
    auto res_t = vector<double>();
    double last_printed = -1000;
    function<void(const vector<double> &, double)> observer =
        [&res_x, &res_t, &res_rho, &rhs_class, &res_vel, &last_printed](
            const vector<double> &x, double t) {
            res_x.push_back(x);
            res_t.push_back(t);
            res_rho.push_back(*(rhs_class.rho_binned));
            res_vel.push_back(*(rhs_class.velocity));
            //if (t >= last_printed + 60 - 1e-13) {
            //    cout << "t=" << t << endl;
            //    last_printed = t;
            //}
        };
    runge_kutta4<vector<double>> rk4;
    integrate_const(rk4, rhs_class, xs, 0.0, tmax, dt, observer);

    return make_tuple(res_t, res_x, res_vel, res_rho);
};
