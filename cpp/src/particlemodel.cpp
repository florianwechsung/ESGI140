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
    public:
        double v_min = 0.1;
        double rho_2 = 2.;
        double rho_1 = .5;
        int rho_bin_shift;
        shared_ptr<vector<double>> rhos;
        vector<double> rhos_sorted;
        vector<double> x_sorted;
        vector<int> x_order;
        vector<double> vmaxs;
        function<double(double)> width_old;
        shared_ptr<vector<double>> rho_binned;
        shared_ptr<vector<double>> velocity;
        vector<double> tstart;
        Barrier width;
        DensitySpeedFunction(vector<double> &_vmaxs, int start_shift, double furthest_distance_possible, Barrier _width, vector<double>& _tstart)
            : rhos(make_shared<vector<double>>(_vmaxs.size(), 0)),
              x_sorted(_vmaxs.size(), 0),
              x_order(_vmaxs.size(), 0),
              vmaxs(_vmaxs),
              rhos_sorted(_vmaxs.size()),
              velocity(make_shared<vector<double>>(_vmaxs.size(), 0)),
              width(_width),
              tstart(_tstart)

        {

            rho_binned = make_shared<vector<double>>(start_shift + int(furthest_distance_possible) + 100, 0.);
            rho_bin_shift = start_shift;
        }
        void operator() (const vector<double> &xs, vector<double> &dxdt, const double t)
        {
            for(int i=0; i<rho_binned->size(); i++)
                (*rho_binned)[i] = 0;
            double xx;
            int xx_floor;
            for(int i=0; i<xs.size(); i++)
            {
                xx = xs[i];
                xx_floor = (int)(xx);
                (*rho_binned)[xx_floor+rho_bin_shift] += (xx_floor + 1 - xx);
                (*rho_binned)[xx_floor+1+rho_bin_shift] += (xx - xx_floor);
            }

            int foresight = 10;
            for(int i=0; i<xs.size(); i++)
            {
                (*rhos)[i] = 0;
                int start = int(xs[i]+1);
                for(int j=start; j<start+foresight; j++)
                {
                    (*rhos)[i] += (*rho_binned)[j + rho_bin_shift]/width(j);
                }
                (*rhos)[i] *= 1./foresight;
            }


            for(int i=0; i<xs.size(); i++)
            {
                if(t < tstart[i])
                    continue;
                dxdt[i] = min(max(vmaxs[i] + ((v_min-vmaxs[i])/(rho_2-rho_1)) * ((*rhos)[i] - rho_1), v_min), vmaxs[i]);
                (*velocity)[i] = dxdt[i];
            }
        }
};

tuple<vector<double>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> solve_particle_model_simple(int num_runners, Barrier width)
{
    double mu = 5.;
    double sd = 0.5;
    auto mersenne_engine = mt19937(1);
    auto normal_dist = normal_distribution<double>(mu, sd);
    auto vmaxs = vector<double>(num_runners);
    generate(begin(vmaxs), end(vmaxs), [&mersenne_engine, &normal_dist](){return normal_dist(mersenne_engine);});
    auto tstarts = vector<double>(num_runners, 0.0);
    return solve_particle_model(vmaxs, tstarts, width);
}

tuple<vector<double>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> solve_particle_model(vector<double>& vmaxs, vector<double>& tstarts, Barrier width)
{
    
    int num_runners = vmaxs.size();
    double start_box_density = 4.;
    double start_box_length = num_runners/(width(0.) * start_box_density);
    double tmax = 40*60;
    //double tmax = 60;
    auto xs = vector<double>(num_runners);

    for(int i=0; i<xs.size(); i++)
        xs[i] = (-i)/(start_box_density * width(0.));



    double furthest_distance_possible = *max_element(begin(vmaxs), end(vmaxs)) * tmax;
    int start_shift = -floor(xs[xs.size()-1]);

    auto rho = vector<double>(num_runners);
    auto rhs_class = DensitySpeedFunction(vmaxs, start_shift, furthest_distance_possible, width, tstarts);


    auto res_x = vector<vector<double>>();
    auto res_rho = vector<vector<double>>();
    auto res_vel = vector<vector<double>>();
    auto res_t = vector<double>();
    function<void(const vector<double>&, double)> observer = [&res_x, &res_t, &res_rho, &rhs_class, &res_vel](const vector<double>& x, double t)
    {
        res_x.push_back(x);
        res_t.push_back(t);
        res_rho.push_back(*(rhs_class.rho_binned));
        res_vel.push_back(*(rhs_class.velocity));
        std::cout << "t=" << t << endl;
    };
    //integrate(rhs_class, xs, 0., tmax, 1., observer); 
    //typedef dense_output_runge_kutta<controlled_runge_kutta<runge_kutta_dopri5<vector<double>>>> stepper_type;
    //integrate_const(stepper_type(), rhs_class, xs, 0.0, tmax, 2.0, observer);
    runge_kutta4< vector<double> > rk4;

    integrate_const(rk4, rhs_class, xs, 0.0, tmax, 2.0, observer);

    return make_tuple(res_t, res_x, res_vel, res_rho);
};
