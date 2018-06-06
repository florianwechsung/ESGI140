#include "particlemodel.hpp"

int main()
{
    solve_particle_model(200, [](double x){return 10.;});
}
