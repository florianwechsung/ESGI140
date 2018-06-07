#include "particlemodel.hpp"

int main()
{
    auto barrier = Barrier(4900, 30, 10, 7.5, 5.0);
    solve_particle_model_simple(1000, barrier);
}
