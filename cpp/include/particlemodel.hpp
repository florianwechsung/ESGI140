#pragma once
#include <functional>
#include <vector>

std::tuple<std::vector<double>, std::vector<std::vector<double>>, std::vector<std::vector<double>>> solve_particle_model(int, std::function<double(double)> width);
