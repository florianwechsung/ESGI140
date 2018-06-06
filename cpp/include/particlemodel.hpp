#pragma once
#include <functional>
#include <vector>

std::pair<std::vector<double>, std::vector<std::vector<double>>> solve_particle_model(int, std::function<double(double)> width);
