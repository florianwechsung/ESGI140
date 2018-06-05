#include "particlemodel.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(pyesgi140, m) {

    m.def("solve_particle_model", &solve_particle_model);
    m.attr("__version__") = "dev";
}
