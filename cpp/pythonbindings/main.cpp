#include "particlemodel.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <memory>

namespace py = pybind11;

PYBIND11_MODULE(pyesgi140, m) {

    m.def("solve_particle_model", &solve_particle_model);
    m.def("solve_particle_model_simple", &solve_particle_model);

    py::class_<Barrier,
               std::shared_ptr<Barrier>>(
        m, "Barrier")
        .def(py::init<double, double, double, double, double>())
        .def("__call__", &Barrier::operator());

    m.attr("__version__") = "dev";
}
