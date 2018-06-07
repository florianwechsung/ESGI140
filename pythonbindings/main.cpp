#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <memory>
#include "measure.hpp"
#include "particlemodel.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pyesgi140, m) {
    m.def("solve_particle_model", &solve_particle_model, py::arg("vmaxs"),
          py::arg("tstarts"), py::arg("width"), py::arg("slope"), py::arg("speedSettings"), py::arg("tmax") = 3600.,
          py::arg("dt") = 1.0, py::arg("rho_start") = 4.);
    m.def("solve_particle_model_simple", &solve_particle_model);
    m.def("calculate_total_wasted_time", &calculate_total_wasted_time);
    m.def("calculate_average_wasted_time", &calculate_average_wasted_time);
    m.def("get_finish_times", &get_finish_times);

    py::class_<Barrier, std::shared_ptr<Barrier>>(m, "Barrier")
        .def(py::init<double, double, double, double, double>())
        .def("__call__", &Barrier::operator());

    py::class_<SpeedSettings, std::shared_ptr<SpeedSettings>>(m, "SpeedSettings")
        .def(py::init<double, double, double>(), py::arg("v_min")=0.1, py::arg("rho_1")=0.5, py::arg("rho_2")=2.);


    m.attr("__version__") = "dev";
}
