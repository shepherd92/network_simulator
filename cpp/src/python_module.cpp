#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "adrcm_model.h"
#include "network.h"

namespace py = pybind11;

PYBIND11_MODULE(cpp, m)
{
    m.doc() = "pybind11 C++ modules";
    py::class_<AdrcmModel>(m, "AdrcmModel")
        .def(py::init<const std::string &>())
        .def("get_facets", &AdrcmModel::create_finite_network);
    py::class_<Network>(m, "Network")
        .def(py::init<const std::string &>())
        .def("get_facets", &Network::get_facets);
}
