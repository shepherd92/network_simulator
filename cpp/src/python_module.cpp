#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>

#include "typedefs.h"
#include "python_module.h"

namespace py = pybind11;

NetworkInterface::NetworkInterface(const dimension max_dimension)
{
    network.max_dimension(max_dimension);
}

void NetworkInterface::add_simplices(const py::list &simplices)
{
    const auto simplices_vector{simplices.cast<std::vector<std::vector<vertex_id>>>()};
    network.add_simplices(simplices_vector);
}

py::array_t<vertex_id> NetworkInterface::get_simplices_by_dimension(const dimension dimension) const
{
    const auto simplices{network.get_simplices_by_dimension(dimension)};

    std::vector<std::vector<vertex_id>> result{};
    for (const auto &simplex : simplices)
    {
        result.push_back(std::vector<vertex_id>(simplex.begin(), simplex.end()));
    }

    return py::array(py::cast(std::move(result)));
}

py::list NetworkInterface::facets()
{
    std::vector<std::vector<vertex_id>> result{};
    for (const auto &facet : network.facets())
    {
        result.push_back(std::vector<vertex_id>(facet.begin(), facet.end()));
    }

    return py::list(py::cast(std::move(result)));
}

PYBIND11_MODULE(cpp_library, m)
{
    m.doc() = "pybind11 C++ modules";

    py::class_<NetworkInterface>(m, "Network")
        .def(py::init<const dimension>())
        .def("add_simplices", &NetworkInterface::add_simplices)
        .def("facets", &NetworkInterface::facets)
        .def("get_simplices_by_dimension", &NetworkInterface::get_simplices_by_dimension);
}
