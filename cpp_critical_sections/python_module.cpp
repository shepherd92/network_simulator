#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "facet_finder.h"

using namespace std;

namespace py = pybind11;

PYBIND11_MODULE(cpp_critical_sections, m)
{
    m.doc() = "pybind11 extract_facets plugin";
    m.def("extract_facets", &extract_facets, "Extract facets from a nested list of simplices.");
    m.def("calc_degree_sequence", &calc_degree_sequence, "Calculate higher order degree sequence from a nested list of simplices.");
}
