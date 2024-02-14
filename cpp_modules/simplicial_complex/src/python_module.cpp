#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "degree_sequence.h"
#include "facet_finder.h"

PYBIND11_MODULE(simplicial_complex, m)
{
    m.doc() = "pybind11 critical_sections plugin";
    m.def("extract_facets", &extract_facets, "Extract facets from a nested list of simplices.");
    m.def("calc_degree_sequence", &calc_degree_sequence, "Calculate higher order degree sequence from a nested list of simplices.");
}
