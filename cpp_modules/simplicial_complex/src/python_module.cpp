#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "degree_sequence.h"
#include "facet_finder.h"
#include "skeleton_extractor.h"

PYBIND11_MODULE(simplicial_complex, m)
{
    m.doc() = "pybind11 critical_sections plugin";
    m.def("calc_degree_sequence", &calc_degree_sequence_interface, "Calculate higher order degree sequence from a nested list of simplices.");
    m.def("create_skeleton", &create_skeleton_interface, "Create skeleton from a list of simplices.");
    m.def("extract_facets", &extract_facets_interface, "Extract facets from a list of simplices.");
}
