#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "numpy_cpp_conversion.h"

namespace py = pybind11;

int main()
{
    // const auto seed{0U};

    const auto max_dimension{3U};
    const auto network_size{10000U};
    const auto num_of_nodes{10000U};
    const auto alpha{0.5};
    const auto beta{1.0};
    const auto gamma{0.7};
    const auto torus_dimension{1U};
    const auto torus_size_in_1_dimension{0.};

    const std::vector<double> model_params_vector{
        max_dimension,
        network_size,
        num_of_nodes,
        alpha,
        beta,
        gamma,
        torus_dimension,
        torus_size_in_1_dimension};

    const auto model_parameters{to_numpy(model_params_vector)};
    return 0;
}