#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "connection_generator.h"
#include "numpy_cpp_conversion.h"

namespace py = pybind11;

int main()
{
    const auto seed{0U};

    const auto num_of_nodes{100U};
    const auto num_of_interactions{10U};
    const auto beta{1.0};
    const auto gamma{0.7};
    const auto gamma_prime{0.2};
    const auto torus_size{100.};

    const std::vector<double> model_params_vector{
        num_of_nodes,
        num_of_interactions,
        beta,
        gamma,
        gamma_prime,
        torus_size};

    const auto model_parameters{vector_to_numpy_1d(model_params_vector)};
    const auto connections{generate_finite_network_connections_interface(model_parameters, seed)};
    return 0;
}