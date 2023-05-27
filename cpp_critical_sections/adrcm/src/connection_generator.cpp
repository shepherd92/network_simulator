#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

#include "connection_generator.h"
#include "model_parameters.h"
#include "numpy_cpp_conversion.h"
#include "point.h"

namespace py = pybind11;

bool is_close(const double first, const double second);
std::vector<Point> create_nodes(
    const py::array_t<double> &birth_times_input,
    const py::array_t<double> &positions_input);
double profile_function(const double argument, const double alpha);

py::array_t<int> generate_connections_default(
    const py::array_t<double> &birth_times_input,
    const py::array_t<double> &positions_input,
    const py::array_t<double> &model_parameters_input)
{
    const auto nodes{create_nodes(birth_times_input, positions_input)};
    const auto num_of_nodes{nodes.size()};
    const auto model_parameters{ModelParameters(model_parameters_input)};

    std::vector<std::pair<int, int>> connections{};
    for (auto source_node_id{0U}; source_node_id < num_of_nodes; ++source_node_id)
    {
        const auto &source_node{nodes[source_node_id]};

        for (auto target_node_id{0U}; target_node_id < source_node_id; ++target_node_id)
        {
            const auto &target_node{nodes[target_node_id]};
            const auto birth_time_ratio{source_node.birth_time() / target_node.birth_time()};
            const auto distance{source_node.torus_distance(target_node, model_parameters.torus_size)};
            const auto profile_function_argument{
                distance * source_node.birth_time() /
                (model_parameters.beta * pow(birth_time_ratio, model_parameters.gamma))};
            if (profile_function_argument < 0.5)
            {
                connections.push_back(std::pair(source_node_id, target_node_id));
            }
        }
    }

    return vector_of_pairs_to_numpy<int>(connections);
}

std::vector<Point> create_nodes(
    const py::array_t<double> &birth_times_input,
    const py::array_t<double> &positions_input)
{
    const auto birth_times{numpy_to_vector_1d<double>(birth_times_input)};
    const auto positions{numpy_to_vector_1d<double>(positions_input)};
    const auto num_of_nodes{birth_times.size()};
    std::vector<Point> nodes;
    nodes.reserve(num_of_nodes);

    for (auto index{0U}; index < num_of_nodes; ++index)
    {
        nodes.emplace_back(Point(birth_times[index], positions[index]));
    }

    return nodes;
}
