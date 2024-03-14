#include "finite_adrcm_model.h"
#include "finite_network.h"
#include "point.h"
#include "typedefs.h"

FiniteAdrcmModel::FiniteAdrcmModel(const py::array_t<double> &parameters_in, const uint32_t seed)
    : Model{seed}, FiniteModel{seed}, AdrcmModel{parameters_in}
{
}

std::tuple<FiniteNetwork, MarkPositionList> FiniteAdrcmModel::generate_network() const
{
    const auto num_of_vertices{std::poisson_distribution<uint32_t>(lambda() * torus_size())(random_number_generator_)};
    const auto vertices{create_points(num_of_vertices)};
    const auto connections{generate_connections(vertices)};

    const auto vertex_ids{convert_to_id_list(vertices)};
    const auto simplices{create_simplices_from_connections(connections)};
    const FiniteNetwork network{max_dimension(), vertex_ids, simplices};
    const auto mark_position_pairs{convert_to_mark_position_pairs(vertices)};

    return std::make_tuple<>(network, mark_position_pairs);
}
