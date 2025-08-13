#include "finite_adrcm_model.hpp"
#include "finite_clique_complex.hpp"
#include "point.hpp"
#include "simplex_list.hpp"
#include "typedefs.hpp"

FiniteAdrcmModel::FiniteAdrcmModel(const std::vector<double> &parameters_in, const uint32_t seed)
    : Model{seed}, FiniteModel{seed}, AdrcmModel{parameters_in}
{
}

std::tuple<FiniteCliqueComplex, MarkPositionList> FiniteAdrcmModel::generate_network() const
{
    const auto num_of_vertices{std::poisson_distribution<uint32_t>(
        lambda() * torus_size)(random_number_generator_)};
    const auto vertices{create_points(num_of_vertices)};
    const auto connections{generate_connections(vertices)};

    const auto vertex_ids{convert_to_id_list(vertices)};
    FiniteCliqueComplex network{max_dimension(), vertex_ids, connections};
    auto mark_position_pairs{convert_to_mark_position_pairs(vertices)};

    return std::tuple<FiniteCliqueComplex, MarkPositionList>(
        std::move(network),
        std::move(mark_position_pairs));
}
