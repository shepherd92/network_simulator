#include "finite_hypergraph.hpp"
#include "simplex.hpp"
#include "simplex_list.hpp"
#include "tools.hpp"

FiniteHypergraph::FiniteHypergraph(
    const Dimension max_dimension,
    const PointIdList &vertices,
    const ISimplexList &interactions,
    const bool weighted)
    : Network{max_dimension, vertices},
      FiniteNetwork{},
      Hypergraph{SimplexList{interactions}},
      weighted_{weighted}
{
}

FiniteHypergraph::FiniteHypergraph(
    const Dimension max_dimension,
    const PointIdList &vertices,
    const SimplexList &interactions,
    const bool weighted)
    : Network{max_dimension, vertices},
      FiniteNetwork(),
      Hypergraph{interactions},
      weighted_{weighted}
{
}

FiniteHypergraph::FiniteHypergraph(FiniteHypergraph &&other) noexcept
    : Network{std::move(other)},
      FiniteNetwork{std::move(other)},
      Hypergraph{std::move(other)},
      weighted_{std::move(other.weighted_)}
{
}

// move assignment defined due to virtual base class
FiniteHypergraph &FiniteHypergraph::operator=(FiniteHypergraph &&other) noexcept
{
    if (this != &other)
    {
        Network::operator=(std::move(other));
        Hypergraph::operator=(std::move(other));
        FiniteNetwork::operator=(std::move(other));
        weighted_ = std::move(other.weighted_);
    }
    return *this;
}

std::vector<uint32_t> FiniteHypergraph::calc_vertex_interaction_degree_distribution() const
{
    // initialize result with zeros
    std::map<PointId, uint32_t> point_id_interaction_count_map{};
    for (const auto vertex_id : vertices_)
    {
        point_id_interaction_count_map.emplace(vertex_id, 0U);
    }

    std::for_each(
        std::execution::seq,
        interactions_.simplices().begin(),
        interactions_.simplices().end(),
        [&](const auto &interaction)
        {
            std::for_each(
                execution_policy, // can execute parallel, different vertices in an interaction
                interaction.vertices().begin(),
                interaction.vertices().end(),
                [&](const auto vertex_id)
                {
                    ++point_id_interaction_count_map[vertex_id];
                });
        });

    std::vector<uint32_t> counts{};
    counts.reserve(vertices_.size());
    for (const auto vertex_id : vertices_)
    {
        counts.emplace_back(point_id_interaction_count_map[vertex_id]);
    }

    return counts;
}

std::vector<std::vector<std::pair<float, float>>> FiniteHypergraph::calc_persistence_intervals()
{
    // assumption: filtration is decreasing, i.e. the weight of the simplices is decreasing
    // we use this custom function instead of Gudhi since we want to have the weights of the vertices as well
    auto &persistent_cohomology{get_persistence()};
    const auto &persistent_pairs{persistent_cohomology.get_persistent_pairs()};
    std::vector<std::vector<std::pair<float, float>>> intervals{
        static_cast<size_t>(max_dimension_),
        std::vector<std::pair<float, float>>{}};

    std::mutex mutex{};
    const auto total{persistent_pairs.size()};
    std::atomic<uint32_t> counter{0U};

    std::unordered_map<Simplex, uint32_t, SimplexHash> simplex_interaction_map{};

    for (Dimension dimension{0}; dimension <= max_dimension_; ++dimension)
    {
        const auto interaction_degree_map_for_dimension{interactions_.calc_degree_sequence(dimension)};
        simplex_interaction_map.insert(
            interaction_degree_map_for_dimension.begin(),
            interaction_degree_map_for_dimension.end());
    }

    std::for_each(
        execution_policy,
        persistent_pairs.begin(),
        persistent_pairs.end(),
        [&](const auto &persistent_interval)
        {
            const Simplex birth_simplex{get_simplex_vertices_simplex_tree(std::get<0>(persistent_interval))};
            const Simplex death_simplex{get_simplex_vertices_simplex_tree(std::get<1>(persistent_interval))};
            assert(simplex_interaction_map.contains(birth_simplex.vertices()));
            const auto birth_weight{simplex_interaction_map[birth_simplex]};
            const auto death_weight{
                simplex_interaction_map.contains(death_simplex.vertices())
                    ? simplex_interaction_map[death_simplex.vertices()]
                    : 0U};
            log_progress(++counter, total, 1000U, "Calc persistence pairs");
            if (birth_weight == death_weight)
            {
                return;
            }

            std::lock_guard<std::mutex> lock{mutex};
            intervals[birth_simplex.dimension()].push_back({birth_weight, death_weight});
        });

    return intervals;
}

std::vector<ISimplexList> FiniteHypergraph::calc_persistence_pairs()
{
    const auto &persistent_cohomology{get_persistence()};
    const auto &persistent_pairs{persistent_cohomology.get_persistent_pairs()};
    std::vector<ISimplexList> result{};
    result.reserve(persistent_pairs.size());

    std::mutex mutex{};
    const auto total{persistent_pairs.size()};
    std::atomic<uint32_t> counter{0U};

    std::for_each(
        execution_policy,
        persistent_pairs.begin(),
        persistent_pairs.end(),
        [&](const auto &persistent_interval)
        {
            const PointIdList birth_simplex{get_simplex_vertices_simplex_tree(std::get<0>(persistent_interval))};
            const PointIdList death_simplex{get_simplex_vertices_simplex_tree(std::get<1>(persistent_interval))};
            std::lock_guard<std::mutex> lock{mutex};
            result.push_back({birth_simplex, death_simplex});
            log_progress(++counter, total, 1000U, "Calc persistence pairs");
        });
    return result;
}

SimplexList FiniteHypergraph::calc_simplices(const Dimension dimension)
{
    return interactions_.faces(dimension);
}

void FiniteHypergraph::fill_simplex_tree()
{
    assert_simplex_tree_is_initialized();

    std::cout << "\rInsert simplices" << std::flush;

    // insert simplices
    for (Dimension dimension{0}; dimension <= max_dimension_; ++dimension)
    {
        const auto &simplices{get_simplices(dimension)};
        const auto total{simplices.size()};
        std::atomic<uint32_t> counter{0U};
        for (auto i{0U}; i < simplices.size(); ++i)
        {
            simplex_tree_->insert_simplex_and_subfaces(
                simplices[i].vertices(), 0.);
            log_progress(++counter, total, 1000U, "Insert simplices");
        }
    }
    reset_persistence();
}

std::vector<uint32_t> FiniteHypergraph::calc_simplex_interaction_degree_sequence(
    const Dimension simplex_dimension)
{
    if (simplex_dimension == 0)
    {
        return calc_vertex_interaction_degree_distribution();
    }

    auto simplex_degree_map{interactions_.calc_degree_sequence(simplex_dimension)};

    // order of the degree values does not matter
    std::vector<uint32_t> result{};
    const auto &simplices{get_simplices(simplex_dimension)};
    result.reserve(simplex_degree_map.size());
    for (const auto &simplex : simplices)
    {
        // simplex has typical vertex, simplex_degree_map is for simplices without typical vertex
        result.emplace_back(simplex_degree_map[simplex]);
    }

    return result;
}

FiniteHypergraph FiniteHypergraph::filter(const PointIdList &vertices)
{
    SimplexList filtered_interactions{interactions_.filter(vertices)};
    return FiniteHypergraph{max_dimension_, vertices, std::move(filtered_interactions), weighted_};
}