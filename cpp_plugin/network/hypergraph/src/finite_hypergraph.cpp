#include "finite_hypergraph.h"
#include "simplex.h"
#include "simplex_list.h"
#include "tools.h"

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

std::vector<uint32_t> FiniteHypergraph::calc_vertex_interaction_degree_distribution() const
{
    // initialize result with zeros
    std::map<PointId, uint32_t> point_id_interaction_count_map{};
    for (const auto vertex_id : get_vertices())
    {
        point_id_interaction_count_map.emplace(vertex_id, 0U);
    }

    std::for_each(
        std::execution::seq,
        get_interactions().simplices().begin(),
        get_interactions().simplices().end(),
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
    counts.reserve(get_vertices().size());
    for (const auto vertex_id : get_vertices())
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
        const auto interaction_degree_map_for_dimension{get_interactions().calc_degree_sequence(dimension)};
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
            const Simplex birth_simplex{get_simplex_vertices(std::get<0>(persistent_interval))};
            const Simplex death_simplex{get_simplex_vertices(std::get<1>(persistent_interval))};
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
            const PointIdList birth_simplex{get_simplex_vertices(std::get<0>(persistent_interval))};
            const PointIdList death_simplex{get_simplex_vertices(std::get<1>(persistent_interval))};
            std::lock_guard<std::mutex> lock{mutex};
            result.push_back({birth_simplex, death_simplex});
            log_progress(++counter, total, 1000U, "Calc persistence pairs");
        });
    return result;
}

SimplexList FiniteHypergraph::calc_simplices(const Dimension dimension)
{
    return get_interactions().faces(dimension);
}

void FiniteHypergraph::fill_simplicial_complex()
{
    assert_simplicial_complex_is_initialized();

    std::cout << "\rInsert simplices" << std::flush;

    // insert simplices
    for (Dimension dimension{0}; dimension <= max_dimension_; ++dimension)
    {
        const auto &simplices{get_simplices(dimension)};
        const auto total{simplices.size()};
        std::atomic<uint32_t> counter{0U};
        for (auto i{0U}; i < simplices.size(); ++i)
        {
            get_simplex_tree()->insert_simplex_and_subfaces(
                simplices[i].vertices(), 0.);
            log_progress(++counter, total, 1000U, "Insert simplices");
        }
    }
    reset_persistence();
}

bool FiniteHypergraph::weighted() const
{
    return weighted_;
}