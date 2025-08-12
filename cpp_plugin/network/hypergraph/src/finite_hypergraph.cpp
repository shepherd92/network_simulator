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
    if (weighted_)
    {
        assert(SimplexTreeOptions::store_filtration);
    }
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

std::vector<std::vector<std::pair<int32_t, int32_t>>> FiniteHypergraph::calc_persistence_intervals()
{
    // assumption: filtration is decreasing, i.e. the weight of the simplices is decreasing
    // we use this custom function instead of Gudhi since we want to have the weights of the vertices as well
    assert_persistence_cohomology_is_calculated();
    const auto &persistent_pairs{persistent_cohomology_->get_persistent_pairs()};
    std::vector<std::vector<std::pair<int32_t, int32_t>>> intervals{
        static_cast<size_t>(max_dimension_),
        std::vector<std::pair<int32_t, int32_t>>{}};

    for (const auto &persistent_interval : persistent_pairs)
    {
        const Dimension dimension{simplex_tree_->dimension(std::get<0>(persistent_interval))};
        const float birth_weight{simplex_tree_->filtration(std::get<0>(persistent_interval))};
        const float death_weight{simplex_tree_->filtration(std::get<1>(persistent_interval))};
        intervals[dimension].emplace_back(
            birth_weight == 0 ? 0 : static_cast<int32_t>(1. / birth_weight),
            std::isinf(death_weight) ? -1 : static_cast<int32_t>(1. / death_weight));
    }

    return intervals;
}

std::vector<ISimplexList> FiniteHypergraph::calc_persistence_pairs()
{
    assert_persistence_cohomology_is_calculated();
    const auto &persistent_pairs{persistent_cohomology_->get_persistent_pairs()};
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
    if (dimension == 0)
    {
        // return vertices as simplices
        std::vector<Simplex> result{};
        result.reserve(vertices_.size());
        for (const auto vertex_id : vertices_)
        {
            result.emplace_back(Simplex{PointIdList{vertex_id}});
        }
        return SimplexList{result};
    }
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
            const auto degrees{
                weighted_
                    ? calc_simplex_interaction_degree_sequence(dimension)
                    : std::vector<uint32_t>{}};
            simplex_tree_->insert_simplex(
                simplices[i].vertices(),
                weighted_ ? 1. / degrees[i] : 0.);
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

bool FiniteHypergraph::is_weighted() const
{
    return weighted_;
}
