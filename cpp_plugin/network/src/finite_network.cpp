#include "finite_network.h"
#include "simplex.h"
#include "simplex_list.h"
#include "tools.h"

FiniteNetwork::FiniteNetwork(
    const Dimension max_dimension,
    const PointIdList &vertices,
    const ISimplexList &interactions,
    const bool weighted)
    : Network{max_dimension, vertices, SimplexList{interactions}},
      weighted_{weighted},
      simplex_tree_{std::nullopt},
      persistent_cohomology_{nullptr}
{
}

FiniteNetwork::FiniteNetwork(
    const Dimension max_dimension,
    const PointIdList &vertices,
    const SimplexList &interactions,
    const bool weighted)
    : Network{max_dimension, vertices, interactions},
      weighted_{weighted},
      simplex_tree_{std::nullopt},
      persistent_cohomology_{nullptr}
{
}

FiniteNetwork::FiniteNetwork(const FiniteNetwork &other)
    : Network{other},
      weighted_{other.weighted()},
      simplex_tree_{std::nullopt},
      persistent_cohomology_{nullptr}
{
}

FiniteNetwork::~FiniteNetwork()
{
    reset_persistence();
}

PointIdList FiniteNetwork::get_vertices() const
{
    return vertices_;
}

SimplexList FiniteNetwork::get_skeleton(const Dimension max_dimension)
{
    SimplexList result{};
    for (auto dimension{0}; dimension <= max_dimension; ++dimension)
    {
        result += get_simplices(dimension);
    }
    return result;
}

SimplexList FiniteNetwork::calc_simplices(const Dimension dimension)
{
    return get_facets().faces(dimension);
}

FiniteNetwork FiniteNetwork::filter(const PointIdList &vertices) const
{
    FiniteNetwork result{*this};
    result.keep_only_vertices(vertices);
    return result;
}

void FiniteNetwork::reset()
{
    Network::reset();
    reset_simplicial_complex();
}

void FiniteNetwork::reset_simplicial_complex()
{
    simplex_tree_ = std::nullopt;
    reset_persistence();
}

void FiniteNetwork::reset_persistence()
{
    delete persistent_cohomology_;
    persistent_cohomology_ = nullptr;
}

void FiniteNetwork::expand()
{
    assert_simplicial_complex_is_built();
    simplex_tree_->expansion(max_dimension_ + 1U);
    reset_persistence();
}

std::vector<ISimplexList> FiniteNetwork::calc_persistence_pairs()
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
            PointIdList birth_simplex{get_simplex_vertices(std::get<0>(persistent_interval))};
            PointIdList death_simplex{get_simplex_vertices(std::get<1>(persistent_interval))};
            std::lock_guard<std::mutex> lock{mutex};
            result.push_back({birth_simplex, death_simplex});
            log_progress(++counter, total, 1000U, "Calc persistence pairs");
        });
    return result;
}

std::vector<std::vector<std::pair<float, float>>> FiniteNetwork::calc_persistence_intervals()
{
    auto &persistent_cohomology{get_persistence()};
    std::vector<std::vector<std::pair<float, float>>> intervals{
        static_cast<size_t>(max_dimension_),
        std::vector<std::pair<float, float>>{}};
    for (Dimension dimension{0}; dimension < max_dimension_; ++dimension)
    {
        intervals[dimension] = persistent_cohomology.intervals_in_dimension(dimension);
    }
    return intervals;
}

FiniteNetwork::PersistentCohomology &FiniteNetwork::get_persistence()
{
    if (!persistent_cohomology_)
    {
        calc_persistent_cohomology();
    }
    return *persistent_cohomology_;
}

void FiniteNetwork::calc_persistent_cohomology()
{
    reset_persistence();
    assert_simplicial_complex_is_built();
    std::cout << "\rCompute persistent cohomology..." << std::flush;
    persistent_cohomology_ = new PersistentCohomology{*simplex_tree_};
    persistent_cohomology_->init_coefficients(2);
    persistent_cohomology_->compute_persistent_cohomology();
    std::cout << "done" << std::flush;
}

std::vector<uint32_t> FiniteNetwork::calc_simplex_interaction_degree_sequence(
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

std::vector<uint32_t> FiniteNetwork::calc_vertex_interaction_degree_distribution() const
{
    // initialize result with zeros
    std::map<PointId, uint32_t> result{};
    for (const auto vertex_id : get_vertices())
    {
        result.emplace(vertex_id, 0U);
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
                    ++result[vertex_id];
                });
        });

    std::vector<uint32_t> counts{};
    counts.reserve(result.size());
    for (std::map<PointId, uint32_t>::iterator it = result.begin(); it != result.end(); ++it)
    {
        counts.push_back(it->second);
    }

    return counts;
}

std::vector<uint32_t> FiniteNetwork::calc_coface_degree_sequence(
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension)
{
    assert(neighbor_dimension > simplex_dimension);

    const auto &cofaces{get_simplices(neighbor_dimension)};
    auto simplex_degree_map{cofaces.calc_degree_sequence(simplex_dimension)};

    // order of the degree values does not matter
    std::vector<uint32_t> result{};
    const auto &simplices{get_simplices(simplex_dimension)};
    result.reserve(simplex_degree_map.size());
    for (const auto &simplex : simplices)
    {
        result.emplace_back(simplex_degree_map[simplex]);
    }

    return result;
}

std::vector<int32_t> FiniteNetwork::calc_betti_numbers()
{
    std::vector<int32_t> result{max_dimension_, 0};
    if (num_simplices(0) == 1U)
    {
        // handle an error in Gudhi
        result[0] = 1;
    }
    else
    {
        const auto &persistent_cohomology{get_persistence()};
        result = persistent_cohomology.betti_numbers();
        result.resize(max_dimension_);
    }
    assert(static_cast<int32_t>(result.size()) == max_dimension_);
    return result;
}

void FiniteNetwork::assert_simplicial_complex_is_built()
{
    if (!is_valid())
    {
        create_simplicial_complex();
    }
}

void FiniteNetwork::assert_simplicial_complex_is_initialized()
{
    if (!is_valid())
    {
        simplex_tree_ = SimplexTree{};
    }
}

bool FiniteNetwork::is_valid() const
{
    return simplex_tree_.has_value();
}

void FiniteNetwork::create_simplicial_complex()
{
    reset_simplicial_complex();
    add_vertices(vertices_);
    fill_simplicial_complex();
}

void FiniteNetwork::fill_simplicial_complex()
{
    assert_simplicial_complex_is_initialized();

    std::cout << "\rInsert simplices" << std::flush;
    std::mutex mutex{};

    for (Dimension dimension{0}; dimension <= max_dimension_; ++dimension)
    {
        const auto &simplices{get_simplices(dimension)};
        std::atomic<uint32_t> counter{0U};
        std::vector<uint32_t> weights{};
        weights.reserve(simplices.size());
        if (weighted())
        {
            auto simplex_weight_map{interactions_.calc_degree_sequence(dimension)};

            for (const auto &simplex : simplices)
            {
                weights.emplace_back(simplex_weight_map[simplex]);
            }

            const auto max_weight{*std::max_element(weights.begin(), weights.end())};
            for (auto &weight : weights)
            {
                weight = max_weight - weight;
            }
        }
        else
        {
            weights = std::vector<uint32_t>{simplices.size(), 0U};
        }

        const auto total{simplices.size()};
        for (size_t i = 0; i < simplices.size(); ++i)
        {
            simplex_tree_->insert_simplex_and_subfaces(simplices[i].vertices(), weights[i]);
            log_progress(counter, total, 1000U, "Insert simplices");
        }
    }

    reset_persistence();
}

void FiniteNetwork::add_vertices(const PointIdList &vertices)
{
    assert_simplicial_complex_is_initialized();
    std::for_each(
        std::execution::seq,
        vertices.begin(),
        vertices.end(),
        [this](const auto vertex)
        { simplex_tree_->insert_simplex_and_subfaces(PointIdList{vertex}); });
}

PointIdList FiniteNetwork::get_simplex_vertices(const SimplexHandle &simplex_handle)
{
    assert_simplicial_complex_is_built();
    PointIdList result{};
    if (simplex_handle != simplex_tree_->null_simplex())
    {
        for (const auto &vertex : simplex_tree_->simplex_vertex_range(simplex_handle))
        {
            result.push_back(vertex);
        }
    }
    return result;
}

bool FiniteNetwork::weighted() const
{
    return weighted_;
}