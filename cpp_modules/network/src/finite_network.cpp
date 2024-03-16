#include "finite_network.h"

FiniteNetwork::FiniteNetwork(const Dimension max_dimension, const PointIdList &vertices, const ISimplexList &interactions)
    : Network{max_dimension, vertices, create_simplices(interactions)},
      simplex_tree_{std::nullopt},
      persistent_cohomology_{nullptr}
{
}

FiniteNetwork::FiniteNetwork(const Dimension max_dimension, const PointIdList &vertices, const SimplexList &interactions)
    : Network{max_dimension, vertices, interactions},
      simplex_tree_{std::nullopt},
      persistent_cohomology_{nullptr}
{
}

FiniteNetwork::FiniteNetwork(const FiniteNetwork &other)
    : Network{other},
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

SimplexList FiniteNetwork::get_simplices(const Dimension dimension)
{
    return get_faces_simplices(get_facets(), dimension);
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

void FiniteNetwork::calc_persistent_cohomology()
{
    reset_persistence();
    assert_simplicial_complex_is_built();
    persistent_cohomology_ = new PersistentCohomology{*simplex_tree_};
    persistent_cohomology_->init_coefficients(2);
    persistent_cohomology_->compute_persistent_cohomology();
}

const FiniteNetwork::PersistentCohomology &FiniteNetwork::get_persistence()
{
    if (!persistent_cohomology_)
    {
        calc_persistent_cohomology();
    }
    return *persistent_cohomology_;
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
            if (++counter % 1000 == 0 && total > 10000U)
            {
                std::cout << "\rCalc persistence pairs ... " << counter << " / " << total;
            }
            result.push_back({birth_simplex, death_simplex});
        });

    if (total > 10000U)
    {
        std::cout << "\rCalc persistence pairs ... " << total << " / " << total;
    }

    return result;
}

std::vector<int32_t> FiniteNetwork::calc_betti_numbers()
{
    std::vector<int32_t> result{max_dimension_, 0};
    if (num_simplices() == 1U)
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
    add_simplices(interactions_);
}

void FiniteNetwork::add_simplices(const SimplexList &simplices)
{
    assert_simplicial_complex_is_initialized();
    const auto representable_simplices{get_skeleton_simplices(simplices, max_dimension_)};
    const auto total{representable_simplices.size()};
    std::atomic<uint32_t> counter{0U};

    std::for_each(
        std::execution::seq,
        representable_simplices.begin(),
        representable_simplices.end(),
        [&](const auto &simplex)
        {
            if (++counter % 1000 == 0 && total > 10000U)
            {
                std::cout << "\rInsert simplices ... " << counter << " / " << total;
            }
            simplex_tree_->insert_simplex_and_subfaces(simplex.vertices());
        });
    reset_persistence();

    if (total > 10000U)
    {
        std::cout << "\rInsert simplices ... " << total << " / " << total;
    }
}

void FiniteNetwork::add_vertices(const PointIdList &vertices)
{
    SimplexList simplices{};
    simplices.reserve(vertices.size());

    std::transform(
        vertices.begin(),
        vertices.end(),
        std::back_inserter(simplices),
        [](const auto vertex)
        { return Simplex{PointIdList{vertex}}; });

    add_simplices(simplices);
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
