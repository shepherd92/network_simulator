#include "finite_network.h"

FiniteNetwork::FiniteNetwork(const VertexList &vertices, const ISimplexList &interactions)
    : Network{vertices, interactions},
      persistent_cohomology_{nullptr}
{
}

FiniteNetwork::~FiniteNetwork()
{
    reset_persistence();
}

void FiniteNetwork::add_vertices(const VertexList &vertices)
{
    Network::add_vertices(vertices);
    reset_persistence();
}

void FiniteNetwork::reset_simplicial_complex()
{
    Network::reset_simplicial_complex();
    reset_persistence();
}

void FiniteNetwork::reset_persistence()
{
    delete persistent_cohomology_;
    persistent_cohomology_ = nullptr;
}

void FiniteNetwork::expand(const Dimension max_dimension)
{
    Network::expand(max_dimension);
    reset_persistence();
}

void FiniteNetwork::calc_persistent_cohomology()
{
    assert(is_valid());
    reset_persistence();
    persistent_cohomology_ = new PersistentCohomology{*simplex_tree_};
    persistent_cohomology_->init_coefficients(2);
    persistent_cohomology_->compute_persistent_cohomology();
}

void FiniteNetwork::add_simplices(const SimplexList &simplices, const Dimension dimension)
{
    Network::add_simplices(simplices, dimension);
    reset_persistence();
}

Network::SimplexHandleList FiniteNetwork::get_simplices()
{
    assert(is_valid());
    return simplex_tree_->filtration_simplex_range();
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
            VertexList birth_simplex{get_vertices(std::get<0>(persistent_interval))};
            VertexList death_simplex{get_vertices(std::get<1>(persistent_interval))};
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
    const auto &persistent_cohomology{get_persistence()};
    const auto result{persistent_cohomology.betti_numbers()};
    return result;
}

uint32_t FiniteNetwork::num_simplices()
{
    assert(is_valid());
    return simplex_tree_->num_simplices();
}

SimplexList FiniteNetwork::get_skeleton_interactions(const Dimension max_dimension)
{
    return get_skeleton_simplices(interactions_, max_dimension);
}

SimplexList FiniteNetwork::get_skeleton_simplicial_complex(const Dimension max_dimension)
{
    const auto &simplices{simplex_tree_->skeleton_simplex_range(max_dimension)};
    SimplexList skeleton_simplices{};
    std::transform(
        simplices.begin(),
        simplices.end(),
        std::back_inserter(skeleton_simplices),
        [this](const auto &simplex_handle)
        {
            return Simplex{get_vertices(simplex_handle)};
        });
    return skeleton_simplices;
}