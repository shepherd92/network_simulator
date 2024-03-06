#ifndef _SIMPLEX_TREE_FINITE_NETWORK_H_
#define _SIMPLEX_TREE_FINITE_NETWORK_H_

#include "simplex_tree_network.h"
#include "typedefs.h"

class SimplexTreeFiniteNetwork : public SimplexTreeNetwork
{
public:
    SimplexTreeFiniteNetwork(const Dimension max_dimension);
    ~SimplexTreeFiniteNetwork();

    void add_vertices(const VertexList &vertices) override;
    void add_simplices_interface(const ISimplexList &simplices) override;

    void expand() override;

    uint32_t num_simplices() override;
    std::vector<ISimplexList> calc_persistence_pairs();
    std::vector<int32_t> calc_betti_numbers();

private:
    using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
    using PersistentCohomology = Gudhi::persistent_cohomology::Persistent_cohomology<SimplexTree, Field_Zp>;

    SimplexHandleList get_simplices() override;
    void reset_simplicial_complex() override;
    void reset_persistence();
    const PersistentCohomology &get_persistence();
    void calc_persistent_cohomology();

    PersistentCohomology *persistent_cohomology_;
};

SimplexTreeFiniteNetwork::SimplexTreeFiniteNetwork(const Dimension max_dimension)
    : SimplexTreeNetwork{max_dimension},
      persistent_cohomology_{nullptr}
{
}

SimplexTreeFiniteNetwork::~SimplexTreeFiniteNetwork()
{
    reset_persistence();
}

void SimplexTreeFiniteNetwork::add_vertices(const VertexList &vertices)
{
    Network::add_vertices(vertices);
    reset_persistence();
}

void SimplexTreeFiniteNetwork::add_simplices_interface(const ISimplexList &simplices)
{
    Network::add_simplices_interface(simplices);
    reset_persistence();
}

void SimplexTreeFiniteNetwork::reset_simplicial_complex()
{
    SimplexTreeNetwork::reset_simplicial_complex();
    reset_persistence();
}

void SimplexTreeFiniteNetwork::reset_persistence()
{
    delete persistent_cohomology_;
    persistent_cohomology_ = nullptr;
}

void SimplexTreeFiniteNetwork::expand()
{
    SimplexTreeNetwork::expand();
    reset_persistence();
}

void SimplexTreeFiniteNetwork::calc_persistent_cohomology()
{
    auto &simplex_tree{get_simplex_tree()};
    reset_persistence();
    persistent_cohomology_ = new PersistentCohomology{simplex_tree};
    persistent_cohomology_->init_coefficients(2);
    persistent_cohomology_->compute_persistent_cohomology();
}

SimplexTreeNetwork::SimplexHandleList SimplexTreeFiniteNetwork::get_simplices()
{
    assert(is_valid());
    return simplex_tree_->filtration_simplex_range();
}

const SimplexTreeFiniteNetwork::PersistentCohomology &SimplexTreeFiniteNetwork::get_persistence()
{
    if (!persistent_cohomology_)
    {
        calc_persistent_cohomology();
    }
    return *persistent_cohomology_;
}

std::vector<ISimplexList> SimplexTreeFiniteNetwork::calc_persistence_pairs()
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
            result.emplace_back({birth_simplex, death_simplex});
        });

    if (total > 10000U)
    {
        std::cout << "\rCalc persistence pairs ... " << total << " / " << total;
    }

    return result;
}

std::vector<int32_t> SimplexTreeFiniteNetwork::calc_betti_numbers()
{
    const auto &persistent_cohomology{get_persistence()};
    const auto result{persistent_cohomology.betti_numbers()};
    return result;
}

#endif