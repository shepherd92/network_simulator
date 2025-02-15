#ifndef _FINITE_NETWORK_H_
#define _FINITE_NETWORK_H_

#include <iostream>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Simplex_tree.h>

#include "network.h"
#include "typedefs.h"

class FiniteNetwork : public Network
{
public:
    FiniteNetwork(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const ISimplexList &interactions,
        const bool weighted);
    FiniteNetwork(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const SimplexList &interactions,
        const bool weighted);
    FiniteNetwork(const FiniteNetwork &other);
    ~FiniteNetwork();
    void create_simplicial_complex();

    bool weighted() const;
    FiniteNetwork filter(const PointIdList &vertices) const;
    void expand();
    void reset() override;

    std::vector<uint32_t> calc_coface_degree_sequence(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension) override;
    std::vector<uint32_t> calc_simplex_interaction_degree_sequence(
        const Dimension simplex_dimension) override;
    PointIdList get_vertices() const;
    std::vector<int32_t> calc_betti_numbers();
    std::vector<ISimplexList> calc_persistence_pairs();
    std::vector<std::vector<std::pair<float, float>>> calc_persistence_intervals();

private:
    struct SimplexTreeOptions
    {
        typedef Gudhi::linear_indexing_tag Indexing_tag;
        typedef int Vertex_handle;
        typedef float Filtration_value;
        typedef uint32_t Simplex_key;
        static const bool store_key = true;
        static const bool store_filtration = true;
        static const bool contiguous_vertices = false;
        static const bool link_nodes_by_label = false;
        static const bool stable_simplex_handles = false;
    };
    using SimplexTree = Gudhi::Simplex_tree<SimplexTreeOptions>;
    using SimplexHandle = SimplexTree::Simplex_handle;
    using SimplexHandleList = std::vector<SimplexHandle>;
    using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
    using PersistentCohomology = Gudhi::persistent_cohomology::Persistent_cohomology<SimplexTree, Field_Zp>;

    SimplexList calc_simplices(const Dimension dimension) override;
    SimplexList get_skeleton(const Dimension max_dimension) override;

    void calc_persistent_cohomology();
    bool is_valid() const;
    void assert_simplicial_complex_is_initialized();
    void assert_simplicial_complex_is_built();
    void add_vertices(const PointIdList &vertices);
    PointIdList get_simplex_vertices(const SimplexHandle &simplex_handle);
    void fill_simplicial_complex();
    void reset_simplicial_complex();
    void reset_persistence();
    PersistentCohomology &get_persistence();
    std::vector<uint32_t> calc_vertex_interaction_degree_distribution() const override;

    const bool weighted_;
    std::optional<SimplexTree> simplex_tree_;
    PersistentCohomology *persistent_cohomology_;
};

#endif