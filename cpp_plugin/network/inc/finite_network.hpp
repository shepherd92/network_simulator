#ifndef _FINITE_NETWORK_HPP_
#define _FINITE_NETWORK_HPP_

#include <iostream>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Simplex_tree.h>

#include "network.hpp"
#include "typedefs.hpp"

class FiniteNetwork : virtual public Network
{
public:
    // Gudhi aliases
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
    using VertexHandle = SimplexTreeOptions::Vertex_handle;
    using SimplexHandle = SimplexTree::Simplex_handle;
    using SimplexHandleList = std::vector<SimplexHandle>;
    using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
    using PersistentCohomology = Gudhi::persistent_cohomology::Persistent_cohomology<SimplexTree, Field_Zp>;

public:
    FiniteNetwork();
    ~FiniteNetwork();
    FiniteNetwork(FiniteNetwork &&other) noexcept;
    FiniteNetwork &operator=(FiniteNetwork &&other) noexcept;

    void reset() override;

    std::vector<uint32_t> calc_coface_degree_sequence(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension) override;
    std::vector<int32_t> calc_betti_numbers();
    SimplexList get_skeleton(const Dimension max_dimension) override;

protected:
    void assert_simplicial_complex_is_initialized();
    void assert_simplicial_complex_is_built();
    PointIdList get_simplex_vertices(const SimplexHandle &simplex_handle);
    PersistentCohomology &get_persistence();
    void reset_persistence();

    std::optional<SimplexTree> simplex_tree_;
    PersistentCohomology *persistent_cohomology_;

private:
    void create_simplicial_complex();

    // simplicial complex methods
    void calc_persistent_cohomology();
    bool is_valid() const;
    void add_vertices();
    virtual void fill_simplicial_complex() = 0;
    void reset_simplicial_complex();
};

#endif