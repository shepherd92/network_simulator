#ifndef _NETWORK_HPP_
#define _NETWORK_HPP_

#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Simplex_tree.h>
#include <map>
#include <optional>

#include "simplex_list.hpp"
#include "typedefs.hpp"

class Network
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
    Network(const Dimension max_dimension, const PointIdList &vertices);
    ~Network();

    Network(Network &&other) noexcept;
    Network &operator=(Network &&other) noexcept;

    virtual void reset();

    std::vector<Dimension> calc_simplex_dimension_distribution();
    virtual std::vector<uint32_t> calc_coface_degree_sequence(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension) = 0;

    ISimplexList get_skeleton_interface(const Dimension max_dimension);

    uint32_t num_vertices();
    const PointIdList &get_vertices() const;
    virtual void set_vertices(const PointIdList &vertices);

    uint32_t num_simplices(const Dimension dimension);
    const SimplexList &get_simplices(const Dimension dimension);
    Dimension get_max_dimension() const;

protected:
    virtual SimplexList get_skeleton(const Dimension max_dimension) = 0;
    virtual SimplexList calc_simplices(const Dimension dimension) = 0;

    void assert_simplex_tree_is_initialized();
    void assert_simplex_tree_is_built();
    PointIdList get_simplex_vertices_simplex_tree(const SimplexHandle &simplex_handle);
    PersistentCohomology &get_persistence();
    void reset_persistence();

    Dimension max_dimension_;
    PointIdList vertices_;
    std::optional<SimplexTree> simplex_tree_;
    PersistentCohomology *persistent_cohomology_;

private:
    // simplicial complex methods
    void create_simplex_tree();
    void calc_persistent_cohomology();
    bool is_simplex_tree_valid() const;
    void add_vertices_to_simplex_tree();
    virtual void fill_simplex_tree() = 0;
    void reset_simplex_tree();

    std::map<PointId, std::vector<PointId>> create_vertex_simplex_map(const SimplexList &simplices) const;
    std::vector<std::optional<SimplexList>> simplices_;
};

#endif
