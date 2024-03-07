#ifndef _NETWORK_H_
#define _NETWORK_H_

#include <gudhi/Simplex_tree.h>
#include <optional>

#include "simplex.h"
#include "typedefs.h"

class Network
{
public:
    Network(const VertexList &vertices, const ISimplexList &interactions);

    void create_simplicial_complex(const Dimension max_dimension);
    virtual uint32_t num_simplices() = 0;

    void keep_only_vertices(const VertexList &vertices);
    void reset();

    std::vector<Dimension> calc_simplex_dimension_distribution();
    std::vector<Dimension> calc_facet_dimension_distribution();
    std::vector<Dimension> calc_interaction_dimension_distribution();

    std::vector<uint32_t> calc_degree_sequence(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension);

    ISimplexList get_skeleton_interface(const Dimension max_dimension);
    virtual void expand(const Dimension max_dimension);

    uint32_t num_vertices();
    const VertexList &get_vertices_interface() const;
    void set_vertices(const VertexList &vertices);

    ISimplexList get_interactions_interface() const;
    void set_interactions(const ISimplexList &interactions);

    ISimplexList get_facets_interface();
    ISimplexList get_simplices_interface();

protected:
    struct SimplexTreeOptions
    {
        typedef Gudhi::linear_indexing_tag Indexing_tag;
        typedef int Vertex_handle;
        typedef float Filtration_value;
        typedef uint32_t Simplex_key;
        static const bool store_key = true;
        static const bool store_filtration = false;
        static const bool contiguous_vertices = false;
    };

    using SimplexTree = Gudhi::Simplex_tree<SimplexTreeOptions>;
    using SimplexHandle = SimplexTree::Simplex_handle;
    using SimplexHandleList = std::vector<SimplexHandle>;

    bool is_valid() const;

    virtual void add_simplices(const SimplexList &simplices, const Dimension dimension);
    virtual void add_vertices(const VertexList &vertices);
    VertexList get_vertices(const SimplexHandle &simplex_handle);
    void calc_facets();

    virtual void reset_simplicial_complex();

    SimplexList get_skeleton(const Dimension max_dimension);

    const SimplexList &get_interactions() const;
    const SimplexList &get_facets();

    VertexList vertices_;
    SimplexList interactions_;
    std::optional<SimplexList> facets_;
    std::optional<SimplexTree> simplex_tree_;

private:
    void initialize_simplicial_complex_if_needed();

    virtual SimplexHandleList get_simplices() = 0;
    virtual SimplexList get_skeleton_simplicial_complex(const Dimension max_dimension) = 0;
    virtual SimplexList get_skeleton_interactions(const Dimension max_dimension) = 0;

    template <typename Iterator>
    ISimplexList convert_to_raw_simplices(const Iterator &simplex_range);

    std::vector<uint32_t> calc_degree_sequence_interactions(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension);
    std::vector<uint32_t> calc_degree_sequence_simplicial_complex(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension);

    void filter_facets(const VertexList &vertices);
    void filter_interactions(const VertexList &vertices);
};

#endif