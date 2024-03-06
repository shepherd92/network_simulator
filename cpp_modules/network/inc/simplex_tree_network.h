#ifndef _SIMPLEX_TREE_NETWORK_H_
#define _SIMPLEX_TREE_NETWORK_H_

#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Simplex_tree.h>

#include "network.h"
#include "typedefs.h"

class SimplexTreeNetwork : public Network
{
public:
    SimplexTreeNetwork(const Dimension max_dimension);
    ~SimplexTreeNetwork();
    virtual void expand();
    Dimension get_max_dimension() const;
    void set_max_dimension(const Dimension max_dimension);

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

    bool is_valid() const override;
    SimplexTree &get_simplex_tree();
    void reset_simplicial_complex() override;

    Dimension max_dimension_;
    std::optional<SimplexTree> simplex_tree_;

private:
    virtual SimplexHandleList get_simplices() = 0;

    SimplexList convert_to_representable_simplices(const SimplexList &simplices_in) const override;
    void add_simplex(const VertexList &simplex) override;

    void calc_facets_simplicial_complex() override;

    void initialize_simplicial_complex_if_needed() override;
    SimplexList get_skeleton_simplicial_complex(const Dimension max_dimension) override;

    void filter_simplicial_complex(const VertexList &vertices) override;
    std::vector<uint32_t> calc_degree_sequence_simplicial_complex(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension) override;
    std::vector<Dimension> calc_simplex_dimension_distribution() override;

    VertexList get_vertices(const SimplexHandle &simplex_handle);
};

#endif // _SIMPLEX_TREE_NETWORK_H_