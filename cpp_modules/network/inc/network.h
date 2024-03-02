#ifndef _NETWORK_H_
#define _NETWORK_H_

#include <gudhi/Simplex_tree.h>
#include <optional>

#include "typedefs.h"

class Network
{
public:
    Network(const Dimension max_dimension);

    void add_simplices(const ISimplexList &simplices);
    void add_simplex(const ISimplex &simplex);

    void keep_only_vertices(const VertexList &vertices);
    void reset();

    void expand();
    uint32_t num_simplices();
    uint32_t num_vertices() const;

    std::vector<std::vector<uint32_t>> calc_persistence_pairs();
    std::vector<uint32_t> calc_simplex_dimension_distribution();
    std::vector<uint32_t> calc_degree_sequence(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension);

    ISimplexList get_simplices();
    ISimplexList get_skeleton(const Dimension max_dimension);
    void set_simplices(const ISimplexList &simplices);

    Dimension get_max_dimension() const;
    void set_max_dimension(const Dimension max_dimension);

    ISimplexList get_facets();
    void set_facets(const ISimplexList &facets);

private:
    template <typename Iterator>
    ISimplexList convert_to_raw_simplices(const Iterator &simplex_range);

    void combination_util(
        const ISimplex &simplex,
        const uint32_t combination_index,
        ISimplexList &result,
        ISimplex &current_combination,
        const uint32_t array_index) const;

    ISimplexList calc_facets();
    void filter_simplex_tree(const VertexList &vertices);
    void filter_facets(const VertexList &vertices);
    void filter_interactions(const VertexList &vertices);

    Dimension max_dimension_;
    SimplexTree simplex_tree_;
    std::optional<ISimplexList> interactions_;
    std::optional<ISimplexList> facets_;
};

#endif