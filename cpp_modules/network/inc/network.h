#ifndef _NETWORK_H_
#define _NETWORK_H_

#include <optional>

#include "simplex.h"
#include "typedefs.h"

class Network
{
public:
    Network(const Dimension max_dimension);
    ~Network();

    void add_vertices(const VertexList &vertices);
    void add_simplices_interface(const ISimplexList &simplices);
    void add_simplices(const SimplexList &simplices);

    void keep_only_vertices(const VertexList &vertices);
    void reset();

    void expand();
    virtual uint32_t num_simplices() = 0;
    uint32_t num_vertices();

    std::vector<int32_t> calc_betti_numbers();
    std::vector<ISimplexList> calc_persistence_pairs();

    std::vector<Dimension> calc_simplex_dimension_distribution();
    std::vector<Dimension> calc_facet_dimension_distribution();
    std::vector<Dimension> calc_interaction_dimension_distribution();

    std::vector<uint32_t> calc_degree_sequence(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension);

    ISimplexList get_simplices_interface();
    virtual SimplexHandleList get_simplices() = 0;
    void set_simplices(const ISimplexList &simplices);
    ISimplexList get_skeleton_interface(const Dimension max_dimension);

    Dimension get_max_dimension() const;
    void set_max_dimension(const Dimension max_dimension);

    const VertexList &get_vertices_interface() const;
    void set_vertices(const VertexList &vertices);

    ISimplexList get_interactions_interface() const;
    void set_interactions(const ISimplexList &interactions);

    ISimplexList get_facets_interface();
    void set_facets(const ISimplexList &facets);

protected:
    template <typename Iterator>
    ISimplexList convert_to_raw_simplices(const Iterator &simplex_range);
    VertexList convert_to_raw_simplex(const SimplexHandle &simplex);

    void create_simplex_tree_from_interactions();

    SimplexTree &get_simplex_tree();
    const SimplexList &get_interactions() const;
    const SimplexList &get_facets();

    Dimension max_dimension_;
    VertexList vertices_;
    std::optional<SimplexTree> simplex_tree_;
    std::optional<SimplexList> interactions_;
    std::optional<SimplexList> facets_;

private:
    void add_simplex(const VertexList &simplex);
    VertexList get_vertices(const SimplexHandle &simplex_handle);
    ISimplexList get_skeleton(const ISimplexList &simplices) const;
    void sort_interactions(const bool ascending);

    void calc_facets_simplex_tree();
    void calc_facets_interactions();

    std::vector<uint32_t> calc_degree_sequence_interactions(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension);
    std::vector<uint32_t> calc_degree_sequence_simplex_tree(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension);

    std::vector<Dimension> calc_dimension_distribution(const SimplexList &simplices) const;
    std::vector<Dimension> calc_dimension_distribution(const ISimplexList &simplices) const;
    std::vector<Dimension> calc_dimension_distribution(const SimplexHandleList &simplices);

    void reset_simplex_tree();
    void reset_persistence();

    void filter_simplex_tree(const VertexList &vertices);
    void filter_facets(const VertexList &vertices);
    void filter_interactions(const VertexList &vertices);

    const PersistentCohomology &get_persistence();
    void calc_persistent_cohomology();

    PersistentCohomology *persistent_cohomology_;
};

#endif