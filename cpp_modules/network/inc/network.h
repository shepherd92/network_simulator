#ifndef _NETWORK_H_
#define _NETWORK_H_

#include <optional>

#include "simplex.h"
#include "typedefs.h"

class Network
{
public:
    Network();

    virtual void add_vertices(const VertexList &vertices);
    virtual void add_simplices_interface(const ISimplexList &simplices);

    void keep_only_vertices(const VertexList &vertices);
    void reset();

    virtual uint32_t num_simplices() = 0;
    uint32_t num_vertices();

    virtual std::vector<Dimension> calc_simplex_dimension_distribution() = 0;
    std::vector<Dimension> calc_facet_dimension_distribution();
    std::vector<Dimension> calc_interaction_dimension_distribution();

    std::vector<uint32_t> calc_degree_sequence(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension);

    ISimplexList get_skeleton_interface(const Dimension max_dimension);

    const VertexList &get_vertices_interface() const;
    void set_vertices(const VertexList &vertices);

    ISimplexList get_interactions_interface() const;
    void set_interactions(const ISimplexList &interactions);

    ISimplexList get_facets_interface();
    void set_facets(const ISimplexList &facets);

protected:
    template <typename Iterator>
    ISimplexList convert_to_raw_simplices(const Iterator &simplex_range);

    void create_simplicial_complex_from_interactions();

    const SimplexList &get_interactions() const;
    const SimplexList &get_facets();

    VertexList vertices_;
    std::optional<SimplexList> interactions_;
    std::optional<SimplexList> facets_;

private:
    virtual bool is_valid() const = 0;
    virtual void initialize_simplicial_complex_if_needed() = 0;

    void add_simplices(const SimplexList &simplices);
    virtual void add_simplex(const VertexList &simplex) = 0;

    virtual SimplexList get_skeleton_simplicial_complex(const Dimension max_dimension) = 0;
    SimplexList get_skeleton(const Dimension max_dimension);
    void sort_interactions(const bool ascending);

    virtual SimplexList convert_to_representable_simplices(const SimplexList &simplices_in) const = 0;

    void calc_facets_interactions();
    virtual void calc_facets_simplicial_complex() = 0;

    std::vector<uint32_t> calc_degree_sequence_interactions(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension);
    virtual std::vector<uint32_t> calc_degree_sequence_simplicial_complex(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension) = 0;

    virtual void reset_simplicial_complex() = 0;
    void reset_persistence();

    virtual void filter_simplicial_complex(const VertexList &vertices) = 0;
    void filter_facets(const VertexList &vertices);
    void filter_interactions(const VertexList &vertices);
};

#endif