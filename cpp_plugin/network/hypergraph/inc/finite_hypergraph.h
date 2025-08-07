#ifndef _FINITE_HYPERGRAPH_H_
#define _FINITE_HYPERGRAPH_H_

#include <iostream>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Simplex_tree.h>

#include "hypergraph.h"
#include "finite_network.h"
#include "typedefs.h"

class FiniteHypergraph : public FiniteNetwork, public Hypergraph
{
public:
    FiniteHypergraph(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const ISimplexList &interactions,
        const bool weighted);
    FiniteHypergraph(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const SimplexList &interactions,
        const bool weighted);
    FiniteHypergraph(const FiniteHypergraph &other);
    ~FiniteHypergraph();

    bool weighted() const;

    std::vector<uint32_t> calc_simplex_interaction_degree_sequence(
        const Dimension simplex_dimension) override;
    std::vector<std::vector<std::pair<float, float>>> calc_persistence_intervals();
    std::vector<ISimplexList> calc_persistence_pairs();

private:
    SimplexList calc_simplices(const Dimension dimension) override;

    void calc_persistent_cohomology();
    void add_vertices(const PointIdList &vertices);
    void fill_simplicial_complex() override;
    std::vector<uint32_t> calc_vertex_interaction_degree_distribution() const override;

    const bool weighted_;
};

#endif