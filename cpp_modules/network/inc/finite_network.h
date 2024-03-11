#ifndef _FINITE_NETWORK_H_
#define _FINITE_NETWORK_H_

#include <iostream>
#include <gudhi/Persistent_cohomology.h>

#include "network.h"
#include "typedefs.h"

class FiniteNetwork : public Network
{
public:
    FiniteNetwork(const Dimension max_dimension, const VertexList &vertices, const ISimplexList &interactions);
    ~FiniteNetwork();

    void add_vertices(const VertexList &vertices) override;
    void expand() override;

    uint32_t num_simplices() override;
    std::vector<ISimplexList> calc_persistence_pairs();
    std::vector<int32_t> calc_betti_numbers();

private:
    using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
    using PersistentCohomology = Gudhi::persistent_cohomology::Persistent_cohomology<SimplexTree, Field_Zp>;

    void add_simplices(const SimplexList &simplices) override;
    SimplexHandleList get_simplices() override;

    SimplexList get_skeleton_interactions(const Dimension max_dimension) override;
    SimplexList get_skeleton_simplicial_complex(const Dimension max_dimension) override;

    void reset_simplicial_complex() override;
    void reset_persistence();
    const PersistentCohomology &get_persistence();
    void calc_persistent_cohomology();

    PersistentCohomology *persistent_cohomology_;
};

#endif