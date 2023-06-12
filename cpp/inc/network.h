#ifndef _NETWORK_H_
#define _NETWORK_H_

#include <cstdint>
#include <vector>

#include "typedefs.h"

class Network
{
public:
    Network(const SimplicialComplex &simplicial_complex);
    Network(const std::vector<Simplex> &interactions);

    auto calc_degree_sequence(
        const dimension simplex_dimension,
        const dimension neighbor_dimension) const;
    auto get_simplices_by_dimension(const dimension dimension) const;

    auto facets() const;
    auto simplices() const;

    auto num_vertices() const;
    auto num_edges() const;
    auto num_triangles() const;
    auto num_simplices() const;

private:
    void combinations(
        const std::vector<int32_t> &elements,
        const uint32_t k,
        std::set<std::vector<int32_t>> &subarrays,
        std::vector<int32_t> &out,
        const uint32_t i);

    SimplicialComplex simplicial_complex_;
    std::vector<Simplex> interactions;
};

#endif