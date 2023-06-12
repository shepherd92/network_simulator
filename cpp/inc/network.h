#ifndef _NETWORK_H_
#define _NETWORK_H_

#include <cstdint>
#include <vector>

#include "typedefs.h"

class Network
{
public:
    Network(const simplicial_complex &simplicial_complex);

    std::vector<int32_t> calc_degree_sequence(
        const dimension simplex_dimension,
        const dimension neighbor_dimension) const;

    auto facets() const;
    auto simplices() const;

private:
    std::vector<simplex> select_simplices_by_dimension(
        const dimension dimension);

    void combinations(
        const std::vector<int32_t> &elements,
        const uint32_t k,
        std::set<std::vector<int32_t>> &subarrays,
        std::vector<int32_t> &out,
        const uint32_t i);

    simplicial_complex simplicial_complex_;
    std::vector<simplex> interactions;
};

#endif