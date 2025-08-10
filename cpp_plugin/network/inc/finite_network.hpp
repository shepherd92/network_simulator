#ifndef _FINITE_NETWORK_HPP_
#define _FINITE_NETWORK_HPP_

#include <iostream>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Simplex_tree.h>

#include "network.hpp"
#include "typedefs.hpp"

class FiniteNetwork : virtual public Network
{
public:
    FiniteNetwork();
    FiniteNetwork(FiniteNetwork &&other) noexcept;
    FiniteNetwork &operator=(FiniteNetwork &&other) noexcept;

    std::vector<uint32_t> calc_coface_degree_sequence(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension) override;
    std::vector<int32_t> calc_betti_numbers();
    SimplexList get_skeleton(const Dimension max_dimension) override;
};

#endif