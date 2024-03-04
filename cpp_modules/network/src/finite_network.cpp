#include "finite_network.h"

FiniteNetwork::FiniteNetwork(const Dimension max_dimension)
    : Network{max_dimension}
{
}

SimplexHandleList FiniteNetwork::get_simplices()
{
    return simplex_tree_->filtration_simplex_range();
}

uint32_t FiniteNetwork::num_simplices()
{
    return simplex_tree_->num_simplices();
}
