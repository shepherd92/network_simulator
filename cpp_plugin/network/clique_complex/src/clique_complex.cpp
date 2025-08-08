#include <atomic>
#include <cassert>
#include <execution>
#include <mutex>

#include "hypergraph.hpp"
#include "simplex.hpp"
#include "simplex_list.hpp"
#include "tools.hpp"
#include "typedefs.hpp"

CliqueComplex::CliqueComplex()
{
}

CliqueComplex::CliqueComplex(CliqueComplex &&other) noexcept
{
}

CliqueComplex &CliqueComplex::operator=(CliqueComplex &&other) noexcept
{
    return *this;
}

void CliqueComplex::keep_only_vertices(const PointIdList &vertices)
{
    Network::keep_only_vertices(vertices);
    edges_ = edges_.filter(vertices);
}
