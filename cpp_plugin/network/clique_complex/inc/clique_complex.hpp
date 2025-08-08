#ifndef _HYPERGRAPH_HPP_
#define _HYPERGRAPH_HPP_

#include "network.hpp"
#include "typedefs.hpp"

class CliqueComplex : virtual public Network
{
public:
    // constructors, assignment operators, and destructors
    CliqueComplex();
    CliqueComplex(CliqueComplex &&other) noexcept;
    CliqueComplex &operator=(CliqueComplex &&other) noexcept;

    void keep_only_vertices(const PointIdList &vertices) override;
};

#endif