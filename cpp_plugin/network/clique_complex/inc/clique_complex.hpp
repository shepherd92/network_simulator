#ifndef _CLIQUE_COMPLEX_HPP_
#define _CLIQUE_COMPLEX_HPP_

#include "network.hpp"
#include "typedefs.hpp"

class CliqueComplex : virtual public Network
{
public:
    CliqueComplex(const ConnectionList &edges);

    // move assignment defined due to virtual base class
    CliqueComplex(CliqueComplex &&other) noexcept;
    CliqueComplex &operator=(CliqueComplex &&other) noexcept;

    void set_edges(const ConnectionList &edges);
    ConnectionList get_edges() const;

    void set_vertices(const PointIdList &vertices) override;

protected:
    void expand();
    void fill_simplex_tree() override;

    ConnectionList filter_edges(const PointIdList &vertices) const;
    ConnectionList edges_; // every vertex is by default connected to the typical vertex
};

#endif