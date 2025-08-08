#ifndef _HYPERGRAPH_HPP_
#define _HYPERGRAPH_HPP_

#include <map>
#include <optional>

#include "network.hpp"
#include "simplex_list.hpp"
#include "typedefs.hpp"

class Hypergraph : virtual public Network
{
public:
    // constructors, assignment operators, and destructors
    Hypergraph(const SimplexList &interactions);
    Hypergraph(Hypergraph &&other) noexcept;
    Hypergraph &operator=(Hypergraph &&other) noexcept;

    virtual std::vector<Dimension> calc_interaction_dimension_distribution() const;
    virtual std::vector<uint32_t> calc_simplex_interaction_degree_sequence(
        const Dimension simplex_dimension) = 0;

    // interface functions exposed to Python
    void set_interactions(const ISimplexList &interactions);
    ISimplexList get_interactions() const;

protected:
    SimplexList interactions_;

private:
    virtual std::vector<uint32_t> calc_vertex_interaction_degree_distribution() const = 0;
};

#endif
