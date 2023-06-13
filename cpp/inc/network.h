#ifndef _NETWORK_H_
#define _NETWORK_H_

#include <cstdint>
#include <vector>

#include "typedefs.h"

class Network
{
public:
    Network() = default;
    explicit Network(const SimplicialComplex &simplicial_complex);
    explicit Network(const std::vector<Vertex> &vertices, const std::vector<Simplex> &interactions);
    explicit Network(const std::vector<vertex_id> &vertices, const std::vector<std::vector<vertex_id>> &interactions);

    void add_simplices(const std::vector<std::vector<vertex_id>> &simplices);

    std::vector<uint32_t> calc_degree_sequence(
        const dimension simplex_dimension,
        const dimension neighbor_dimension) const;
    std::vector<Simplex> get_simplices_by_dimension(const dimension dimension) const;

    SimplexIterator simplices() const;
    const std::vector<Simplex> &facets();
    const std::vector<Simplex> &interactions();
    void interactions(const std::vector<Simplex> &interactions);
    void max_dimension(const dimension &max_dimension);

    uint32_t num_vertices() const;
    uint32_t num_edges() const;
    uint32_t num_triangles() const;
    uint32_t num_simplices() const;

private:
    void combinations(
        const std::vector<int32_t> &elements,
        const uint32_t k,
        std::set<std::vector<int32_t>> &subarrays,
        std::vector<int32_t> &out,
        const uint32_t i);

    std::vector<Simplex> get_simplex_skeleton_for_max_dimension(const Simplex &simplex) const;

    SimplicialComplex simplicial_complex_;
    std::vector<Simplex> interactions_;
    std::vector<Simplex> facets_;
    dimension max_dimension_;
};

#endif