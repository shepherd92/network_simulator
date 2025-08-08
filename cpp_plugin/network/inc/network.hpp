#ifndef _NETWORK_HPP_
#define _NETWORK_HPP_

#include <map>
#include <optional>

#include "simplex_list.hpp"
#include "typedefs.hpp"

class Network
{
public:
    Network(const Dimension max_dimension, const PointIdList &vertices);

    Network(Network &&other) noexcept
    {
        max_dimension_ = std::move(other.max_dimension_);
        vertices_ = std::move(other.vertices_);
        simplices_ = std::move(other.simplices_);
    }
    Network &operator=(Network &&other) noexcept
    {
        if (this != &other)
        {
            max_dimension_ = std::move(other.max_dimension_);
            vertices_ = std::move(other.vertices_);
            simplices_ = std::move(other.simplices_);
        }
        return *this;
    }

    uint32_t num_simplices(const Dimension dimension);

    virtual void reset();

    std::vector<Dimension> calc_simplex_dimension_distribution();
    virtual std::vector<uint32_t> calc_coface_degree_sequence(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension) = 0;

    virtual void keep_only_vertices(const PointIdList &vertices);
    ISimplexList get_skeleton_interface(const Dimension max_dimension);

    uint32_t num_vertices();

protected:
    virtual SimplexList get_skeleton(const Dimension max_dimension) = 0;
    virtual SimplexList calc_simplices(const Dimension dimension) = 0;
    const SimplexList &get_simplices(const Dimension dimension);
    void set_simplices(const Dimension dimension, const SimplexList &simplices);

    Dimension max_dimension_;
    PointIdList vertices_;

private:
    std::map<PointId, std::vector<PointId>> create_vertex_simplex_map(const SimplexList &simplices) const;
    std::vector<std::optional<SimplexList>> simplices_;
};

#endif
