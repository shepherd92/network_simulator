#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#include <bitset>

#include "typedefs.h"

class Simplex;
class SimplexHash;

using SimplexList = std::vector<Simplex>;
using SimplexSet = std::unordered_set<Simplex, SimplexHash>;

constexpr uint32_t BLOOM_FILTER_SIZE{64};

class Simplex
{
public:
    Simplex(const PointIdList &vertices);
    const PointIdList &vertices() const;
    bool is_valid() const;
    std::vector<uint32_t> hash() const;

    inline bool is_face(const Simplex &other) const;
    inline std::bitset<BLOOM_FILTER_SIZE> bloom_filter() const;
    SimplexList get_skeleton(const Dimension dimension) const;
    Dimension dimension() const;
    Simplex operator-(const Simplex &other) const;
    bool operator==(const Simplex &other) const;
    friend std::ostream &operator<<(std::ostream &os, const Simplex &simplex);

private:
    void combination_util(
        const Dimension dimension,
        const uint32_t combination_index,
        SimplexList &result,
        PointIdList &current_combination,
        const uint32_t array_index) const;

    PointIdList vertices_;
    std::bitset<BLOOM_FILTER_SIZE> bloom_filter_;
};

struct SimplexHash
{
    std::size_t operator()(const Simplex &simplex) const;
};

inline bool Simplex::is_face(const Simplex &other) const
{
    // if (!is_valid())
    // {
    //     return false;
    // }
    if ((bloom_filter_ & other.bloom_filter()) != bloom_filter_)
    {
        return false;
    }
    return std::includes(other.vertices().begin(), other.vertices().end(), vertices_.begin(), vertices_.end());
}

inline std::bitset<BLOOM_FILTER_SIZE> Simplex::bloom_filter() const
{
    return bloom_filter_;
}

SimplexList create_simplices(const std::vector<PointIdList> &simplices_in);
ISimplexList create_raw_simplices(const SimplexList &simplices_in);
SimplexList create_simplices_from_connections(const ConnectionList &connections);
ISimplexList filter_simplices_interface(const ISimplexList &simplices, const PointIdList &vertices_to_keep);
SimplexList filter_simplices(const SimplexList &simplices, const PointIdList &vertices_to_keep);
SimplexList get_skeleton_simplices(const SimplexList &simplices_in, const Dimension dimension);
SimplexList select_simplices_by_dimension(const SimplexList &simplices, const Dimension dimension);
SimplexList select_higher_dimensional_simplices(const SimplexList &simplices, const Dimension dimension);
void sort_simplices(SimplexList &simplices, const bool ascending);
std::vector<Dimension> calc_dimension_distribution(const ISimplexList &simplices_in);
std::vector<Dimension> calc_dimension_distribution(const SimplexList &simplices);

#endif