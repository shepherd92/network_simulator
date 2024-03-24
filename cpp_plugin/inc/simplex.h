#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#include <bitset>

#include "typedefs.h"

class Simplex;
struct SimplexHash;

constexpr auto BLOOM_FILTER_SIZE{64U};

class Simplex
{
public:
    Simplex(const PointIdList &vertices);
    const PointIdList &vertices() const;

    inline bool is_face(const Simplex &other) const;
    SimplexList faces(const Dimension dimension) const;
    SimplexList skeleton(const Dimension max_dimension) const;
    std::bitset<BLOOM_FILTER_SIZE> bloom_filter() const;
    Dimension dimension() const;

    Simplex operator-(const Simplex &other) const;
    bool operator==(const Simplex &other) const;
    bool operator<(const Simplex &other) const;
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
    uint64_t operator()(const Simplex &simplex) const;
};

#include "simplex.inl"

SimplexList create_simplices(const std::vector<PointIdList> &simplices_in);
ISimplexList create_raw_simplices(const SimplexList &simplices_in);
SimplexList filter_simplices(const SimplexList &simplices, const PointIdList &vertices_to_keep);
SimplexList get_faces_simplices(const SimplexList &simplices_in, const Dimension dimension);
SimplexList get_cofaces(const SimplexList &simplices_in, const Simplex &simplex);
SimplexList get_skeleton_simplices(const SimplexList &simplices, const Dimension max_dimension);
SimplexList select_simplices_by_dimension(const SimplexList &simplices, const Dimension dimension);
SimplexList select_higher_dimensional_simplices(const SimplexList &simplices, const Dimension dimension);
void sort_simplices(SimplexList &simplices, const bool ascending);
std::vector<Dimension> calc_dimension_distribution(const ISimplexList &simplices_in);
std::vector<Dimension> calc_dimension_distribution(const SimplexList &simplices);

#endif