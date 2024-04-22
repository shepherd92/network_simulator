#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#include <bitset>

#include "typedefs.h"

struct SimplexHash;

constexpr auto BLOOM_FILTER_SIZE{64U};

class Simplex
{
public:
    Simplex(const PointIdList &vertices);
    const PointIdList &vertices() const;

    inline bool is_face(const Simplex &other) const;
    std::vector<Simplex> faces(const Dimension dimension) const;
    std::vector<Simplex> skeleton(const Dimension max_dimension) const;
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
        std::vector<Simplex> &result,
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

#endif