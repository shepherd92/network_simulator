#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#include "typedefs.h"

class Simplex;
class SimplexHash;

using SimplexList = std::vector<Simplex>;
using SimplexSet = std::unordered_set<Simplex, SimplexHash>;

class Simplex
{
public:
    Simplex(const VertexList &vertices);
    const VertexList &vertices() const;
    bool is_valid() const;
    std::vector<uint32_t> hash() const;

    bool is_face(const Simplex &other) const;
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
        VertexList &current_combination,
        const uint32_t array_index) const;

    VertexList vertices_;
};

struct SimplexHash
{
    std::size_t operator()(const Simplex &simplex) const;
};

SimplexList create_simplices(const std::vector<VertexList> &simplices_in);
ISimplexList create_raw_simplices(const SimplexList &simplices_in);
SimplexList get_skeleton_simplices(const SimplexList &simplices_in, const Dimension dimension);
SimplexList select_simplices_by_dimension(const SimplexList &simplices, const Dimension dimension);
SimplexList select_higher_dimensional_simplices(const SimplexList &simplices, const Dimension dimension);
void sort_simplices(SimplexList &simplices, const bool ascending);
std::vector<Dimension> calc_dimension_distribution(const ISimplexList &simplices_in);
std::vector<Dimension> calc_dimension_distribution(const SimplexList &simplices);

#endif