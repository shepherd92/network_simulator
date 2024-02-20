#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#include "typedefs.h"

class Simplex
{
public:
    Simplex(const VertexList &vertices);
    const VertexList &vertices() const;

    bool is_face(const Simplex &other) const;
    SimplexList get_skeleton(const Dimension dimension) const;
    Dimension dimension() const;
    void print() const;
    Simplex operator-(const Simplex &other) const;
    bool operator==(const Simplex &other) const;

private:
    void combination_util(
        const Dimension dimension,
        const uint32_t combination_index,
        SimplexList &result,
        VertexList &current_combination,
        const uint32_t array_index) const;

    VertexList vertices_;
};

struct Hash
{
    std::size_t operator()(const Simplex &simplex) const;
};

SimplexList create_simplices(const std::vector<VertexList> &simplices_in);
SimplexList select_simplices_by_dimension(const SimplexList &simplices, const Dimension dimension);
SimplexList select_higher_dimensional_simplices(const SimplexList &simplices, const Dimension dimension);
SimplexList sort_simplices(const SimplexList &simplices, const bool ascending);

#endif