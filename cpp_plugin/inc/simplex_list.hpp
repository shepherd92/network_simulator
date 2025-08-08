#ifndef _SIMPLEX_LIST_HPP_
#define _SIMPLEX_LIST_HPP_

#include <map>
#include <vector>

#include "simplex.hpp"
#include "typedefs.hpp"

class SimplexList
{
public:
    SimplexList(const ISimplexList &input = {});
    SimplexList(const std::vector<Simplex> &input);
    SimplexList(SimplexSet &input);

    const std::vector<Simplex> &simplices() const;
    uint32_t size() const;
    PointIdList vertices() const;
    SimplexList facets() const;
    SimplexList faces(const Dimension dimension) const;
    SimplexList cofaces(const Simplex &simplex) const;
    SimplexList skeleton(const Dimension max_dimension) const;

    ISimplexList raw() const;
    SimplexList filter(const PointIdList &vertices_to_keep) const;
    std::map<PointId, std::vector<int32_t>> vertex_simplex_map() const;
    std::vector<uint32_t> calc_degree_sequence(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension) const;

    std::unordered_map<Simplex, uint32_t, SimplexHash>
    calc_degree_sequence(const Dimension face_dimension) const;

    std::vector<Dimension> calc_dimension_distribution() const;
    SimplexList simplices_by_dimension(const Dimension dimension) const;
    SimplexList select_higher_dimensional_simplices(const Dimension min_dimension) const;

    typedef typename std::vector<Simplex>::iterator iterator;
    typedef typename std::vector<Simplex>::const_iterator const_iterator;
    iterator begin();
    const_iterator begin() const;
    const_iterator cbegin() const;
    iterator end();
    const_iterator end() const;
    const_iterator cend() const;

    const Simplex &operator[](const int32_t index) const;
    Simplex &operator[](const int32_t index);
    SimplexList operator+(const SimplexList &other) const;
    void operator+=(const SimplexList &other);
    bool operator<(const SimplexList &other) const;
    bool operator==(const SimplexList &other) const;

private:
    void sort(const bool ascending);
    std::vector<Simplex> simplices_;
};

#endif