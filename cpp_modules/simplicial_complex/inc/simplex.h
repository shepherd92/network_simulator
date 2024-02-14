#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#include <bits/stdc++.h>

class Simplex
{
public:
    Simplex(const std::vector<int32_t> &vertices);
    const std::set<int32_t> &vertices() const;

    bool is_face(const Simplex &other) const;
    uint32_t dimension() const;
    std::vector<int32_t> operator-(const Simplex &other) const;

private:
    std::set<int32_t> vertices_;
};

std::vector<Simplex> create_simplices(const std::vector<std::vector<int32_t>> &simplices_in);
std::vector<Simplex> select_simplices_by_dimension(
    const std::vector<Simplex> &simplices,
    const uint32_t dimension);
std::vector<Simplex> select_higher_dimensional_simplices(
    const std::vector<Simplex> &simplices,
    const uint32_t dimension);
std::vector<Simplex> sort_simplices(
    const std::vector<Simplex> &simplices,
    const bool ascending);

#endif