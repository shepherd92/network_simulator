#ifndef _COMBINATION_H_
#define _COMBINATION_H_

#include <bits/stdc++.h>
#include <vector>

typedef std::unordered_set<int32_t> combination_type;

struct hash_on_sum : private std::hash<int32_t>
{
    typedef std::hash<int32_t> base;
    std::size_t operator()(const combination_type &comb) const
    {
        return base::operator()(std::accumulate(comb.begin(), comb.end(), int32_t()));
    }
};

typedef std::unordered_set<combination_type, hash_on_sum> set_of_combinations;
void combination(const std::vector<int32_t> &array, const uint32_t k, set_of_combinations &result);

#endif