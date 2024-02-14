#include "combination.h"

void combination_util(
    const std::vector<int32_t> &array,
    const uint32_t k,
    const uint32_t combination_index,
    set_of_combinations &result,
    std::vector<int32_t> &current_combination,
    const uint32_t array_index);

uint32_t choose(const uint32_t n, const uint32_t k);

void combination(const std::vector<int32_t> &array, const uint32_t k, set_of_combinations &result)
{
    std::vector<int32_t> current_combination(k);
    combination_util(array, k, 0U, result, current_combination, 0);
}

void combination_util(
    const std::vector<int32_t> &array,
    const uint32_t k,
    const uint32_t combination_index,
    set_of_combinations &result,
    std::vector<int32_t> &current_combination,
    const uint32_t array_index)
{
    if (combination_index == k)
    {
        // combination ready
        result.insert(combination_type(current_combination.begin(), current_combination.end()));
        return;
    }

    if (array_index >= array.size())
        return;

    current_combination[combination_index] = array[array_index];
    combination_util(array, k, combination_index + 1U, result, current_combination, array_index + 1U);
    combination_util(array, k, combination_index, result, current_combination, array_index + 1U);
}

uint32_t choose(const uint32_t n, const uint32_t k)
{
    return k != 0U ? (n * choose(n - 1U, k - 1U)) / k : 1U;
}

std::vector<std::vector<int32_t>> sort_simplices(const std::vector<std::vector<int32_t>> &simplices)
{
    auto sorted_simplices{simplices};

    // sort simplices themselves
    for (auto &simplex : sorted_simplices)
    {
        std::sort(simplex.begin(), simplex.end());
    }

    // sort list of simplices
    std::sort(sorted_simplices.begin(), sorted_simplices.end(), [](const std::vector<int32_t> &a, const std::vector<int32_t> &b)
              { return a.size() < b.size(); });

    auto unique_simplices{sorted_simplices};
    unique_simplices.erase(unique(unique_simplices.begin(), unique_simplices.end()), unique_simplices.end());

    return unique_simplices;
}
