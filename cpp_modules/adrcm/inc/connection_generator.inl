#ifndef _CONNECTION_GENERATOR_INL_
#define _CONNECTION_GENERATOR_INL_

#include <cmath>

bool is_close(const double first, const double second);

inline bool is_close(const double first, const double second)
{
    constexpr auto relative_difference_factor = 1e-4;
    const auto greater_magnitude = std::max(fabs(first), fabs(second));
    return (fabs(first - second) < relative_difference_factor * greater_magnitude);
}

#endif