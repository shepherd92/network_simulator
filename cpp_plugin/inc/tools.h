#ifndef _TOOLS_H_
#define _TOOLS_H_

#include <cmath>
#include <cstdint>
#include <string>

bool is_close(const double first, const double second);
void log_progress(
    const uint32_t current,
    const uint32_t total,
    const uint32_t condition,
    const std::string &message = "");

inline bool is_close(const double first, const double second)
{
    constexpr auto relative_difference_factor = 1e-4;
    const auto greater_magnitude = std::max(fabs(first), fabs(second));
    return (fabs(first - second) < relative_difference_factor * greater_magnitude);
}

#endif