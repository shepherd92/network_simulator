#include <cassert>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <string>

#include "tools.h"

void log_progress(const uint32_t current, const uint32_t total, const uint32_t condition, const std::string &message)
{
    static std::mutex mutex{};
    if (current % condition != 0)
    {
        return;
    }
    constexpr auto window_width{80};
    const auto progress_bar_width{window_width - static_cast<int32_t>(message.size()) - 10};
    assert(progress_bar_width > 0 && "Progress bar width must be greater than 0; Message too long.");
    const auto progress{static_cast<double>(current) / total};
    const auto position{static_cast<uint32_t>(progress * progress_bar_width)};
    const auto percent{static_cast<uint32_t>(progress * 100.)};

    std::lock_guard<std::mutex> lock_guard(mutex);
    std::cout << "\r" << message << ": [" << std::string(position, '=') << '>'
              << std::string(progress_bar_width - position, ' ') << "] "
              << std::setfill(' ') << std::setw(3) << percent << "%";

    std::cout.flush();
}