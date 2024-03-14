#include "infinite_model.h"
#include "infinite_network.h"
#include "point.h"

std::vector<InfiniteNetwork> InfiniteModel::generate_networks(const uint32_t num_of_infinite_networks) const
{
    std::vector<InfiniteNetwork> result{};
    result.reserve(num_of_infinite_networks);

    for (auto network_index{0U}; network_index < num_of_infinite_networks; ++network_index)
    {
        result.push_back(generate_network());
        std::cout << "\rGenerating infinite networks: " << network_index + 1U << "/" << num_of_infinite_networks;
    }
    return result;
}

float InfiniteModel::distance(const Point &first, const Point &second) const
{
    return fabs(first.position() - second.position());
}