#ifndef _ADRCM_MODEL_H_
#define _ADRCM_MODEL_H_

#include <random>

#include <pybind11/numpy.h>

#include "network.h"
#include "typedefs.h"

namespace py = pybind11;

class AdrcmModel
{
public:
    struct Parameters
    {
        Parameters(const py::array_t<double> &model_parameters_input);
        dimension max_dimension;
        uint32_t num_of_nodes;
        double alpha;
        double beta;
        double gamma;
        dimension torus_dimension;
    };

    AdrcmModel(const Parameters &parameters, const uint32_t seed);
    auto generate_finite_network() const;
    auto generate_infinite_networks(const uint32_t num_of_infinite_networks) const;

private:
    class Point
    {
    public:
        Point(const vertex_id id, const float birth_time, const float position);

        inline auto distance(const Point &other) const;
        inline auto torus_distance(const Point &other, const double torus_size) const;

        inline auto id() const;
        inline const auto &birth_time() const;
        inline const auto &position() const;

        inline auto operator<(const Point &other) const;

    private:
        vertex_id id_;
        float birth_time_;
        float position_;
    };

    std::vector<Point> create_vertices() const;
    connections AdrcmModel::generate_network_connections_default(
        const std::vector<AdrcmModel::Point> &vertices,
        const bool is_finite) const;
    auto create_finite_network(const std::vector<Point> &vertices, const connections &connections) const;

    auto create_birth_times() const;
    inline auto profile_function(const double argument) const;

    const Parameters parameters;
    mutable std::mt19937 random_number_generator;
};

#include "adrcm_model.inl"

#endif