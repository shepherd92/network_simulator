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

    class Vertex
    {
    public:
        Vertex(const vertex_id id, const float birth_time, const float position);

        inline double distance(const Vertex &other) const;
        inline double torus_distance(const Vertex &other, const double torus_size) const;

        inline const vertex_id &id() const;
        inline const float &birth_time() const;
        inline const float &position() const;

        inline bool operator<(const Vertex &other) const;

    private:
        vertex_id id_;
        float birth_time_;
        float position_;
    };

    AdrcmModel(const Parameters &parameters, const uint32_t seed);
    Network generate_finite_network() const;

private:
    std::vector<Vertex> create_vertices() const;
    connections AdrcmModel::generate_network_connections_default(
        const std::vector<AdrcmModel::Vertex> &vertices,
        const bool is_finite) const;
    Network create_finite_network(const std::vector<Vertex> &vertices, const connections &connections) const;

    std::vector<float> create_birth_times() const;
    inline double profile_function(const double argument) const;

    const Parameters parameters;
    mutable std::mt19937 random_number_generator;
};

#include "adrcm_model.inl"

#endif