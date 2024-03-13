#ifndef _HYPERGRAPH_MODEL_H_
#define _HYPERGRAPH_MODEL_H_

#include <random>

#include <pybind11/numpy.h>

#include "typedefs.h"

namespace py = pybind11;

class HypergraphModel
{
public:
    struct Parameters
    {
        Parameters(const py::array_t<double> &parameters_in);
        Dimension max_dimension;
        float network_size;
        float interaction_intensity;
        float beta;
        float gamma;
        float gamma_prime;
    };
    HypergraphModel(const py::array_t<double> &parameters_in, const uint32_t seed);

protected:
    MarkList generate_marks(const size_t num_nodes, const Mark min_mark = 0.) const;

    ConnectionList generate_connections(const PointList &vertices, const PointList &interactions) const;

    Dimension max_dimension() const;
    float beta() const;
    float gamma() const;
    float gamma_prime() const;
    float lambda() const;
    float lambda_prime() const;

    bool connects(const Rectangle &vertex_rectangle, const Rectangle &interaction_rectangle) const;
    bool connects(const Point &vertex, const Point &interaction) const;

    mutable std::mt19937 random_number_generator_;

private:
    RectangleList create_rectangles(const PointList &points, const float exponent) const;
    float determine_network_size(const PointList &points) const;
    ConnectionList calc_connected_point_pairs(
        const Rectangle &vertex_rectangle,
        const Rectangle &interaction_rectangle) const;
    virtual float distance(const Point &first, const Point &second) const = 0;

    Parameters parameters_;
};

#endif