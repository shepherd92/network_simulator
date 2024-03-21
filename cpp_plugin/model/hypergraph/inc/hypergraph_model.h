#ifndef _HYPERGRAPH_MODEL_H_
#define _HYPERGRAPH_MODEL_H_

#include "model.h"
#include "typedefs.h"

class HypergraphModel : virtual public Model
{
public:
    struct Parameters
    {
        Parameters(const std::vector<double> &parameters_in);
        Dimension max_dimension;
        float network_size;
        float interaction_intensity;
        float beta;
        float gamma;
        float gamma_prime;
    };
    HypergraphModel(const std::vector<double> &parameters_in);

protected:
    ConnectionList generate_connections(const PointList &vertices, const PointList &interactions) const;
    SimplexList create_simplices_from_connections(const ConnectionList &connections) const;

    inline Dimension max_dimension() const;
    inline float beta() const;
    inline float gamma() const;
    inline float gamma_prime() const;
    inline float lambda() const;
    inline float lambda_prime() const;

    bool connects(const Rectangle &vertex_rectangle, const Rectangle &interaction_rectangle) const;
    bool connects(const Point &vertex, const Point &interaction) const;

private:
    RectangleList create_rectangles(const PointList &points, const float exponent) const;
    ConnectionList calc_connected_point_pairs(
        const Rectangle &vertex_rectangle,
        const Rectangle &interaction_rectangle) const;

    Parameters parameters_;
};

#include "hypergraph_model.inl"

#endif