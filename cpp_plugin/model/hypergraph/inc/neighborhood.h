#ifndef _NEIGHBORHOOD_H_
#define _NEIGHBORHOOD_H_

#include "typedefs.h"

class InfiniteHypergraphModel;

class NeighborhoodPart
{
public:
    NeighborhoodPart(const Position left, const Position right);
    virtual PointList create_points(const HypergraphModel::Parameters &parameters, std::mt19937 &rng) const = 0;

protected:
    Position left;
    Position right;

private:
};

class Hyperbola : public NeighborhoodPart
{
public:
    Hyperbola(
        const Position left,
        const Position right,
        const Position position,
        const Mark transformed_mark);
    PointList create_points(const HypergraphModel::Parameters &parameters, std::mt19937 &rng) const override;
    float operator()(const Position x, const float gamma) const;

private:
    float integral(const float exponent) const;

    Position position;
    Mark transformed_mark;
};

class Center : public NeighborhoodPart
{
public:
    Center(
        const Position left,
        const Position right);
    PointList create_points(const HypergraphModel::Parameters &parameters, std::mt19937 &rng) const override;
};

#endif // _NEIGHBORHOOD_H_