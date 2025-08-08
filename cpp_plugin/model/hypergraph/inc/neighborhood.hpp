#ifndef _NEIGHBORHOOD_HPP_
#define _NEIGHBORHOOD_HPP_

#include <optional>

#include "typedefs.hpp"

class InfiniteHypergraphModel;

class NeighborhoodPart
{
public:
    NeighborhoodPart(const Position left, const Position right);
    virtual PointList create_points(const HypergraphModel::Parameters &parameters, std::mt19937 &rng) const = 0;

    Position left() const;
    Position right() const;
    void set_left(const Position value);
    void set_right(const Position value);

protected:
private:
    Position left_;
    Position right_;
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
    std::optional<Hyperbola> intersect_domain(const Position min_, const Position max_) const;
    std::vector<Hyperbola> get_dominating_hyperbola_parts(const Hyperbola &other) const;

    float operator()(const Position x, const float gamma) const;
    float integral(const float exponent) const;
    bool is_left_tail() const;
    bool less_than(const Hyperbola &other, const Position position) const;

    Position position() const;
    Mark transformed_mark() const;

private:
    Position position_;
    Mark transformed_mark_;
};

class Center : public NeighborhoodPart
{
public:
    Center(
        const Position left,
        const Position right);
    PointList create_points(const HypergraphModel::Parameters &parameters, std::mt19937 &rng) const override;
};

#endif // _NEIGHBORHOOD_HPP_