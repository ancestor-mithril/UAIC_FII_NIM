#pragma once

#include "../../KDTree/KDTree.hpp"
#include "../utils/Utils.h"

#include <algorithm>
#include <concepts>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <vector>

namespace function_layer::cache_layer {

namespace {

// TODO: maybe move somewhere else (utils?, new file with all concepts used?)

template <typename T>
concept IndexComparator = std::predicate<T, std::size_t, std::size_t> &&
    std::move_constructible<T> && std::copy_constructible<T>;
// cannot be called without a predicate (invokable that returns bool) that takes
// two std::size_t arguments
// cannot be called without a move-constructible and copy-constructible type

} // namespace

// TODO: Use more concepts to accept more types (std::is_array<T> or ...)
class KDTreeCache
{
  public:
    KDTree kdtree;
    std::vector<point_t> points;
    std::vector<double> values;

    KDTreeCache() = delete;

    KDTreeCache(int maxFES)
    {
        points.reserve(maxFES);
        values.reserve(maxFES);
    }

    void insert(const point_t& point, double value)
    {
        points.push_back(point);
        values.push_back(value);
        kdtree.insertPoint(point);
    }

    std::optional<double> retrieve(const point_t& point, double epsilon)
    {
        try {
            const auto [index, value] = kdtree.nearestIndexAndValue(point);
            if (value < epsilon) {
                return values[index];
            }
        } catch (const std::logic_error& e) {
            // root is empty
        }
        return std::nullopt;
    }

    std::optional<double> retrievePoints(const point_t& point, double epsilon,
                                         IndexComparator auto&& func)
    {
        const auto indices = kdtree.neighborhood(point, epsilon);
        if (indices.empty()) {
            return std::nullopt;
        }
        const auto bestIndex =
            std::min_element(indices.begin(), indices.end(), func);
        return values[*bestIndex];
    }

    std::optional<double>
    retrievePointsWorst(const point_t& point, double epsilon)
    {
        return retrievePoints(point, epsilon,
                              [this](const std::unsigned_integral auto a,
                                     const std::unsigned_integral auto b) {
                                  return values[a] > values[b];
                              });
    }

    std::optional<double>
    retrievePointsBest(const point_t& point, double epsilon)
    {
        return retrievePoints(point, epsilon,
                              [this](const std::unsigned_integral auto a,
                                     const std::unsigned_integral auto b) {
                                  return values[a] < values[b];
                              });
    }

    std::optional<double>
    retrieveFirstNeighbor(const point_t& point, double epsilon)
    {
        return kdtree.firstNeighbor(point, epsilon);
    }
};

} // namespace function_layer::cache_layer
