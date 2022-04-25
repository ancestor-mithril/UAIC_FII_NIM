#pragma once

#include "../../KDTree/KDTree.hpp"
#include "../utils/Timer.h"
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
    enum class CacheRetrievalStrategy
    {
        Nearest,
        BestNeigbor,
        WorstNeighbor,
        FirstNeighbor,
    };

    KDTree kdtree;
    std::vector<point_t> points;
    std::vector<double> values;
    std::function<std::optional<double>(const point_t& point, double epsilon)>
        retrievalStrategy;

    KDTreeCache() = delete;

    KDTreeCache(int maxFES, CacheRetrievalStrategy type)
    {
        points.reserve(maxFES);
        values.reserve(maxFES);

        retrievalStrategy =
            [this, type]() -> std::function<std::optional<double>(
                               const point_t& point, double epsilon)> {
            if (type == CacheRetrievalStrategy::Nearest) {
                return [this](const point_t& point, double epsilon) {
                    return retrieve(point, epsilon);
                };
            }
            if (type == CacheRetrievalStrategy::BestNeigbor) {
                return [this](const point_t& point, double epsilon) {
                    return retrievePointsBest(point, epsilon);
                };
            }
            if (type == CacheRetrievalStrategy::WorstNeighbor) {
                return [this](const point_t& point, double epsilon) {
                    return retrievePointsWorst(point, epsilon);
                };
            }
            return [this](const point_t& point, double epsilon) {
                return retrieveFirstNeighbor(point, epsilon);
            };
        }();
    }

    void recreate()
    {
        const auto timer = utils::timer::Timer{"Recreate KDTree"};
        kdtree = KDTree(points);
    }

    void insert(const point_t& point, double value)
    {
        const auto timer = utils::timer::Timer{"Insert into KDTree"};
        points.push_back(point);
        values.push_back(value);
        kdtree.insertPoint(point);
    }

    std::optional<double> retrieve(const point_t& point, double epsilon)
    {
        const auto timer = utils::timer::Timer{"KDTree nearest index"};
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
        const auto timer = utils::timer::Timer{"KDTree all neighbors"};
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
        const auto timer = utils::timer::Timer{"KDTree first neighbor"};
        return kdtree.firstNeighbor(point, epsilon);
    }
};

} // namespace function_layer::cache_layer
