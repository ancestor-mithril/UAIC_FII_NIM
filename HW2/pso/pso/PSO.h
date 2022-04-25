#pragma once

#include "../functions/FunctionManager.h"
#include "../utils/Constants.h"

#include <limits>
#include <random>
#include <string>
#include <string_view>
#include <vector>

namespace pso {

using cacheStrategy =
    function_layer::cache_layer::KDTreeCache::CacheRetrievalStrategy;

enum class topology
{
    StaticRing,
    MiniBatchRing, // TODO: Implement
    Star,
    // Random,
    // Grid,
    // Full
};

class PSO
{
  public:
    // clang-format off
    PSO(std::string_view functionName,
        int dimensions,
        int populationSize,
        int resetThreshold,
        double inertia,
        double cognition,
        double social,
        double chaosCoef,
        cacheStrategy cacheRetrievalStrategy,
        topology topology,
        bool augment,
        bool shiftFlag,
        bool rotateFlag);
    // clang-format on

    double run();

    int getCacheHits() const;
    std::string getBestVector() const;

  private:
    void runInternal();
    void resetPopulation();
    void updateVelocity();
    void mutate();
    void evaluate();
    void updateBest();
    void updateInertia();
    bool stop() const;

    double getStaticRingBest(std::size_t index, std::size_t dimension) const;
    double getStarBest(std::size_t index, std::size_t dimension) const;

    // TODO: use seed from file
    std::random_device seed;
    std::mt19937_64 gen{seed()};
    std::uniform_real_distribution<double> randomDouble{0.0, 1.0};
    std::uniform_real_distribution<double> randomFromDomain{
        utils::constants::minimum, utils::constants::maximum};
    std::uniform_real_distribution<double> randomFromDomainRange{
        -utils::constants::valuesRange, utils::constants::valuesRange};

    function_layer::FunctionManager functionManager;
    std::function<double(std::size_t, std::size_t)> getVisibleBest;

    const int dimensions;
    const int resetThreshold;
    const int populationSize;
    int currentEpoch = 0;
    int lastImprovement = 0;

    double inertia;
    double cognition;
    double social;
    double chaosCoef;

    std::vector<std::vector<double>> population;
    std::vector<std::vector<double>> aux;
    std::vector<std::vector<double>> populationVelocity;
    std::vector<std::vector<double>> populationPastBests;
    std::vector<double> populationInertia;
    std::vector<double> evaluations;
    std::vector<double> populationPastBestEval;
    std::vector<double> globalBest;
    std::vector<std::size_t> indices;
    std::vector<std::size_t> neighbors;
    double globalBestEval = std::numeric_limits<double>::infinity();

    const bool augment;
};

PSO getDefault(std::string_view functionName, int dimensions);

} // namespace pso
