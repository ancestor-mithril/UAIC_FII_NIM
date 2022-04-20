#pragma once

#include "../functions/FunctionManager.h"
#include "../utils/Constants.h"

#include <limits>
#include <random>
#include <string>
#include <string_view>
#include <vector>

namespace pso {

class PSO
{
  public:
    // clang-format off
    PSO(std::string_view functionName,
        int dimensions,
        int populationSize,
        double inertia,
        double cognition,
        double social,
        double chaosCoef,
        bool augment,
        bool shiftFlag,
        bool rotateFlag);
    // clang-format on

    double run();

    int getCacheHits() const;
    std::string getBestVector() const;

  private:
    void runInternal();
    void updateVelocity(int i);
    void updateBest(int i);
    bool stop() const;

    // TODO: use seed from file
    std::random_device seed;
    std::mt19937_64 gen{seed()};
    std::uniform_real_distribution<double> randomDouble{0.0, 1.0};
    std::uniform_real_distribution<double> randomFromDomain{
        utils::constants::minimum, utils::constants::maximum};
    std::uniform_real_distribution<double> randomFromDomainRange{
        -utils::constants::valuesRange, utils::constants::valuesRange};

    function_layer::FunctionManager functionManager;

    std::vector<std::vector<double>> population;
    std::vector<std::vector<double>> aux;
    std::vector<std::vector<double>> populationVelocity;
    std::vector<std::vector<double>> populationPastBests;
    std::vector<double> populationPastBestEval;
    std::vector<double> globalBest;
    double globalBestEval = std::numeric_limits<double>::infinity();

    const int dimensions;
    const int populationSize;
    int currentEpoch = 0;
    int lastImprovement = 0;

    double inertia;
    double cognition;
    double social;
    double chaosCoef;

    const bool augment;
};

PSO getDefault(std::string_view functionName, int dimensions);

} // namespace pso
