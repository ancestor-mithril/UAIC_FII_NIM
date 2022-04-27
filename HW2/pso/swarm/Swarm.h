#pragma once

#include "../functions/FunctionManager.h"
#include "../utils/Constants.h"

#include <limits>
#include <random>
#include <string>
#include <string_view>
#include <vector>

namespace swarm {

enum class topology
{
    StaticRing,
    MiniBatchRing, // TODO: Implement
    Star,
    // Random,
    // Grid,
    // Full
};

class Swarm
{
  public:
    // clang-format off
    
    Swarm(
        int dimensions,
        int populationSize,
        int resetThreshold,
        double inertia,
        double cognition,
        double social,
        double swarmAttraction,
        double chaosCoef,
        topology topology,
        bool augment);
    // clang-format on

    void updatePopulation(const std::vector<double>& swarmsBest);
    double getBestEvaluation();
    const std::vector<double>& getBestParticle();
    std::string getBestVector() const;
    void initialize(std::shared_ptr<function_layer::FunctionManager> sharedFunctionManager);

  private:
    void checkForPopulationReset();
    void resetPopulation();
    void updateVelocity(const std::vector<double>& swarmsBest);
    void mutate();
    void evaluate();
    void updateBest();
    void updateInertia();
    void selectNewPopulation();
    void endIteration();

    double getStaticRingBest(std::size_t index, std::size_t dimension) const;
    double getStarBest(std::size_t index, std::size_t dimension) const;

    // TODO: use seed from file
    std::mt19937_64 gen;
    std::uniform_real_distribution<double> randomDouble{0.0, 1.0};
    std::uniform_real_distribution<double> randomFromDomain{
        utils::constants::minimum, utils::constants::maximum};
    std::uniform_real_distribution<double> randomFromDomainRange{
        -utils::constants::valuesRange, utils::constants::valuesRange};

    std::shared_ptr<function_layer::FunctionManager> functionManager;
    double getVisibleBest(int index, int dimensions);

    const int dimensions;
    const int resetThreshold;
    const int populationSize;
    int currentEpoch = 0;
    int lastImprovement = 0;

    double inertia;
    double cognition;
    double social;
    double swarmAttraction;
    double chaosCoef;

    topology swarmTopology;

    std::vector<std::vector<double>> population;
    std::vector<std::vector<double>> aux;
    std::vector<std::vector<double>> populationVelocity;
    std::vector<std::vector<double>> populationPastBests;
    std::vector<double> populationInertia;
    std::vector<double> evaluations;
    std::vector<double> populationFitness;
    std::vector<double> selectionProbability;
    std::vector<double> populationPastBestEval;
    std::vector<double> globalBest;
    std::vector<std::size_t> indices;
    std::vector<std::size_t> neighbors;
    double globalBestEval = std::numeric_limits<double>::infinity();

    const bool augment;
};

Swarm getDefault(int dimensions);

} // namespace swarm
