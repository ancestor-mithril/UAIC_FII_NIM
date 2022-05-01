#pragma once

#include "../functions/FunctionManager.h"
#include "../utils/Constants.h"

#include <limits>
#include <random>
#include <string>
#include <string_view>
#include <vector>

namespace pso::swarm {

enum class topology
{
    StaticRing,
    MiniBatchRing, // TODO: Implement
    Star,
    // Random,
    // Grid,
    // Full
};

struct SwarmParameters {
    int populationSize = 100;
    int resetThreshold = 100;
    double inertia = 0.3;
    double cognition = 1.0;
    double social = 3.0;
    double swarmAttraction = 0.0;
    double chaosCoef = 0.0;
    topology topology_ = topology::Star;
    bool selection = false;
    bool jitter = false;
};

class Swarm
{
  public:
    // clang-format off
    
    Swarm(int dimensions, const SwarmParameters& parameters, std::random_device& seed, function_layer::FunctionManager& function);
    // clang-format on

    void updatePopulation(const std::vector<double>& swarmsBest);
    double getBestEvaluation() const;
    std::vector<double> getBestParticle() const;
    std::string getBestVector() const;

  private:
    void checkForPopulationReset();
    void resetPopulation();
    void updateVelocity(const std::vector<double>& swarmsBest);
    void mutate();
    void evaluate();
    void updateBest();
    void updateInertia();
    void selectNewPopulation();
    void mutatePopulation();
    void crossOver(int indexPair1, int indexPair2);
    void crossOverPopulation();
    void endIteration();

    double getStaticRingBest(std::size_t index, std::size_t dimension) const;
    double getStarBest(std::size_t index, std::size_t dimension) const;

    // TODO: use seed from file
    std::mt19937_64 gen;
    std::uniform_real_distribution<double> randomDouble{0.0, 1.0};
    std::uniform_int_distribution<int> randomInt{0, 1};
    std::uniform_int_distribution<int> randomFromDimensions;
    std::uniform_real_distribution<double> randomFromDomain{
        utils::constants::minimum, utils::constants::maximum};
    std::uniform_real_distribution<double> randomFromDomainRange{
        -utils::constants::valuesRange, utils::constants::valuesRange};

    function_layer::FunctionManager& function;
    std::function<double()> getJitter;

    double getVisibleBest(int index, int dimensions);

    const int dimensions;
    const int resetThreshold;
    const int populationSize;
    int currentEpoch = 0;
    int lastImprovement = 0;

    // const ??
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
    std::vector<std::vector<bool>> topologyChromosomes;
    std::vector<double> populationInertia;
    std::vector<double> evaluations;
    std::vector<double> populationFitness;
    std::vector<double> selectionProbability;
    std::vector<double> populationPastBestEval;
    std::vector<double> globalBest;
    std::vector<std::size_t> indices;
    std::vector<std::size_t> neighbors;
    double globalBestEval = std::numeric_limits<double>::infinity();

    const bool selection;
};

} // namespace pso::swarm
