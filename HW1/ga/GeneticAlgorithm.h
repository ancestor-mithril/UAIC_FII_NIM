#include "FunctionManager.h"

#include <functional>
#include <random>
#include <string>

namespace ga {

using const_bool_it = std::vector<bool>::const_iterator;

// TODO: Maybe use template for population size
class GeneticAlgorithm
{
  public:
    GeneticAlgorithm(double crossoverProbability, double mutationProbability,
                     double hypermutationRate, double elitesPercentage,
                     double selectionPressure, double encodingChangeRate,
                     int populationSize, int dimensions,
                     int stepsToHypermutation, int maxNoImprovementSteps,
                     std::string&& functionName);
    void sanityCheck();
    void run();

  private:
    void randomizePopulationAndInitBest();
    /// Decoding chromozome and returning reference to vector to avoid
    /// creating new vector for each call.
    std::vector<double>& decodeChromozome(std::size_t index);
    ///
    double evaluateChromozome(std::size_t index);
    double evaluateChromozomeAndUpdateBest(std::size_t index);
    bool stop() const;

    // const?
    std::random_device seed;
    std::mt19937_64 gen{seed()};
    std::bernoulli_distribution randomBool;
    std::uniform_real_distribution<double> randomDouble{0.0, 1.0};

    // TODO: find good values
    // TODO: use const where we should use const
    double crossoverProbability = 0.3;  // can this modify?
    double mutationProbability = 0.005; // can this modify?
    double hypermutationRate = 0.1;     // can this modify?
    double elitesPercentage = 0.04;     // can this modify?
    double selectionPressure = 1.01;    // can this modify?
    double encodingChangeRate = 0.1;    // can this modify?

    const int maxSteps = 200'000;
    const int populationSize = 100; // should we always use 100?
    const int dimensions = 10;
    const int bitsPerChromozome;
    const int stepsToHypermutation = 10;
    const int maxNoImprovementSteps = 100;

    int elitesNumber = 4; // can this modify?

    int epoch = 0;
    int lastImprovement = 0;
    double bestValue;
    std::vector<bool> bestChromozome;

    // TODO: test against char, might be faster because std::vector<double> is
    // space efficient but not time efficient
    std::vector<std::vector<bool>> population;
    std::vector<std::vector<double>> decodings;
    // considering population.size() == decoding.size() == population,
    // it might be an idea to use std array

    /// Takes a std::vector<bool::const_iterator expects to iterate thorugh
    /// bitsPerVariable variables, without bound checking.
    std::function<double(const const_bool_it begin)> decodingStrategy;
    FunctionManager function;
};

GeneticAlgorithm getDefault(std::string&& functionName);

} // namespace ga
