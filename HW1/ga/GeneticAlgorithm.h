#pragma once
#include "FunctionManager.h"

#include <functional>
#include <random>
#include <string>

namespace ga {

using chromozome_cit = std::vector<bool>::const_iterator;
using chromozome = std::vector<bool>;

enum class CrossoverType
{
    Chaotic,
    Classic,
    Sorted,
};

// TODO: Maybe use template for population size
class GeneticAlgorithm
{
  public:
    GeneticAlgorithm(double crossoverProbability, double mutationProbability,
                     double hypermutationRate, double elitesPercentage,
                     double selectionPressure, double encodingChangeRate,
                     CrossoverType crossoverType, int populationSize,
                     int dimensions, int stepsToHypermutation,
                     int maxNoImprovementSteps, std::string&& functionName,
                     bool applyShift, bool applyRotation);
    void sanityCheck();
    void run();
    void printBest() const; // TODO: also add stream to print to

  private:
    void randomizePopulationAndInitBest();
    /// Decoding chromozome and returning reference to vector to avoid
    /// creating new vector for each call.
    std::vector<double>& decodeChromozome(std::size_t index);
    /// Decoding version for chromozome by reprezentation
    std::vector<double> decodeChromozome(const chromozome& chromozome) const;

    ///
    double evaluateChromozome(std::size_t index);
    double evaluateChromozomeAndUpdateBest(std::size_t index);
    /// this has to be done sequentially
    void evaluatePopulation();
    /// name is misleading, does not normalize, but computes fitness values in a
    /// single iteration and returns their sum
    double normalizeFitness(double min, double max);
    /// fitness normalization is done while computing selection probabilities,
    /// dividing by the total sum
    void computeSelectionProbabilities(double total);
    // TODO: return by value then assign, or return by reference then assign?
    chromozome selectChromozome();
    /// copies in newPopulation selected chromozomes, then swaps vectors
    void selectNewPopulation();
    bool stop() const;

    /// we mutate all population except half the elites
    void mutatePopulation();
    /// doing chromozome mutation in a separate method comes with the overhead
    /// of a function call
    void mutateChromozome(chromozome& chromozome);

    /// selection unique chormozomes for crossover and for each of them do
    /// crossover with one (any) random chromozome (could even be itself, is
    /// this ok?)
    void crossoverPopulationChaotic();
    /// doing crossover only between unique chromozomes in a single iteration
    void crossoverPopulationClassic();
    /// sorting indices and doing crossover only for indices with value lower
    /// than crossoverProbability
    void crossoverPopulationSorted();
    void crossoverChromozomes(std::size_t i, std::size_t j);

    /// applies one iteration of hillclimbing to all population
    void hillclimbPopulation(); // TODO: test std::threads vs execution::unseq
    void hillclimbChromozome(std::size_t index);
    void hillclimbChromozome(chromozome& chromozome);

    /// ctor stuff
    void initContainers();
    void initStrategies(CrossoverType crossoverType);
    void initDistributions(
        int populationSize); // we will not need this if template Size

    // const?
    std::random_device seed;
    std::mt19937_64 gen{seed()};
    std::bernoulli_distribution randomBool;
    std::uniform_real_distribution<double> randomDouble{0.0, 1.0};
    std::uniform_int_distribution<> radomChromozome; // initialized in ctor
    std::uniform_int_distribution<> randomSlice;     // initialized in ctor

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
    chromozome bestChromozome;

    // TODO: test against char, might be faster because std::vector<double> is
    // space efficient but not time efficient
    std::vector<chromozome> population;
    std::vector<chromozome> newPopulation;
    std::vector<std::vector<double>> decodings;
    std::vector<double> fitnesses;
    std::vector<double> selectionProbabilities;
    std::vector<std::size_t> indices; // [0, ..populationSize)
    // considering population.size() == decoding.size() == fitnesses.size() ==
    // selectionProbabilities.size() == population, it might be an idea to use
    // std array

    /// Takes a std::vector<bool::const_iterator and expects to iterate through
    /// bitsPerVariable variables, without bound checking.
    std::function<double(const chromozome_cit begin)> decodingStrategy;
    std::function<void()> crossoverPopulationStrategy;
    FunctionManager function;
};

GeneticAlgorithm getDefault(std::string&& functionName);

} // namespace ga
