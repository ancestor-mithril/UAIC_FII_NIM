#pragma once
#include "FunctionManager.h"

#include <functional>
#include <random>
#include <string>

namespace ga {

using chromozome = std::vector<bool>;
using chromozome_cit = chromozome::const_iterator; // rename?, or remove?

enum class CrossoverType
{
    Chaotic,
    Classic,
    Sorted,
};

enum class HillclimbingType
{
    BestImprovement,
    FirstImprovement,
    // TODO: add the others
};

// TODO: Maybe use template for population size
class GeneticAlgorithm
{
  public:
    // TODO: Change clang format to one parameter per line
    GeneticAlgorithm(double crossoverProbability, double mutationProbability,
                     double hypermutationRate, double elitesPercentage,
                     double selectionPressure, CrossoverType crossoverType,
                     HillclimbingType hillclimbingType, int populationSize,
                     int dimensions, int stepsToHypermutation,
                     int encodingChangeRate, int maxNoImprovementSteps,
                     std::string&& functionName, bool applyShift,
                     bool applyRotation);
    void sanityCheck();
    double run();
    void printBest() const; // TODO: also add stream to print to

  private:
    void randomizePopulationAndInitBest();

    /// Decoding chromozome and returning reference to vector to avoid
    /// creating new vector for each call.
    std::vector<double>& decodeChromozome(std::size_t index);
    /// uses decodings[index]
    std::vector<double>&
    decodeChromozome(const chromozome& chromozome, std::size_t index);
    /// Decoding version for chromozome which creates new vector
    std::vector<double> decodeChromozome(const chromozome& chromozome) const;

    /// decodes one dimension from the chromozome, given by bounds
    double
    decodeDimension(const chromozome_cit begin, const chromozome_cit end) const;

    /// convert from one encoding to another using an auxiliar
    void binaryToGray(chromozome& binary, chromozome& gray);
    void grayToBinary(chromozome& gray, chromozome& binary);
    /// convert from one encoding to another without an auxiliar
    void binaryToGray(chromozome& binary);
    void grayToBinary(chromozome& gray);
    void binaryToGreyPopulation();
    void grayToBinaryPopulation();

    /// evaluating chromozomes outside of population
    /// uses decodeChromozome 2nd overload
    double evaluateChromozome(const chromozome& chromozome) const;
    double evaluateChromozomeAndUpdateBest(const chromozome& chromozome);
    /// evaluating population members by index
    double evaluateChromozome(std::size_t index);
    /// evaluates chromozome outside of population, but uses decodings[index]
    /// for decoding
    double evaluateChromozome(const chromozome& chromozome, std::size_t index);
    double evaluateChromozomeAndUpdateBest(std::size_t index);

    void updateBestChromozome(double newValue, const chromozome& newBest);
    void updateBestChromozome(double newValue, std::size_t index);

    /// this has to be done sequentially
    void evaluatePopulation();
    /// only updates best
    void updateBestFromPopulation();

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

    /// evaluates and selects the next generation
    void evaluateAndSelect();

    bool stop() const;

    /// we mutate all population except half the elites
    void mutatePopulation();
    void mutateChromozome(chromozome& chromozome);

    /// select unique chormozomes for crossover and for each of them do
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
    void hillclimbChromozome(chromozome& chromozome, std::size_t index);
    void hillclimbBest();
    void applyHillclimbing(chromozome& chromozome, std::size_t index);
    /// this version of first improvement hillclimbing doesn't do any sorting on
    /// indices
    bool
    firstImprovementHillclimbing(chromozome& chromozome, std::size_t index);
    // TODO: add first improvement with random sorting, best improvement, worst
    // improvement (?)

    /// Adaptation of hyperparameters depending on various factors
    void adapt();
    // TODO: adapt better

    /// ctor stuff
    void initContainers();
    void initStrategies(CrossoverType crossoverType,
                        HillclimbingType hillclimbingType);
    void initDistributions(
        int populationSize); // we will not need this if template Size

    void printChromozome(const chromozome& chromozome) const;
    void printChromozomeRepr(const chromozome& chromozome) const;
    void printPopulation() const;

    // TODO: test against char, might be faster because std::vector<double> is
    // space efficient but not time efficient
    std::vector<chromozome> population;
    std::vector<chromozome> newPopulation;
    std::vector<std::vector<double>> decodings;
    // used to optimize rotate operation
    std::vector<std::vector<double>> auxiliars;
    std::vector<double> fitnesses;
    std::vector<double> selectionProbabilities;
    std::vector<std::size_t> indices; // [0, ..populationSize)
    // considering population.size() == decoding.size() == fitnesses.size() ==
    // selectionProbabilities.size() == population, it might be an idea to use
    // std array

    chromozome bestChromozome;
    double bestValue;

    // TODO: find good values
    // TODO: use const where we should use const
    const double crossoverProbability;
    double mutationProbability;
    double hypermutationRate;
    const double elitesPercentage;
    const double selectionPressure;

    const int maxSteps;
    const int populationSize;    // TODO: make template Size
    const int dimensions;        // TODO: use constexpr
    const int bitsPerChromozome; // TODO: make template bitsPerChromozome
    const int stepsToHypermutation;
    const int encodingChangeRate;
    const int maxNoImprovementSteps;
    const int elitesNumber;

    int epoch = 0;
    int lastImprovement = 0;

    bool isBinary = true;

    std::random_device seed;
    std::mt19937_64 gen{seed()};
    std::bernoulli_distribution randomBool;
    std::uniform_real_distribution<double> randomDouble{0.0, 1.0};
    std::uniform_int_distribution<> radomChromozome; // initialized in ctor
    std::uniform_int_distribution<> randomSlice;     // initialized in ctor

    /// applies decoding within bounds
    std::function<long long(const chromozome_cit, const chromozome_cit)>
        decodingStrategy;
    std::function<void()> crossoverPopulationStrategy;
    std::function<bool(chromozome&, std::size_t)> hillclimbingStrategy;
    FunctionManager function;
};

GeneticAlgorithm getDefault(std::string&& functionName);

} // namespace ga
