#pragma once
#include "FunctionManager.h"

#include <functional>
#include <random>
#include <string>

namespace ga {

using gene = bool;
using chromosome = std::vector<gene>;
using chromosome_cit = chromosome::const_iterator; // rename?, or remove?

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
    FirstImprovementRandom,
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
                     const std::string& functionName, bool applyShift,
                     bool applyRotation);
    void sanityCheck();
    double run();
    void printBest() const; // TODO: also add stream to print to
    std::string toString() const;

  private:
    void randomizePopulationAndInitBest();

    /// Decoding chromosome and returning reference to vector to avoid
    /// creating new vector for each call.
    std::vector<double>& decodeChromosome(std::size_t index);
    /// uses decodings[index]
    std::vector<double>&
    decodeChromosome(const chromosome& chromosome, std::size_t index);
    /// Decoding version for chromosome which creates new vector
    std::vector<double> decodeChromosome(const chromosome& chromosome) const;

    /// decodes one dimension from the chromosome, given by bounds
    double
    decodeDimension(const chromosome_cit begin, const chromosome_cit end) const;

    /// convert from one encoding to another using an auxiliar
    void binaryToGray(chromosome& binary, chromosome& gray);
    void grayToBinary(chromosome& gray, chromosome& binary);
    /// convert from one encoding to another without an auxiliar
    void binaryToGray(chromosome& binary);
    void grayToBinary(chromosome& gray);
    void binaryToGreyPopulation();
    void grayToBinaryPopulation();

    /// evaluating chromosomes outside of population
    /// uses decodeChromosome 2nd overload
    double evaluateChromosome(const chromosome& chromosome);
    double evaluateChromosomeAndUpdateBest(const chromosome& chromosome);
    /// evaluating population members by index
    double evaluateChromosome(std::size_t index);
    /// evaluates chromosome outside of population, but uses decodings[index]
    /// for decoding
    double evaluateChromosome(const chromosome& chromosome, std::size_t index);
    double evaluateChromosomeAndUpdateBest(std::size_t index);

    void updateBestChromosome(double newValue, const chromosome& newBest);
    void updateBestChromosome(double newValue, std::size_t index);

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
    chromosome selectChromosome();
    /// copies in newPopulation selected chromosomes, then swaps vectors
    void selectNewPopulation();

    /// evaluates and selects the next generation
    void evaluateAndSelect();

    bool stop() const;

    /// we mutate all population except half the elites
    void mutatePopulation();
    void mutateChromosome(chromosome& chromosome);

    /// select unique chormozomes for crossover and for each of them do
    /// crossover with one (any) random chromosome (could even be itself, is
    /// this ok?)
    void crossoverPopulationChaotic();
    /// doing crossover only between unique chromosomes in a single iteration
    void crossoverPopulationClassic();
    /// sorting indices and doing crossover only for indices with value lower
    /// than crossoverProbability
    void crossoverPopulationSorted();
    void crossoverChromosomes(std::size_t i, std::size_t j);

    /// applies one iteration of hillclimbing to all population
    void hillclimbPopulation(); // TODO: test std::threads vs execution::unseq
    void hillclimbChromosome(std::size_t index);
    void hillclimbChromosome(chromosome& chromosome, std::size_t index);
    void hillclimbBest();
    void applyHillclimbing(chromosome& chromosome, std::size_t index);
    bool
    firstImprovementHillclimbing(chromosome& chromosome, std::size_t index);
    bool firstImprovementRandomHillclimbing(chromosome& chromosome,
                                            std::size_t index);
    bool bestImprovementHillclimbing(chromosome& chromosome, std::size_t index);

    /// Adaptation of hyperparameters depending on various factors
    void adapt();
    // TODO: adapt better

    /// ctor stuff
    void initContainers();
    void initStrategies(CrossoverType crossoverType,
                        HillclimbingType hillclimbingType);
    void initDistributions(
        int populationSize); // we will not need this if template Size

    void printChromosome(const chromosome& chromosome) const;
    void printChromosomeRepr(const chromosome& chromosome) const;
    void printPopulation() const;

    // TODO: test against char, might be faster because std::vector<double> is
    // space efficient but not time efficient
    std::vector<chromosome> population;
    std::vector<chromosome> newPopulation;
    std::vector<std::vector<double>> decodings;
    // used to optimize rotate operation
    std::vector<std::vector<double>> auxiliars;
    std::vector<double> fitnesses;
    std::vector<double> selectionProbabilities;
    std::vector<std::size_t> indices; // [0, ..populationSize)
    // considering population.size() == decoding.size() == fitnesses.size() ==
    // selectionProbabilities.size() == population, it might be an idea to use
    // std array

    chromosome bestChromosome;
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
    const int bitsPerChromosome; // TODO: make template bitsPerChromosome
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
    std::uniform_int_distribution<> radomChromosome; // initialized in ctor
    std::uniform_int_distribution<> randomBitIndex;  // initialized in ctor

    /// applies decoding within bounds
    std::function<long long(const chromosome_cit, const chromosome_cit)>
        decodingStrategy;
    std::function<void()> crossoverPopulationStrategy;
    std::function<bool(chromosome&, std::size_t)> hillclimbingStrategy;
    FunctionManager function;
};

GeneticAlgorithm getDefault(const std::string& functionName);

} // namespace ga
