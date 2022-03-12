#include "GeneticAlgorithm.h"

#include <cmath>
#include <iostream>

namespace ga {

namespace {

// Problem specific
constexpr auto minimum = -100.0;
constexpr auto maximum = 100.0;
constexpr auto valuesRange = maximum - minimum;
constexpr auto precision = 8;
const auto bitsPerVariable = static_cast<int>(
    std::ceil(std::log2(valuesRange * std::pow(10, precision))));
const auto discriminator = (1ll << bitsPerVariable) - 1.0;
// maybe use another name

double decodeVariable(std::vector<bool>::const_iterator begin)
{
    // Nice, except the formatting
    return std::accumulate(begin, std::next(begin, bitsPerVariable), 0LL,
                           [](auto f, auto elem) { return f * 2 + elem; }) /
               (discriminator) * (maximum - minimum) +
           minimum;
}

} // namespace

GeneticAlgorithm getDefault()
{
    return {0.3, 0.005, 0.1, 0.04, 1.01, 0.1, 200'000, 100, 10, 10, 100};
}

GeneticAlgorithm::GeneticAlgorithm(
    double crossoverProbability, double mutationProbability,
    double hypermutationRate, double elitesPercentage, double selectionPressure,
    double encodingChangeRate, int maxSteps, int populationSize, int dimensions,
    int stepsToHypermutation, int maxNoImprovementSteps)
    // clang-format off
    : crossoverProbability{crossoverProbability}
    , mutationProbability{mutationProbability}
    , hypermutationRate{hypermutationRate}
    , elitesPercentage{elitesPercentage}
    , selectionPressure{selectionPressure}
    , encodingChangeRate{encodingChangeRate}
    , maxSteps{maxSteps}
    , populationSize{populationSize}
    , dimensions{dimensions}
    , bitsPerChromozome{dimensions * bitsPerVariable}
    , stepsToHypermutation{stepsToHypermutation}
    , maxNoImprovementSteps{maxNoImprovementSteps}
    , elitesNumber{static_cast<int>(elitesPercentage * populationSize)}
// clang-format on
{
    std::cout << "Using " << bitsPerVariable << " bits per variable\n";
    std::cout << "Using " << discriminator << " discriminator\n";
    std::cout << "Using " << bitsPerChromozome << " bits per chromozome\n";
    for (auto i = 0; i < populationSize; ++i) {
        population.push_back(std::vector<bool>(bitsPerChromozome, true));
    }
    // will be randomized at each run call
}

void GeneticAlgorithm::sanityCheck() const
{
    std::cout << "GeneticAlgorithm::sanityCheck" << '\n';
}

void GeneticAlgorithm::randomizePopulation()
{
    for (auto& chromozome : population) {
        for (auto i = 0; i < bitsPerChromozome; ++i) {
            chromozome[i] = randomBool(gen);
        }
    }
}

void GeneticAlgorithm::run()
{
    randomizePopulation();
}

} // namespace ga
