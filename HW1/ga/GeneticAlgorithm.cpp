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
const auto count = valuesRange * std::pow(10, precision);
const auto bitsPerVariable = static_cast<int>(std::ceil(std::log2(count)));
const auto discriminator = (1ll << bitsPerVariable) - 1.0;
// maybe use another name

double decodeBinaryVariable(const const_bool_it begin)
{
    // Nice, except the formatting
    return std::accumulate(begin, std::next(begin, bitsPerVariable), 0LL,
                           [](auto f, auto elem) { return f * 2 + elem; }) /
               (discriminator) * (maximum - minimum) +
           minimum;
}

} // namespace

GeneticAlgorithm getDefault(std::string&& functionName)
{
    return {0.3,
            0.005,
            0.1,
            0.04,
            1.01,
            0.1,
            100,
            10,
            10,
            100,
            std::move(functionName)};
}

GeneticAlgorithm::GeneticAlgorithm(
    double crossoverProbability, double mutationProbability,
    double hypermutationRate, double elitesPercentage, double selectionPressure,
    double encodingChangeRate, int populationSize, int dimensions,
    int stepsToHypermutation, int maxNoImprovementSteps,
    std::string&& functionName)
    // clang-format off
    : crossoverProbability{crossoverProbability}
    , mutationProbability{mutationProbability}
    , hypermutationRate{hypermutationRate}
    , elitesPercentage{elitesPercentage}
    , selectionPressure{selectionPressure}
    , encodingChangeRate{encodingChangeRate}
    , maxSteps{dimensions == 10 ? 200'000 : 1'000'000}
    , populationSize{populationSize}
    , dimensions{dimensions}
    , bitsPerChromozome{dimensions * bitsPerVariable}
    , stepsToHypermutation{stepsToHypermutation}
    , maxNoImprovementSteps{maxNoImprovementSteps}
    , elitesNumber{static_cast<int>(elitesPercentage * populationSize)}
    , function{std::move(functionName), dimensions}
// clang-format on
{
    std::cout << "Using " << bitsPerVariable << " bits per variable\n";
    std::cout << "Using " << discriminator << " discriminator\n";
    std::cout << "Using " << bitsPerChromozome << " bits per chromozome\n";
    for (auto i = 0; i < populationSize; ++i) {
        population.push_back(std::vector<bool>(bitsPerChromozome, true));
        // will be randomized at each run call
        decodings.push_back(std::vector<double>(dimensions, 0.0));
    }

    decodingStrategy = decodeBinaryVariable;
}

void GeneticAlgorithm::sanityCheck()
{
    std::cout << "GeneticAlgorithm::sanityCheck" << '\n';
    std::cout << evaluateChromozome(0) << '\n';
}

void GeneticAlgorithm::randomizePopulationAndInitBest()
{
    for (auto& chromozome : population) {
        for (auto i = 0; i < bitsPerChromozome; ++i) {
            chromozome[i] = randomBool(gen);
        }
    }
    bestChromozome = population[0];
    bestValue = evaluateChromozome(0);
}

std::vector<double>& GeneticAlgorithm::decodeChromozome(std::size_t index)
{
    auto it = population[index].cbegin();
    for (auto i = 0; i < dimensions; ++i) {
        decodings[index][i] = decodingStrategy(it);
        it = std::next(it, bitsPerChromozome);
    }
    // TODO: Refactor to use std algorithm
    return decodings[index];
}

double GeneticAlgorithm::evaluateChromozome(std::size_t index)
{
    return function(decodeChromozome(index));
}

double GeneticAlgorithm::evaluateChromozomeAndUpdateBest(std::size_t index)
{
    auto ret = evaluateChromozome(index);
    if (ret < bestValue) {
        bestValue = ret;
        bestChromozome = population[index];
        lastImprovement = epoch;
    }
    return ret;
}

bool GeneticAlgorithm::stop() const
{
    return (epoch - lastImprovement > maxNoImprovementSteps);
    // TODO: Add condition to check that global optimum has been achieved.
}

void GeneticAlgorithm::run()
{
    randomizePopulationAndInitBest();
    // TODO: Check what's the maximum number of steps
    for (epoch = 0; epoch < maxSteps / populationSize; ++epoch) {
        if (stop()) {
            break;
        }
    }
}

} // namespace ga
