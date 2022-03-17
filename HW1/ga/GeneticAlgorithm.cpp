#include "GeneticAlgorithm.h"
#include "Constants.h"

#include <cmath>
#include <execution>
#include <fstream>
#include <iostream>
#include <ranges>

namespace ranges = std::ranges;
namespace exec = std::execution;
namespace cst = ga::constants;

namespace ga {

namespace {

double
decodeBinaryVariable(const chromozome_cit begin, const chromozome_cit end)
{
    // Nice, except the formatting
    return std::accumulate(begin, end, 0LL,
                           [](auto f, auto elem) { return f * 2 + elem; }) /
               (cst::discriminator) * (cst::maximum - cst::minimum) +
           cst::minimum;
}

} // namespace

GeneticAlgorithm getDefault(std::string&& functionName)
{
    return {0.7,   // crossoverProbability
            0.001, // mutationProbability
            0.01,  // hypermutationRate
            0.04,  // elitesPercentage
            10,    // selectionPressure
            // TODO: Try a lower selection pressure
            0.1,                                // encodingChangeRate
            CrossoverType::Chaotic,             // crossoverType
            HillclimbingType::FirstImprovement, // hillclimbingType
            cst::populationSize,                // populationSize
            10,                                 // dimensions
            10,                                 // stepsToHypermutation
            1000,                                // maxNoImprovementSteps
            std::move(functionName),
            true,  // applyShift
            true}; // applyRotation
}

GeneticAlgorithm::GeneticAlgorithm(
    double crossoverProbability, double mutationProbability,
    double hypermutationRate, double elitesPercentage, double selectionPressure,
    double encodingChangeRate, CrossoverType crossoverType,
    HillclimbingType hillclimbingType, int populationSize, int dimensions,
    int stepsToHypermutation, int maxNoImprovementSteps,
    std::string&& functionName, bool applyShift, bool applyRotation)
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
    , bitsPerChromozome{dimensions * cst::bitsPerVariable}
    , stepsToHypermutation{stepsToHypermutation}
    , maxNoImprovementSteps{maxNoImprovementSteps}
    , elitesNumber{static_cast<int>(elitesPercentage * populationSize)}
    , function{std::move(functionName), dimensions, applyShift, applyRotation}
// clang-format on
{
    std::cout << "Using " << cst::bitsPerVariable << " bits per variable\n";
    std::cout << "Using " << cst::discriminator << " discriminator\n";
    std::cout << "Using " << bitsPerChromozome << " bits per chromozome\n";

    initContainers();
    initStrategies(crossoverType, hillclimbingType);
    initDistributions(populationSize);
}

void GeneticAlgorithm::sanityCheck()
{
    std::cout << "GeneticAlgorithm::sanityCheck" << '\n';
    std::cout << evaluateChromozome(0) << '\n';
    std::vector<double> ourCheck(dimensions, 0.0);
    std::cout << function(ourCheck);
}

void GeneticAlgorithm::randomizePopulationAndInitBest()
{
    for (auto& chromozome : population) {
        for (auto bit : chromozome) {
            bit = randomBool(gen); // bit is ref
        }
    }
    bestChromozome = population[0];
    bestValue = evaluateChromozome(0);
}

std::vector<double>& GeneticAlgorithm::decodeChromozome(std::size_t index)
{
    auto it = population[index].cbegin();
    for (auto i = 0; i < dimensions; ++i) {
        const auto end = std::next(it, cst::bitsPerVariable);
        decodings[index][i] = decodingStrategy(it, end);
        it = end;
    }
    // TODO: Refactor to use std algorithm
    return decodings[index];
}

std::vector<double>
GeneticAlgorithm::decodeChromozome(const chromozome& chromozome) const
{
    std::vector<double> x;
    x.reserve(dimensions);
    auto it = chromozome.cbegin();
    for (auto i = 0; i < dimensions; ++i) {
        const auto end = std::next(it, cst::bitsPerVariable);
        x.push_back(decodingStrategy(it, end));
        it = end;
    }
    return x;
}

double GeneticAlgorithm::evaluateChromozome(const chromozome& chromozome) const
{
    auto decoded = decodeChromozome(chromozome);
    return function(decoded);
}

double
GeneticAlgorithm::evaluateChromozomeAndUpdateBest(const chromozome& chromozome)
{
    auto decoded = decodeChromozome(chromozome);
    auto ret = function(decoded);
    if (ret < bestValue) {
        bestValue = ret;
        bestChromozome = chromozome;
        lastImprovement = epoch;
    }
    return ret;
}

double GeneticAlgorithm::evaluateChromozome(std::size_t index)
{
    return function(decodeChromozome(index));
}

double GeneticAlgorithm::evaluateChromozomeAndUpdateBest(std::size_t index)
{
    // might be better to use iterator instead of index
    auto ret = evaluateChromozome(index);
    if (ret < bestValue) {
        bestValue = ret;
        bestChromozome = population[index];
        lastImprovement = epoch;
    }
    return ret;
}

void GeneticAlgorithm::evaluatePopulation()
{
    fitnesses[0] = evaluateChromozomeAndUpdateBest(0);
    auto min = fitnesses[0];
    auto max = fitnesses[0];

    // should not be parallelized because it has sided effects (setting best,
    // min, max)
    // TODO: do apdate best with min, not for every index
    std::transform(std::next(indices.begin()), indices.end(),
                   std::next(fitnesses.begin()), [&](auto i) {
                       auto value = evaluateChromozomeAndUpdateBest(i);
                       if (value < min) {
                           min = fitnesses[i];
                       }
                       if (value > max) {
                           max = fitnesses[i];
                       }
                       return value;
                   });

    computeSelectionProbabilities(normalizeFitness(min, max));
}

void GeneticAlgorithm::updateBestFromPopulation()
{
    // TODO: just evaluate all and do update only with the best
    std::for_each(indices.begin(), indices.end(),
                  [&](auto i) { evaluateChromozomeAndUpdateBest(i); });
}

double GeneticAlgorithm::normalizeFitness(double min, double max)
{
    constexpr auto epsilon = 0.00001;
    auto total = 0.0;
    std::for_each(indices.begin(), indices.end(), [&](auto i) {
        fitnesses[i] =
            std::pow((max - fitnesses[i]) / (max - min + epsilon) + 1,
                     selectionPressure);
        total += fitnesses[i];
    });
    return total;
}

void GeneticAlgorithm::computeSelectionProbabilities(double total)
{
    // we reduce the scope of prev
    std::transform(fitnesses.begin(), fitnesses.end(),
                   selectionProbabilities.begin(),
                   [prev = 0.0, total](auto elem) mutable {
                       prev += elem / total;
                       return prev;
                   });
}

chromozome GeneticAlgorithm::selectChromozome()
{
    const auto random = randomDouble(gen);
    for (auto i = 0; i < populationSize; ++i) {
        if (random <= selectionProbabilities[i]) {
            return population[i];
        }
    }
    // this is returned by value hoping for RVO
    // TODO: test
    return population[populationSize - 1];
}

void GeneticAlgorithm::selectNewPopulation()
{

    if (elitesNumber > 0) {
        // using indices for partial sorting
        const auto elitesEnd = std::next(indices.begin(), elitesNumber);
        std::nth_element(
            indices.begin(), elitesEnd, indices.end(), [this](auto i, auto j) {
                return selectionProbabilities[i] > selectionProbabilities[j];
            });

        // moving best elitesNumber chromozomes to elites
        std::transform(indices.begin(), elitesEnd, newPopulation.begin(),
                       [this](auto i) { return population[i]; });

        // reseting indices
        std::iota(indices.begin(), indices.end(), 0);
    }
    // skipping elites number for both iterators
    std::transform(
        std::next(population.begin(), elitesNumber), population.end(),
        std::next(newPopulation.begin(), elitesNumber),
        [this]([[maybe_unused]] auto& elem) { return selectChromozome(); });
    // Swapping back
    population.swap(newPopulation);
}

bool GeneticAlgorithm::stop() const
{
    return (epoch - lastImprovement > maxNoImprovementSteps);
    // TODO: Add condition to check that global optimum has been achieved.
}

void GeneticAlgorithm::mutatePopulation()
{
    // skipping half the elites
    std::for_each(std::next(population.begin(), elitesNumber / 2),
                  population.end(), [&](auto& x) { mutateChromozome(x); });
    // This can be vectorized (std::execution::unseq), and tested if it brings
    // any benefit for such a small population. Might do in the long run.
    // Side effects due to using generator in multiple threads?
}

void GeneticAlgorithm::mutateChromozome(chromozome& chromozome)
{
    for (auto bit : chromozome) {
        if (randomDouble(gen) < mutationProbability) {
            bit.flip();
        }
    }
    // This solution would not work if using vector of char
    // TODO: Add template specialization for bool and char for these kind of methods
}

void GeneticAlgorithm::crossoverPopulationChaotic()
{
    std::for_each(indices.begin(), indices.end(), [this](auto i) {
        if (randomDouble(gen) < crossoverProbability / 2) {
            crossoverChromozomes(i, radomChromozome(gen));
        }
    });
    // should not be parallelized because it has side effects
}

void GeneticAlgorithm::crossoverPopulationClassic()
{
    // TODO
}

void GeneticAlgorithm::crossoverPopulationSorted()
{
    // TODO
}

void GeneticAlgorithm::crossoverChromozomes(std::size_t i, std::size_t j)
{
    // TODO: Maybe we could do crossover by using different operations,
    // for example xor and or nor
    // Try to add crossover strategies and use them

    const auto slicePosition = randomSlice(gen);
    for (auto k = 0; k < slicePosition; ++k) {
        population[i].swap(population[i][k], population[j][k]);
    }

    // TODO: should be a template specialization
    // Swap version for vector of char, might not work for bool because
    // iterator is not LegacyForwardIterator
    //
    // const auto end = std::next(population[i].begin(), randomSlice(gen));
    // std::swap_ranges(population[i].begin(), end, population[j].begin());
    //
    // this can and might be worth vectorizing for char, because char has
    // LegacyContiguousIterator and each chromozome has 350 bits
}

void GeneticAlgorithm::hillclimbPopulation()
{
    // unsequential execution
    std::for_each(
        exec::unseq, population.begin(), population.end(),
        [this](auto& chromozome) { hillclimbChromozome(chromozome); });
        // TODO: test all execution contexts to check if there is any improvement
}

void GeneticAlgorithm::hillclimbChromozome(std::size_t index)
{
    // if not used, remove
    hillclimbChromozome(population[index]);
}

void GeneticAlgorithm::hillclimbBest()
{
    // TODO: here we would actually need to change encoding until no possible improvement can be done
    hillclimbChromozome(bestChromozome);
    evaluateChromozomeAndUpdateBest(bestChromozome);
}

void GeneticAlgorithm::hillclimbChromozome(chromozome& chromozome)
{
    auto best = std::move(chromozome); // moving from chromozome
    applyHillclimbing(best);           // doing complex operations in here
    chromozome = std::move(best);      // moving back to chromozome
}

void GeneticAlgorithm::applyHillclimbing(chromozome& chromozome) const
{
    for (auto progress = true; progress;) {
        progress = hillclimbingStrategy(chromozome);
    }
}

bool GeneticAlgorithm::firstImprovementHillclimbing(
    chromozome& chromozome) const
{
    // we create a new vector with this method, instead of using the already
    // created decodings positions
    // it may be a good idea to provide a new oveload which takes a chromozome
    // and an index and modifies position[index] and then return reference to it
    auto bestValue = evaluateChromozome(chromozome);

    // this is a std::_Bit_iterator::reference
    for (auto bit : chromozome) {
        bit.flip();
        const auto value = evaluateChromozome(chromozome);
        if (value < bestValue) {
            return true;
            // returning before flipping back
        }
        bit.flip();
    }
    return false; // no improvement could be done
}

void GeneticAlgorithm::adapt()
{
    // hypermutation
    if (epoch % stepsToHypermutation == 0 or epoch % stepsToHypermutation == 1) {
        std::swap(hypermutationRate, mutationProbability);
    }

    // TODO: add encoding change
    // changing encodings is good to go over hamming walls
}

void GeneticAlgorithm::printBest() const
{
    std::cout << "Best: " << bestValue << '\n';
    printChromozome(bestChromozome);
}

void GeneticAlgorithm::printChromozome(const chromozome& chromozome) const
{
    const auto decoded = decodeChromozome(chromozome);
    for (const auto x : decoded) {
        std::cout << x << ' ';
    }
    std::cout << '\n';
}

void GeneticAlgorithm::printChromozomeRepr(const chromozome& chromozome) const
{
    const auto decoded = decodeChromozome(chromozome);
    int i = 0;
    for (const auto bit : chromozome) {
        if (i++ % cst::bitsPerVariable == 0) {
            std::cout << '\n';
            std::cout << decoded[(i - 1) / cst::bitsPerVariable] << '\n';
        }
        std::cout << bit;
    }
    std::cout << '\n';
}

void GeneticAlgorithm::printPopulation() const
{
    std::for_each(
        population.begin(), population.end(),
        [this](const auto& chromozome) { printChromozome(chromozome); });
}

void GeneticAlgorithm::run()
{
    randomizePopulationAndInitBest();
    // printPopulation();
    hillclimbPopulation();
    updateBestFromPopulation();

    for (epoch = 0; epoch < maxSteps / populationSize; ++epoch) {
        std::cout << "Epoch: " << epoch << "\tBest: " << bestValue << '\n';
        if (stop()) {
            break;
        }
        adapt();

        mutatePopulation();
        crossoverPopulationStrategy();
        evaluatePopulation();
        selectNewPopulation();
    }
    hillclimbPopulation();
    // printPopulation();
    updateBestFromPopulation();
    hillclimbBest(); // this is supposed to also change encodings to exploit all
                     // improvements in all implemented representations
    printBest();
}

void GeneticAlgorithm::initContainers()
{
    for (auto i = 0; i < populationSize; ++i) {
        population.push_back(chromozome(bitsPerChromozome, true));
        // population will be randomized at each run call
        newPopulation.push_back(chromozome(bitsPerChromozome, true));
        decodings.push_back(std::vector<double>(dimensions, 0.0));
    }

    fitnesses.resize(populationSize);
    selectionProbabilities.resize(populationSize);
    indices.resize(populationSize);
    std::iota(indices.begin(), indices.end(), 0);
}

void GeneticAlgorithm::initStrategies(CrossoverType crossoverType,
                                      HillclimbingType hillclimbingType)
{
    // TODO: add gray decoding strategy
    decodingStrategy = decodeBinaryVariable;

    crossoverPopulationStrategy = [&]() -> std::function<void()> {
        if (crossoverType == CrossoverType::Chaotic) {
            return [this]() { crossoverPopulationChaotic(); };
        }
        if (crossoverType == CrossoverType::Classic) {
            return [this]() { crossoverPopulationClassic(); };
        }
        if (crossoverType == CrossoverType::Sorted) {
            return [this]() { crossoverPopulationSorted(); };
        }
        throw std::runtime_error{"Unknown CrossoverType"};
    }();

    hillclimbingStrategy =
        [&]() -> std::function<bool(chromozome & chromozome)> {
        if (hillclimbingType == HillclimbingType::FirstImprovement) {
            return [this](chromozome& chromozome) {
                return firstImprovementHillclimbing(chromozome);
            };
        }
        throw std::runtime_error{"Implement the others"};
    }();
}

void GeneticAlgorithm::initDistributions(int populationSize)
{
    radomChromozome = std::uniform_int_distribution<>{0, populationSize - 1};
    randomSlice = std::uniform_int_distribution<>{0, bitsPerChromozome - 1};
}

} // namespace ga
