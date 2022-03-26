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

long long
decodeBinaryVariable(const chromozome_cit begin, const chromozome_cit end)
{
    return std::accumulate(begin, end, 0LL,
                           [](auto f, auto elem) { return f * 2 + elem; });
}

long long
decodeGrayVariable(const chromozome_cit begin, const chromozome_cit end)
{
    auto prev = *begin;
    auto f = 0LL + prev;
    for (auto it = std::next(begin); it != end; ++it) {
        const auto temp = (*it ? not prev : prev);
        f = f * 2 + temp;
        prev = temp;
    }
    return f;
}

} // namespace

GeneticAlgorithm getDefault(const std::string& functionName)
{
    // TODO: Make tests
    return {0.9,                               // crossoverProbability
            0.005,                             // mutationProbability
            0.1,                               // hypermutationRate
            0.02,                              // elitesPercentage
            5.0,                               // selectionPressure
            CrossoverType::Classic,            // crossoverType
            HillclimbingType::BestImprovement, // hillclimbingType
            cst::populationSize,               // populationSize
            10,                                // dimensions
            10,                                // stepsToHypermutation
            5,                                 // encodingChangeRate
            1'000'000,                         // maxNoImprovementSteps
            functionName,
            true,  // applyShift
            true}; // applyRotation
}

GeneticAlgorithm::GeneticAlgorithm(
    double crossoverProbability, double mutationProbability,
    double hypermutationRate, double elitesPercentage, double selectionPressure,
    CrossoverType crossoverType, HillclimbingType hillclimbingType,
    int populationSize, int dimensions, int stepsToHypermutation,
    int encodingChangeRate, int maxNoImprovementSteps,
    const std::string& functionName, bool applyShift, bool applyRotation)
    // clang-format off
    : crossoverProbability{crossoverProbability}
    , mutationProbability{mutationProbability}
    , hypermutationRate{hypermutationRate}
    , elitesPercentage{elitesPercentage}
    , selectionPressure{selectionPressure}
    , maxSteps{dimensions == 10 ? 200'000 : 1'000'000}
    , populationSize{populationSize}
    , dimensions{dimensions}
    , bitsPerChromozome{dimensions * cst::bitsPerVariable}
    , stepsToHypermutation{stepsToHypermutation}
    , encodingChangeRate{encodingChangeRate}
    , maxNoImprovementSteps{maxNoImprovementSteps}
    , elitesNumber{static_cast<int>(elitesPercentage * populationSize)}
    , function{functionName, dimensions, applyShift, applyRotation}
// clang-format on
{
    // std::cout << "Using " << cst::bitsPerVariable << " bits per variable\n";
    // std::cout << "Using " << cst::discriminator << " discriminator\n";
    // std::cout << "Using " << bitsPerChromozome << " bits per chromozome\n";

    initContainers();
    initStrategies(crossoverType, hillclimbingType);
    initDistributions(populationSize);
    // TODO: Ofast might break some wrong assuptions about concurrency if done
    // wrong
}

void GeneticAlgorithm::sanityCheck()
{
    std::cout << "GeneticAlgorithm::sanityCheck" << '\n';
    population[0][2] = 0;
    population[0][3] = 1;
    auto firstVal = evaluateChromozome(0);
    std::cout << firstVal << '\n';
    binaryToGray(population[0], newPopulation[0]);

    decodingStrategy = decodeGrayVariable;
    auto secondVal = evaluateChromozome(0);
    if (secondVal != firstVal) {
        throw std::runtime_error{"Gray code conversion is not equivalent"};
    }

    grayToBinary(population[0], newPopulation[0]);
    decodingStrategy = decodeBinaryVariable;
    if (firstVal != evaluateChromozome(0)) {
        throw std::runtime_error{"Binary to gray not working"};
    }

    std::vector<double> ourCheck(dimensions, 0.0);
    std::vector<double> aux(dimensions, 0.0);
    std::cout << function(ourCheck, aux) << '\n';
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
    return decodeChromozome(population[index], index);
}

std::vector<double>&
GeneticAlgorithm::decodeChromozome(const chromozome& chromozome,
                                   std::size_t index)
{
    auto it = chromozome.cbegin();
    for (auto i = 0; i < dimensions; ++i) {
        const auto end = std::next(it, cst::bitsPerVariable);
        decodings[index][i] = decodeDimension(it, end);
        it = end;
    }
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
        x.push_back(decodeDimension(it, end));
        it = end;
    }
    return x;
}

double GeneticAlgorithm::decodeDimension(const chromozome_cit begin,
                                         const chromozome_cit end) const
{
    return decodingStrategy(begin, end) / cst::discriminator *
               (cst::maximum - cst::minimum) +
           cst::minimum;
}

void GeneticAlgorithm::binaryToGray(chromozome& binary, chromozome& gray)
{
    // test this against the method bellow

    // this has the advantage of doing std::vector::swap (3 * std::move), and
    // only once, while the one bellow uses swap_ranges which actually iterates
    // trough given range dimensions time, therefore O(chromozome.size())
    for (auto i = 0; i < dimensions; ++i) {
        const auto begin = i * cst::bitsPerVariable;
        const auto end = begin + cst::bitsPerVariable;

        gray[begin] = binary[begin];
        for (auto j = begin + 1; j < end; ++j) {
            gray[j] = (binary[j - 1] != binary[j]); // xor
        }
    }
    std::swap(binary, gray);
}

void GeneticAlgorithm::binaryToGray(chromozome& binary)
{
    std::array<gene, cst::bitsPerVariable> aux;
    for (auto i = 0; i < dimensions; ++i) {
        const auto begin = i * cst::bitsPerVariable;
        const auto end = begin + cst::bitsPerVariable;

        aux[0] = binary[begin];
        for (auto j = begin + 1; j < end; ++j) {
            aux[j - begin] = (binary[j - 1] != binary[j]);
        }
        std::swap_ranges(aux.begin(), aux.end(), binary.begin() + begin);
    }
    // doesn't seem to be used
}

void GeneticAlgorithm::grayToBinary(chromozome& gray, chromozome& binary)
{
    for (auto i = 0; i < dimensions; ++i) {
        const auto begin = i * cst::bitsPerVariable;
        const auto end = begin + cst::bitsPerVariable;

        binary[begin] = gray[begin];
        for (auto j = begin + 1; j < end; ++j) {
            binary[j] = gray[j] ? not binary[j - 1] : binary[j - 1];
        }
    }
    std::swap(binary, gray);
}

void GeneticAlgorithm::grayToBinary(chromozome& gray)
{
    std::array<gene, cst::bitsPerVariable> aux;
    for (auto i = 0; i < dimensions; ++i) {
        const auto begin = i * cst::bitsPerVariable;
        const auto end = begin + cst::bitsPerVariable;

        aux[0] = gray[begin];
        for (auto j = begin + 1; j < end; ++j) {
            const auto auxIndex = j - begin;
            aux[auxIndex] = gray[j] ? not aux[auxIndex - 1] : aux[auxIndex - 1];
        }
        std::swap_ranges(aux.begin(), aux.end(), gray.begin() + begin);
    }
}

void GeneticAlgorithm::binaryToGreyPopulation()
{
    // Test
    std::for_each(exec::unseq, indices.begin(), indices.end(),
                  [this](auto index) {
                      binaryToGray(population[index], newPopulation[index]);
                  });
}

void GeneticAlgorithm::grayToBinaryPopulation()
{
    // Test
    std::for_each(exec::unseq, indices.begin(), indices.end(),
                  [this](auto index) {
                      grayToBinary(population[index], newPopulation[index]);
                  });
}

double GeneticAlgorithm::evaluateChromozome(const chromozome& chromozome)
{
    auto decoded = decodeChromozome(chromozome);
    // copy used to cache result of rotate operation
    // doesn't bring any benefit for this method, however it does for the index
    // case
    auto aux = decoded;
    return function(decoded, aux);
}

double
GeneticAlgorithm::evaluateChromozomeAndUpdateBest(const chromozome& chromozome)
{
    auto decoded = decodeChromozome(chromozome);
    auto aux = decoded;
    // copy used to cache result of rotate operation
    auto ret = function(decoded, aux);

    if (ret < bestValue) {
        updateBestChromozome(ret, chromozome);
    }
    return ret;
}

double GeneticAlgorithm::evaluateChromozome(std::size_t index)
{
    return function(decodeChromozome(index), auxiliars[index]);
}

double GeneticAlgorithm::evaluateChromozome(const chromozome& chromozome,
                                            std::size_t index)
{
    return function(decodeChromozome(chromozome, index), auxiliars[index]);
}

double GeneticAlgorithm::evaluateChromozomeAndUpdateBest(std::size_t index)
{
    // might be better to use iterator instead of index
    auto ret = evaluateChromozome(index);
    if (ret < bestValue) {
        updateBestChromozome(ret, index);
    }
    return ret;
}

void GeneticAlgorithm::updateBestChromozome(double newValue, std::size_t index)
{
    bestValue = newValue;
    bestChromozome = population[index];
    if (not isBinary) {
        grayToBinary(bestChromozome, newPopulation[index]);
    }
    lastImprovement = epoch;
}

void GeneticAlgorithm::updateBestChromozome(double newValue,
                                            const chromozome& newBest)
{
    bestValue = newValue;
    bestChromozome = newBest;
    if (not isBinary) {
        grayToBinary(bestChromozome);
    }
    lastImprovement = epoch;
}

void GeneticAlgorithm::evaluatePopulation()
{
    fitnesses[0] = evaluateChromozome(0);
    auto minIndex = 0;
    auto min = fitnesses[0];
    auto max = fitnesses[0];

    // TODO: Check what execution context can support these side effects
    // (setting min, max, index)
    std::transform(std::next(indices.begin()), indices.end(),
                   std::next(fitnesses.begin()), [&](auto i) {
                       auto value = evaluateChromozome(i);
                       if (value < min) {
                           minIndex = i;
                           min = value;
                       }
                       if (value > max) {
                           max = value;
                       }
                       return value;
                   });

    // update best
    if (min < bestValue) {
        updateBestChromozome(min, minIndex);
    }

    computeSelectionProbabilities(normalizeFitness(min, max));
}

void GeneticAlgorithm::updateBestFromPopulation()
{
    std::for_each(indices.begin(), indices.end(),
                  [&](auto i) { evaluateChromozomeAndUpdateBest(i); });
}

double GeneticAlgorithm::normalizeFitness(double min, double max)
{
    constexpr auto epsilon = 0.00001;
    auto total = 0.0;
    // Execution context ?
    std::for_each(indices.begin(), indices.end(), [&](auto i) {
        fitnesses[i] =
            std::pow((max - fitnesses[i]) / (max - min + epsilon) + 1,
                     selectionPressure);
        total += fitnesses[i]; // read + write
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

void GeneticAlgorithm::evaluateAndSelect()
{
    evaluatePopulation();
    selectNewPopulation();
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
    // TODO: Add template specialization for bool and char for these kind of
    // methods
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
    bool pair = false;
    auto pairIndex = 0;
    for (auto i = 0; i < populationSize; ++i) {
        if (randomDouble(gen) < crossoverProbability) {
            if (pair) {
                crossoverChromozomes(i, pairIndex);
            } else {
                pairIndex = i;
            }
            pair = not pair;
        }
    }
}

void GeneticAlgorithm::crossoverPopulationSorted()
{
    std::transform(
        selectionProbabilities.begin(), selectionProbabilities.end(),
        selectionProbabilities.begin(),
        [this]([[maybe_unused]] auto elem) { return randomDouble(gen); });

    // sorting by p
    std::sort(indices.begin(), indices.end(), [this](auto a, auto b) {
        return selectionProbabilities[a] < selectionProbabilities[b];
    });
    // crossover between first x and last x chromozomes
    // should have taken pairs of 2 instead?
    for (auto i = 0; i < populationSize; ++i) {
        if (selectionProbabilities[indices[i]] > crossoverProbability / 2) {
            break;
        }
        crossoverChromozomes(indices[i], indices[populationSize - i - 1]);
    }
    // reseting indices
    std::iota(indices.begin(), indices.end(), 0);
}

void GeneticAlgorithm::crossoverChromozomes(std::size_t i, std::size_t j)
{
    // TODO: Maybe we could do crossover by using different operations,
    // for example xor and or nor
    // Try to add crossover strategies and use them

    const auto slicePosition = randomBitIndex(gen);
    for (auto k = 0; k < slicePosition; ++k) {
        population[i].swap(population[i][k], population[j][k]);
    }

    // TODO: should be a template specialization
    // Swap version for vector of char, might not work for bool because
    // iterator is not LegacyForwardIterator
    //
    // const auto end = std::next(population[i].begin(), randomBitIndex(gen));
    // std::swap_ranges(population[i].begin(), end, population[j].begin());
    //
    // this can and might be worth vectorizing for char, because char has
    // LegacyContiguousIterator and each chromozome has 350 bits
}

void GeneticAlgorithm::hillclimbPopulation()
{
    // unsequential execution
    std::for_each(exec::unseq, indices.begin(), indices.end(),
                  [this](auto index) { hillclimbChromozome(index); });
    // TODO: test all execution contexts to check if there is any improvement
}

void GeneticAlgorithm::hillclimbChromozome(std::size_t index)
{
    hillclimbChromozome(population[index], index);
}

void GeneticAlgorithm::hillclimbBest()
{
    auto previousBest = bestValue;
    isBinary = true;
    decodingStrategy = decodeBinaryVariable;
    // default encoding and values for current best

    auto isFirst = true;
    auto& best = bestChromozome;

    while (true) {
        if (isBinary) {
            decodingStrategy = decodeGrayVariable;
            binaryToGray(best, newPopulation[0]);
        } else {
            decodingStrategy = decodeBinaryVariable;
            grayToBinary(best, newPopulation[0]);
        }
        isBinary = not isBinary;

        // using 1st index because its free
        hillclimbChromozome(best, 0);

        auto decoded = decodeChromozome(best);
        auto aux = decoded;
        // copy used to cache result of rotate operation
        auto ret = function(decoded, aux);
        if (ret < bestValue) {
            bestValue = ret;
        }

        if (bestValue == previousBest and not isFirst) {
            break;
        }

        previousBest = bestValue;
        isFirst = false;
    }
}

void GeneticAlgorithm::hillclimbChromozome(chromozome& chromozome,
                                           std::size_t index)
{
    auto best = std::move(chromozome); // moving from chromozome
    applyHillclimbing(best, index);    // doing complex operations in here
    chromozome = std::move(best);      // moving back to chromozome
}

void GeneticAlgorithm::applyHillclimbing(chromozome& chromozome,
                                         std::size_t index)
{
    for (auto progress = true; progress;) {
        progress = hillclimbingStrategy(chromozome, index);
    }
}

bool GeneticAlgorithm::firstImprovementHillclimbing(chromozome& chromozome,
                                                    std::size_t index)
{
    const auto bestValue = evaluateChromozome(chromozome, index);

    // this is a std::_Bit_iterator::reference
    for (auto bit : chromozome) {
        bit.flip();
        const auto value = evaluateChromozome(chromozome, index);
        if (value < bestValue) {
            return true;
            // returning before flipping back
        }
        bit.flip();
    }
    return false; // no improvement could be done
}

bool GeneticAlgorithm::firstImprovementRandomHillclimbing(
    chromozome& chromozome, std::size_t index)
{
    // cannot use member indices to sort them because of concurrency issues
    const auto bestValue = evaluateChromozome(chromozome, index);

    // should we do more or less tries?
    for (std::size_t tries = 0; tries < chromozome.size(); ++tries) {
        const auto i = randomBitIndex(gen);
        chromozome[i].flip();

        const auto value = evaluateChromozome(chromozome, index);
        if (value < bestValue) {
            return true;
            // returning before flipping back
        }
        chromozome[i].flip();
    }
    return false;
}

bool GeneticAlgorithm::bestImprovementHillclimbing(chromozome& chromozome,
                                                   std::size_t index)
{
    // cannot use member indices because of concurrency issues
    auto bestValue = evaluateChromozome(chromozome, index);
    auto bestIndex = 0;
    auto updated = false;

    for (std::size_t i = 0; i < chromozome.size(); ++i) {
        chromozome[i].flip();

        const auto value = evaluateChromozome(chromozome, index);
        if (value < bestValue) {
            updated = true;
            bestIndex = i;
            bestValue = value;
        }

        chromozome[i].flip();
    }

    if (not updated) {
        return false;
    }
    chromozome[bestIndex].flip();
    return true;
}

void GeneticAlgorithm::adapt()
{
    // hypermutation
    if (epoch % stepsToHypermutation == 0 or
        epoch % stepsToHypermutation == 1) {
        std::swap(hypermutationRate, mutationProbability);
    }

    // changing encodings is good to go over hamming walls
    if (epoch % encodingChangeRate == 0) {
        if (isBinary) {
            decodingStrategy = decodeGrayVariable;
            binaryToGreyPopulation();
        } else {
            decodingStrategy = decodeBinaryVariable;
            grayToBinaryPopulation();
        }
        isBinary = not isBinary;
    }
}

std::string GeneticAlgorithm::toString() const
{
    return function.toString() + "Best: " + std::to_string(bestValue) + '\n';
}

void GeneticAlgorithm::printBest() const
{
    std::cout << "Best: " << bestValue << '\n';
    printChromozome(bestChromozome);
}

void GeneticAlgorithm::printChromozome(const chromozome& chromozome) const
{
    // TODO: make template for bool
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

double GeneticAlgorithm::run()
{
    randomizePopulationAndInitBest();
    // hillclimbPopulation();
    updateBestFromPopulation();

    for (epoch = 0; epoch < maxSteps / populationSize; ++epoch) {
        if (stop()) {
            break;
        }
        adapt();

        mutatePopulation();
        crossoverPopulationStrategy();
        evaluateAndSelect();
    }
    // hillclimbPopulation();
    // printPopulation();
    updateBestFromPopulation();
    // hillclimbBest();
    // printBest();
    return bestValue;
}

void GeneticAlgorithm::initContainers()
{
    for (auto i = 0; i < populationSize; ++i) {
        population.push_back(chromozome(bitsPerChromozome, true));
        // population will be randomized at each run call
        newPopulation.push_back(chromozome(bitsPerChromozome, true));
        decodings.push_back(std::vector<double>(dimensions, 0.0));
        auxiliars.push_back(std::vector<double>(dimensions, 0.0));
    }

    fitnesses.resize(populationSize);
    selectionProbabilities.resize(populationSize);
    indices.resize(populationSize);
    std::iota(indices.begin(), indices.end(), 0);
}

void GeneticAlgorithm::initStrategies(CrossoverType crossoverType,
                                      HillclimbingType hillclimbingType)
{
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
        [&]() -> std::function<bool(chromozome&, std::size_t)> {
        if (hillclimbingType == HillclimbingType::FirstImprovement) {
            return [this](chromozome& chromozome, std::size_t index) {
                return firstImprovementHillclimbing(chromozome, index);
            };
        }
        if (hillclimbingType == HillclimbingType::FirstImprovementRandom) {
            return [this](chromozome& chromozome, std::size_t index) {
                return firstImprovementRandomHillclimbing(chromozome, index);
            };
        }
        if (hillclimbingType == HillclimbingType::BestImprovement) {
            return [this](chromozome& chromozome, std::size_t index) {
                return bestImprovementHillclimbing(chromozome, index);
            };
        }
        throw std::runtime_error{"Implement the others"};
    }();
}

void GeneticAlgorithm::initDistributions(int populationSize)
{
    radomChromozome = std::uniform_int_distribution<>{0, populationSize - 1};
    randomBitIndex = std::uniform_int_distribution<>{0, bitsPerChromozome - 1};
}

} // namespace ga
