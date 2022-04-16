#include "PSO.h"

#include <iostream>
#include <stdexcept>

namespace constants = utils::constants;

namespace pso {

namespace {

// TODO: use templates and concepts
void randomizeVector(std::vector<double>& v,
                     std::uniform_real_distribution<double>& dist,
                     std::mt19937_64& gen, double l)
{
    std::generate(v.begin(), v.end(), [&, l]() { return dist(gen) * l; });
}

void randomizeVector(std::vector<double>& v,
                     std::uniform_real_distribution<double>& dist,
                     std::mt19937_64& gen)
{
    randomizeVector(v, dist, gen, 1);
}

} // namespace

PSO getDefault(std::string_view functionName, int dimensions)
{
    return PSO{
        functionName,
        dimensions,
        100,  // populationSize
        0.3,  // inertia
        1,    // cognition
        3,    // social
        true, // augment
        true, // shiftFlag
        true  // rotateFlag
    };
}

// clang-format off
PSO::PSO(std::string_view functionName,
        int dimensions,
        int populationSize,
        double inertia,
        double cognition,
        double social,
        bool augment,
        bool shiftFlag,
        bool rotateFlag)
    : functionManager{functionName, dimensions, shiftFlag, rotateFlag}
    , dimensions{dimensions}
    , populationSize{populationSize}
    , inertia{inertia}
    , cognition{cognition}
    , social{social}
    , augment{augment}
// clang-format on
{
    population = std::vector<std::vector<double>>(
        populationSize, std::vector<double>(dimensions));
    aux = std::vector<std::vector<double>>(populationSize,
                                           std::vector<double>(dimensions));
    populationVelocity = std::vector<std::vector<double>>(
        populationSize, std::vector<double>(dimensions));
    populationPastBests = std::vector<std::vector<double>>(
        populationSize, std::vector<double>(dimensions));
    populationPastBestEval.resize(populationSize);
    globalBest.resize(dimensions);

    for (auto i = 0; i < populationSize; ++i) {
        randomizeVector(population[i], randomFromDomain, gen);
        randomizeVector(populationVelocity[i], randomFromDomainRange, gen);

        populationPastBests[i] = population[i];
        populationPastBestEval[i] = functionManager(population[i], aux[i]);

        if (populationPastBestEval[i] < globalBestEval) {
            globalBestEval = populationPastBestEval[i];
            globalBest = population[i];
        }
    }
    population2 = population;
}

bool PSO::stop() const
{
    return globalBestEval <= constants::best;
}

double PSO::run()
{
    try {
        runInternal();
    } catch (const std::out_of_range& err) {
        // max function calls reached
    }
    std::cout << "Epochs done: " << currentEpoch << std::endl;
    std::cout << "Cache hits: " << functionManager.hitCount() << std::endl;
    std::cout << "Least difference higher than epsilon: "
              << functionManager.getMinimum() << std::endl;
    return globalBestEval;
}

void PSO::runInternal()
{
    while (not stop()) {

        for (auto i = 0; i < populationSize; ++i) {
            updateVelocity(i);
            updateBest(i);
        }

        ++currentEpoch;
        ++lastImprovement; // if we had improvement, it was already reset to 0,
                           // now it's 1
    }
}

void PSO::updateBest(int i)
{
    std::copy(population[i].begin(), population[i].end(),
              population2[i].begin());
    const auto current = functionManager(
        population2[i], aux[i]); // calling function with copied value
    if (current < populationPastBestEval[i]) {
        populationPastBestEval[i] = current;
        populationPastBests[i] = population[i];

        if (current < globalBestEval) {
            // std::cout << "BEST: " << current << '\n';
            globalBestEval = current;
            globalBest = population[i];

            lastImprovement = 0;
        }
    }
}

void PSO::updateVelocity(int i)
{
    const auto rp = randomDouble(gen);
    const auto rg = randomDouble(gen);

    for (auto d = 0; d < dimensions; ++d) {
        // TODO: use parameter for 0.001
        if (augment and randomDouble(gen) < 0.001) {
            populationVelocity[i][d] = randomFromDomainRange(gen);
        }

        populationVelocity[i][d] =
            inertia * populationVelocity[i][d] +
            cognition * rp * (populationPastBests[i][d] - population[i][d]) +
            social * rg * (globalBest[d] - population[i][d]);

        if (populationVelocity[i][d] > constants::valuesRange) {
            populationVelocity[i][d] = constants::valuesRange;
        } else if (populationVelocity[i][d] < -constants::valuesRange) {
            populationVelocity[i][d] = -constants::valuesRange;
        }

        population[i][d] += populationVelocity[i][d];

        // TODO: Add strategy (clipping to domain or reflection)

        // TODO: See if modulo arithmetic can be used in this case.
        // How would this work: use an usigned to represent [minimum, maximum]
        // and do operations for unsigneds then convert to double
        while (population[i][d] < constants::minimum or
               population[i][d] > constants::maximum) {
            if (population[i][d] < constants::minimum) {
                population[i][d] = 2 * constants::minimum - population[i][d];
            }
            if (population[i][d] > constants::maximum) {
                population[i][d] = 2 * constants::maximum - population[i][d];
            }
        }
    }
}

} // namespace pso
