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

std::string vecToString(const std::vector<double>& v)
{
    using namespace std::string_literals;
    if (v.empty()) {
        return "[]"s;
    }
    auto ret = "["s + std::to_string(v[0]);
    return std::accumulate(std::next(v.begin()), v.end(), ret,
                           [](auto&& f, auto x) {
                               return std::move(f) + ","s + std::to_string(x);
                           }) +
           "]"s;
}

} // namespace

PSO getDefault(std::string_view functionName, int dimensions)
{
    return PSO{
        functionName,
        dimensions,
        100,   // populationSize
        100,   // resetThreshold
        0.3,   // inertia
        1,     // cognition
        3,     // social
        0.001, // chaosCoef
        cacheStrategy::Nearest,
        true, // augment
        true, // shiftFlag
        true  // rotateFlag
    };
}

// clang-format off
PSO::PSO(std::string_view functionName,
        int dimensions,
        int populationSize,
        int resetThreshold,
        double inertia,
        double cognition,
        double social,
        double chaosCoef,
        cacheStrategy cacheRetrievalStrategy,
        bool augment,
        bool shiftFlag,
        bool rotateFlag)
    : functionManager{functionName, dimensions, cacheRetrievalStrategy, shiftFlag, rotateFlag}
    , dimensions{dimensions}
    , resetThreshold{resetThreshold}
    , populationSize{populationSize}
    , inertia{inertia}
    , cognition{cognition}
    , social{social}
    , chaosCoef{chaosCoef}
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
    populationInertia = std::vector<double>(populationSize);
    evaluations = std::vector<double>(populationSize);
    populationPastBestEval = std::vector<double>(
        populationSize, std::numeric_limits<double>::infinity());
    globalBest.resize(dimensions);

    resetPopulation();
}

void PSO::resetPopulation()
{
    for (auto i = 0; i < populationSize; ++i) {
        randomizeVector(population[i], randomFromDomain, gen);
        randomizeVector(populationVelocity[i], randomFromDomainRange, gen);
        const auto particleValue = functionManager(population[i], aux[i]);

        if (particleValue < populationPastBestEval[i]) {
            populationPastBests[i] = population[i];
            populationPastBestEval[i] = particleValue;
        }

        if (particleValue < globalBestEval) {
            globalBestEval = particleValue;
            globalBest = population[i];
        }

        populationInertia[i] = inertia;
    }
}

bool PSO::stop() const
{
    return globalBestEval <= constants::best or currentEpoch > 10000;
}

int PSO::getCacheHits() const
{
    return functionManager.hitCount();
}

std::string PSO::getBestVector() const
{
    return vecToString(globalBest);
}

double PSO::run()
{
    try {
        runInternal();
    } catch (const std::out_of_range& err) {
        // max function calls reached
    }
    // std::cout << "Epochs done: " << currentEpoch << std::endl;
    //           << functionManager.getMinimum() << std::endl;
    // std::cout << "Cache hits: " << getCacheHits() << std::endl;
    return globalBestEval;
}

void PSO::runInternal()
{
    // TODO: use a better stopping criterion
    while (not stop()) {
        // if (currentEpoch % 200 == 0) {
        //     std::cout << currentEpoch << std::endl;
        // }

        if (lastImprovement > resetThreshold) {
            std::cout << "Reset at epoch: " << currentEpoch << std::endl;
            resetPopulation();
        }

        // TODO: do updateVelocity and update best in separate loop
        // parallelize update velocity
        updateVelocity();
        evaluate();
        updateBest();
        updateInertia();
        updateInertia();
        // std::cout << currentEpoch << ' ' << functionManager.hitCount() << ' '
        // << functionManager.getEpsilon() << std::endl;
        ++currentEpoch;
        ++lastImprovement; // if we had improvement, it was already reset to 0,
                           // now it's 1
    }
}

void PSO::updateVelocity()
{
    for (auto i = 0; i < populationSize; ++i) {
        const auto rp = randomDouble(gen);
        const auto rg = randomDouble(gen);

        for (auto d = 0; d < dimensions; ++d) {
            if (augment and randomDouble(gen) < chaosCoef) {
                populationVelocity[i][d] = randomFromDomainRange(gen);
            } else {
                populationVelocity[i][d] =
                    populationInertia[i] * populationVelocity[i][d] +
                    cognition * rp *
                        (populationPastBests[i][d] - population[i][d]) +
                    social * rg * (globalBest[d] - population[i][d]);
            }

            if (populationVelocity[i][d] > constants::valuesRange) {
                populationVelocity[i][d] = constants::valuesRange;
            } else if (populationVelocity[i][d] < -constants::valuesRange) {
                populationVelocity[i][d] = -constants::valuesRange;
            }

            population[i][d] += populationVelocity[i][d];

            // TODO: Add strategy (clipping to domain or reflection)

            // TODO: See if modulo arithmetic can be used in this case.
            // How would this work: use an usigned to represent [minimum,
            // maximum] and do operations for unsigneds then convert to double
            while (population[i][d] < constants::minimum or
                   population[i][d] > constants::maximum) {
                if (population[i][d] < constants::minimum) {
                    population[i][d] =
                        2 * constants::minimum - population[i][d];
                }
                if (population[i][d] > constants::maximum) {
                    population[i][d] =
                        2 * constants::maximum - population[i][d];
                }
            }
        }
    }
}

void PSO::evaluate()
{
    std::transform(population.begin(), population.end(), aux.begin(),
                   evaluations.begin(),
                   [this](const auto& particle, auto& aux) {
                       return functionManager(particle, aux);
                   });
}

void PSO::updateBest()
{
    auto min = std::numeric_limits<double>::infinity();
    for (auto i = 0; i < populationSize; ++i) {
        if (evaluations[i] < populationPastBestEval[i]) {
            populationPastBestEval[i] = evaluations[i];
            populationPastBests[i] = population[i];

            if (evaluations[i] < min) {
                min = evaluations[i];

                // std::cout << functionManager.getFunctionName()
                //           << " Epoch: " << currentEpoch << " BEST: " <<
                //           current
                //           << '\n';
                globalBestEval = evaluations[i];
                globalBest = population[i];

                lastImprovement = 0;
            }
        }
    }
}

void PSO::updateInertia()
{
    std::transform(
        evaluations.begin(), evaluations.end(), populationInertia.begin(),
        [this, globalBestFitness = 1.0 / globalBestEval](auto evaluation) {
            const auto particleFitness = 1.0 / evaluation;
            return 1 - (inertia + randomDouble(gen) * particleFitness) /
                           (globalBestFitness + 0.1);
        });
}

} // namespace pso
