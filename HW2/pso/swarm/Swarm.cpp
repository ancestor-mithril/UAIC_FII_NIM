#include "Swarm.h"

#include <execution>
#include <iostream>
#include <stdexcept>

namespace constants = utils::constants;

namespace pso::swarm {

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
                           [](auto&& f, const auto x) {
                               return std::move(f) + ","s + std::to_string(x);
                           }) +
           "]"s;
}

constexpr auto epsilon = 1e-6;

} // namespace

// clang-format off
Swarm::Swarm(
        int dimensions,
        const SwarmParameters& parameters,
        std::random_device& seed, 
        function_layer::FunctionManager& function)
    : gen{std::mt19937_64(seed())}
    , function{function}
    , dimensions{dimensions}
    , resetThreshold{parameters.resetThreshold}
    , populationSize{parameters.populationSize}
    , inertia{parameters.inertia}
    , cognition{parameters.cognition}
    , social{parameters.social}
    , swarmAttraction{parameters.swarmAttraction}
    , chaosCoef{parameters.chaosCoef}
    , swarmTopology{parameters.topology}
    , selection{parameters.selection}
// clang-format on
{
    getJitter = [&]() -> std::function<double()> {
        if (parameters.jitter) {
            return [&]() { return randomFromDomain(gen) * 0.0001; };
        }
        return []() { return 0.0; };
    }();

    population = std::vector<std::vector<double>>(
        populationSize, std::vector<double>(dimensions));

    // Maybe only one aux for every swarm is enough
    aux = std::vector<std::vector<double>>(populationSize,
                                           std::vector<double>(dimensions));
    populationVelocity = std::vector<std::vector<double>>(
        populationSize, std::vector<double>(dimensions));
    populationPastBests = std::vector<std::vector<double>>(
        populationSize, std::vector<double>(dimensions));
    populationInertia = std::vector<double>(populationSize);
    evaluations = std::vector<double>(populationSize);
    populationFitness = std::vector<double>(populationSize);
    selectionProbability = std::vector<double>(populationSize);
    populationPastBestEval = std::vector<double>(
        populationSize, std::numeric_limits<double>::infinity());
    globalBest.resize(dimensions);

    indices.resize(populationSize);
    std::iota(indices.begin(), indices.end(), 0);
    // indices are used in here
    resetPopulation();

    neighbors.resize(populationSize + 2);

    // creating neighbors for ring topology
    std::iota(neighbors.begin(), neighbors.end(), -1);
    neighbors[0] = populationSize - 1;
    neighbors[neighbors.size() - 1] = 0;
}

double Swarm::getVisibleBest(int index, int dimensions)
{
    // TODO : use strategy
    if (swarmTopology == topology::StaticRing) {
        return getStaticRingBest(index, dimensions);
    }
    if (swarmTopology == topology::Star) {
        return getStarBest(index, dimensions);
    }
    throw std::runtime_error("Not implemented topology");
}

void Swarm::resetPopulation()
{
    std::for_each(indices.begin(), indices.end(), [&](const auto i) {
        randomizeVector(population[i], randomFromDomain, gen);
        randomizeVector(populationVelocity[i], randomFromDomainRange, gen);
        const auto particleValue = function(population[i], aux[i]);
        evaluations[i] = particleValue;

        if (particleValue < populationPastBestEval[i]) {
            populationPastBests[i] = population[i];
            populationPastBestEval[i] = particleValue;
        }

        if (particleValue < globalBestEval) {
            globalBestEval = particleValue;
            globalBest = population[i];
        }

        populationInertia[i] = inertia;
    });
}

std::string Swarm::getBestVector() const
{
    return vecToString(globalBest);
}

void Swarm::updatePopulation(const std::vector<double>& swarmsBest)
{
    checkForPopulationReset();
    selectNewPopulation();
    mutate();
    updateVelocity(swarmsBest);
    evaluate();
    updateBest();
    updateInertia();
    endIteration();
}

void Swarm::selectNewPopulation()
{
    if (not selection) {
        return;
    }
    // Calculate selection probability
    const auto [minIt, maxIt] =
        std::minmax_element(evaluations.begin(), evaluations.end());
    const auto min = *minIt;
    const auto max = *maxIt;

    // TODO: use accumulate
    auto totalFitness = 0.0f;

    for (auto i = 0; i < populationSize; ++i) {
        populationFitness[i] =
            pow((max - evaluations[i]) / (max - min + epsilon) + 1, 10);
        totalFitness += populationFitness[i];
    }

    auto prevProb = 0.0f;

    for (auto i = 0; i < populationSize; ++i) {
        selectionProbability[i] =
            prevProb + (populationFitness[i] / totalFitness);
        prevProb = selectionProbability[i];
    }

    // Do selection
    auto elites = 0.5 * populationSize;
    std::vector<std::pair<double, int>> sortedPopulation;
    std::vector<std::vector<double>> newPopulation =
        std::vector<std::vector<double>>(populationSize,
                                         std::vector<double>(dimensions));
    std::vector<std::vector<double>> newVelocity =
        std::vector<std::vector<double>>(populationSize,
                                         std::vector<double>(dimensions));

    for (auto i = 0; i < populationSize; ++i) {
        sortedPopulation.push_back(std::make_pair(populationFitness[i], i));
    }

    std::sort(sortedPopulation.begin(), sortedPopulation.end());

    for (auto i = 0; i < elites; ++i) {
        newPopulation[i] =
            population[sortedPopulation[sortedPopulation.size() - i - 1]
                           .second];
        newVelocity[i] =
            populationVelocity[sortedPopulation[sortedPopulation.size() - i - 1]
                                   .second];
    }

    for (auto i = elites; i < populationSize; ++i) {
        auto r = randomDouble(gen);
        auto selected = populationSize - 1;

        for (auto j = 0; j < populationSize; ++j) {
            if (r < selectionProbability[j]) {
                selected = j;
                break;
            }
        }
        newPopulation[i] = population[selected];
        newVelocity[i] = populationVelocity[selected];
    }

    population = newPopulation;
    populationVelocity = newVelocity;
}

void Swarm::checkForPopulationReset()
{
    if (lastImprovement > resetThreshold) {
        // std::cout << "Reset at epoch: " << currentEpoch << std::endl;
        resetPopulation();
    }
}

void Swarm::endIteration()
{
    ++currentEpoch;
    ++lastImprovement; // if we had improvement, it was already reset to 0,
                       // now it's 1
}

void Swarm::mutate()
{
    // TODO: generate positions
    if (chaosCoef <= 0.0) {
        return;
    }
    std::for_each(indices.begin(), indices.end(), [&](auto i) {
        std::transform(populationVelocity[i].begin(),
                       populationVelocity[i].end(),
                       populationVelocity[i].begin(), [&](const auto x) {
                           if (randomDouble(gen) < chaosCoef) {
                               return randomFromDomainRange(gen);
                           }
                           return x;
                       });
    });
}

void Swarm::updateVelocity(const std::vector<double>& swarmsBest)
{
    // par_unseq or unseq?
    std::for_each(
        std::execution::par_unseq, indices.begin(), indices.end(),
        [&](const auto i) {
            const auto rCognition = randomDouble(gen);
            const auto rSocial = randomDouble(gen);
            const auto rInertia = randomDouble(gen);
            const auto rSwarm = randomDouble(gen);

            for (auto d = 0; d < dimensions; ++d) {
                // TODO: this can be faster if we only do the else and apply the
                // mutation outside when applying the mutation it is not
                // necessary to iterate through all particles all dimensions, we
                // can generate the positions that are going to be mutated
                populationVelocity[i][d] =
                    rInertia * populationInertia[i] * populationVelocity[i][d] +
                    cognition * rCognition *
                        (populationPastBests[i][d] - population[i][d]) +
                    social * rSocial *
                        (getVisibleBest(i, d) - population[i][d]) +
                    getJitter() +
                    swarmAttraction * rSwarm *
                        (swarmsBest[d] - population[i][d]);

                // TODO: Use modulo arithmetics
                if (populationVelocity[i][d] > constants::valuesRange) {
                    populationVelocity[i][d] = constants::valuesRange;
                } else if (populationVelocity[i][d] < -constants::valuesRange) {
                    populationVelocity[i][d] = -constants::valuesRange;
                }

                // TODO: Add strategy (clipping to domain or reflection)

                // TODO: See if modulo arithmetic can be used in this case.
                // How would this work: use an usigned to represent [minimum,
                // maximum] and do operations for unsigneds then convert to
                // double
                population[i][d] += populationVelocity[i][d];
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
        });
}

void Swarm::evaluate()
{
    // cannot be parallelized because exception triggers std::terminate
    std::transform(population.begin(), population.end(), aux.begin(),
                   evaluations.begin(), [&](const auto& particle, auto& aux) {
                       return function(particle, aux);
                   });
}

void Swarm::updateBest()
{
    std::for_each(indices.begin(), indices.end(), [&](const auto i) {
        if (evaluations[i] < populationPastBestEval[i]) {
            populationPastBestEval[i] = evaluations[i];
            populationPastBests[i] = population[i];

            // TODO: maybe do update best outside loop
            if (evaluations[i] < globalBestEval) {
                // std::cout << functionManager.getFunctionName()
                //           << " Epoch: " << currentEpoch << " BEST: " <<
                //           current
                //           << '\n';
                globalBestEval = evaluations[i];
                globalBest = population[i];

                lastImprovement = 0;
            }
        }
    });
}

void Swarm::updateInertia()
{
    std::transform(std::execution::unseq, evaluations.begin(),
                   evaluations.end(), populationInertia.begin(),
                   [&](const auto evaluation) {
                       return (inertia + (1.0 - (globalBestEval / evaluation)) *
                                             (1.0 - inertia)) *
                              randomDouble(gen);
                   });
}

double Swarm::getStarBest([[maybe_unused]] std::size_t index,
                          std::size_t dimension) const
{
    return globalBest[dimension];
}

double Swarm::getStaticRingBest(std::size_t index, std::size_t dimension) const
{
    const auto leftIndex = neighbors[index];
    const auto rightIndex = neighbors[index + 2];
    const auto leftBest = populationPastBestEval[leftIndex];
    const auto rightBest = populationPastBestEval[rightIndex];
    const auto currentBest = populationPastBestEval[index];

    if (leftBest < currentBest and leftBest < rightBest) {
        return populationPastBests[leftIndex][dimension];
    } else if (rightBest < currentBest and rightBest < leftBest) {
        return populationPastBests[rightIndex][dimension];
    }
    return populationPastBests[index][dimension];
}

double Swarm::getBestEvaluation() const
{
    return globalBestEval;
}

std::vector<double> Swarm::getBestParticle() const
{
    return globalBest;
}

} // namespace pso::swarm
