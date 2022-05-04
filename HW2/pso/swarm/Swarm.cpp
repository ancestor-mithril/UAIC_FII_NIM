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

void randomizeVector(std::vector<bool>& v,
                     std::uniform_int_distribution<int>& dist,
                     std::mt19937_64& gen)
{
    std::generate(v.begin(), v.end(), [&]() { return dist(gen); });
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
    , swarmTopology{parameters.topology_}
    , selection{parameters.selection}
    , randomFromDimensions{0, dimensions}
// clang-format on
{
    getJitter = [&]() -> std::function<double()> {
        if (parameters.jitter) {
            return [&]() { return randomFromDomain(gen) * 0.00005; };
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
    topologyChromosomes = std::vector<std::vector<bool>>(
        populationSize, std::vector<bool>(dimensions));
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
    resetParticles();

    neighbors.resize(populationSize + 2);

    // creating neighbors for ring topology
    std::iota(neighbors.begin(), neighbors.end(), -1);
    neighbors[0] = populationSize - 1;
    neighbors[neighbors.size() - 1] = 0;

    for(auto i = 0; i < populationSize; ++i) {
        randomizeVector(topologyChromosomes[i], randomInt, gen);
    }
}

double Swarm::getVisibleBest(int index, int dimensions)
{
    // TODO : use strategy
    if (topologyChromosomes[index][dimensions] == 0) {
        return getStaticRingBest(index, dimensions);
    }
    if (topologyChromosomes[index][dimensions] == 1) {
        return getStarBest(index, dimensions);
    }
    throw std::runtime_error("Not implemented topology");
}

void Swarm::resetParticles()
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
    checkForParticlesReset();
    selectNewPopulation();
    mutateTopologies();
    //crossOverTopologies();
    mutateParticles();
    crossOverParticles();
    updateVelocity(swarmsBest);
    evaluate();
    updateBest();
    //updateInertia();
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
    std::vector<std::vector<bool>> newTopology =
        std::vector<std::vector<bool>>(populationSize,
                                         std::vector<bool>(dimensions));

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
        newTopology[i] =
            topologyChromosomes[sortedPopulation[sortedPopulation.size() - i - 1]
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
        newTopology[i] = topologyChromosomes[selected];
    }

    population = newPopulation;
    populationVelocity = newVelocity;
    topologyChromosomes = newTopology;
}

void Swarm::mutateTopologies()
{
    auto elites = 0.0 * populationSize;
    auto mutationProbability = 0.001;

    for (int i = elites / 2; i < populationSize; ++i)
	{
		for (int j = 0; j < dimensions; ++j)
		{
			if (randomDouble(gen) < mutationProbability)
			{
				topologyChromosomes[i][j] = 1 - topologyChromosomes[i][j];
			}
		}
	}
}

void Swarm::crossOverTwoTopologies(int indexPair1, int indexPair2)
{
	auto cutOff = randomFromDimensions(gen);

	for (auto i = cutOff; i < dimensions; ++i)
	{
		std::swap(topologyChromosomes[indexPair1][i], topologyChromosomes[indexPair2][i]);
	}
}

void Swarm::crossOverTopologies()
{
    auto crossOverProbability = 0.7;
    auto availablePair = false;
	auto indexPair1 = 0;
    auto indexPair2 = 0;

	for (auto i = 0; i < populationSize; ++i)
	{
		if (randomDouble(gen) < crossOverProbability)
		{
			if (availablePair)
			{
				availablePair = false;
				indexPair2 = i;
				this->crossOverTwoTopologies(indexPair1, indexPair2);
			}
			else
			{
				availablePair = true;
				indexPair1 = i;
			}
		}
	}
}

void Swarm::crossOverTwoParticles(int indexPair1, int indexPair2)
{
    auto cutOff = randomFromDimensions(gen);

	for (auto i = cutOff; i < dimensions; ++i)
	{
		std::swap(population[indexPair1][i], population[indexPair2][i]);
        std::swap(populationVelocity[indexPair1][i], populationVelocity[indexPair2][i]);
        std::swap(topologyChromosomes[indexPair1][i], topologyChromosomes[indexPair2][i]);
	}
}

void Swarm::crossOverParticles()
{
    auto crossOverProbability = 0.1;
    auto availablePair = false;
	auto indexPair1 = 0;
    auto indexPair2 = 0;

	for (auto i = 0; i < populationSize; ++i)
	{
		if (randomDouble(gen) < crossOverProbability)
		{
			if (availablePair)
			{
				availablePair = false;
				indexPair2 = i;
				this->crossOverTwoParticles(indexPair1, indexPair2);
			}
			else
			{
				availablePair = true;
				indexPair1 = i;
			}
		}
	}
}

void Swarm::checkForParticlesReset()
{
    if (lastImprovement > resetThreshold) {
        // std::cout << "Reset at epoch: " << currentEpoch << std::endl;
        resetParticles();
    }
}

void Swarm::endIteration()
{
    ++currentEpoch;
    ++lastImprovement; // if we had improvement, it was already reset to 0,
                       // now it's 1
    inertia = 0.001 * exp(1 / (1 + 7 * function.missCount() /function.getMaxFes()));
}

void Swarm::mutateParticles()
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
                    inertia * populationVelocity[i][d] +
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
                    //    return (inertia + (1.0 - (globalBestEval / evaluation)) *
                    //                          (1.0 - inertia)) *
                    //           randomDouble(gen);
                        return inertia;
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
