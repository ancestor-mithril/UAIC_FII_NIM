#include "PSO.h"

#include <execution>
#include <iostream>
#include <stdexcept>

namespace constants = utils::constants;

namespace pso {

namespace {

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

} // namespace

PSO getDefault(std::string_view functionName, int dimensions)
{
    return PSO{
        {{}, {}},
        functionName,
        dimensions,
        cacheStrategy::Nearest, // cacheRetrievalStrategy
        true,                   // shiftFlag
        true                    // rotateFlag
    };
}

// clang-format off
PSO::PSO(
        std::vector<swarm::SwarmParameters> swarms,
        std::string_view functionName,
        int dimensions,
        cacheStrategy cacheRetrievalStrategy,
        bool shiftFlag,
        bool rotateFlag
    )
    : functionManager{
        functionName, 
        dimensions, 
        cacheRetrievalStrategy, 
        shiftFlag, 
        rotateFlag}
// clang-format on
{
    // TODO: read from file
    std::random_device rd;

    for (auto& swarm : swarms) {
        populations.push_back(
            swarm::Swarm{dimensions, swarm, rd, functionManager});
    }
}

bool PSO::stop() const
{
    return globalBestEval <= constants::best;
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
        // std::cout << "Error: " << err.what() << std::endl;
        // max function calls reached
    }
    // std::cout << "Epochs done: " << currentEpoch << std::endl;
    //           << functionManager.getMinimum() << std::endl;
    // std::cout << "Cache hits: " << getCacheHits() << std::endl;
    return globalBestEval;
}

void PSO::runInternal()
{
    while (not stop()) {
        retrieveBestAmongSwarms();

        // TODO: refactor this to vectorize
        for (auto& swarm : populations) {
            swarm.updatePopulation(globalBest);
        }

        ++currentEpoch;
    }
}

void PSO::retrieveBestAmongSwarms()
{
    for (const auto& swarm : populations) {
        const auto swarmEvaluation = swarm.getBestEvaluation();
        if (swarmEvaluation < globalBestEval) {
            globalBestEval = swarmEvaluation;
            globalBest = swarm.getBestParticle();
        }
    }
}

} // namespace pso
