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
        {
            swarm::Swarm(
                    dimensions, 
                    500, 
                    150, 
                    0.3, 
                    1.0, 
                    3.0,
                    0.1, 
                    0.001, 
                    swarm::topology::Star, 
                    true
                )
        },
        functionName,
        dimensions,
        cacheStrategy::Nearest, // cacheRetrievalStrategy
        true,                   // shiftFlag
        true                    // rotateFlag
    };
}

// clang-format off
PSO::PSO(
        std::vector<swarm::Swarm> swarms,
        std::string_view functionName,
        int dimensions,
        cacheStrategy cacheRetrievalStrategy,
        bool shiftFlag,
        bool rotateFlag
    )
    : populations{swarms}
    , functionManager{
        std::make_shared<function_layer::FunctionManager>(
            functionName, 
            dimensions, 
            cacheRetrievalStrategy, 
            shiftFlag, 
            rotateFlag
        )
    }
    , dimensions{dimensions}
// clang-format on
{
    //Add the functionManager to each swarm. Refactor this in the future

    for(auto& swarm : populations) {
        swarm.initialize(functionManager);
    }

    retrieveBestAmongSwarms();
}

bool PSO::stop() const
{
    return globalBestEval <= constants::best or currentEpoch > 10000;
}

int PSO::getCacheHits() const
{
    return functionManager->hitCount();
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
    // TODO: use a better stopping criterion
    while (not stop()) {
        for(auto& swarm : populations) {
            swarm.updatePopulation(globalBest);
        }

        retrieveBestAmongSwarms();
        ++currentEpoch;
    }
}

void PSO::retrieveBestAmongSwarms()
{
    for(auto& swarm : populations) {
        double swarmEvaluation = swarm.getBestEvaluation();

        if(swarmEvaluation < globalBestEval) {
            globalBestEval = swarmEvaluation;
            globalBest = swarm.getBestParticle();
        }
    }
}

} // namespace pso
