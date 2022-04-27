#pragma once

#include "../functions/FunctionManager.h"
#include "../swarm/Swarm.h"
#include "../utils/Constants.h"

#include <limits>
#include <random>
#include <string>
#include <string_view>
#include <vector>

namespace pso {

using cacheStrategy =
    function_layer::cache_layer::KDTreeCache::CacheRetrievalStrategy;

class PSO
{
  public:
    // clang-format off
    PSO(std::vector<swarm::SwarmParameters> swarms,
        std::string_view functionName,
        int dimensions,
        cacheStrategy cacheRetrievalStrategy,
        bool shiftFlag,
        bool rotateFlag);
    // clang-format on

    double run();

    int getCacheHits() const;
    std::string getBestVector() const;

  private:
    void runInternal();
    bool stop() const;
    void retrieveBestAmongSwarms();

    function_layer::FunctionManager functionManager;
    std::vector<swarm::Swarm> populations;

    int currentEpoch = 0;

    std::vector<double> globalBest;
    double globalBestEval = std::numeric_limits<double>::infinity();
};

PSO getDefault(std::string_view functionName, int dimensions);

} // namespace pso
