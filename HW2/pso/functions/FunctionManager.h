#pragma once

#include "../utils/Timer.h"
#include "CacheLayer.h"

#include <functional>
#include <limits>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace function_layer {

class FunctionManager
{
  public:
    FunctionManager();

    FunctionManager(std::string_view functionName, int dimensions,
                    cache_layer::KDTreeCache::CacheRetrievalStrategy
                        cacheRestrievalStrategy,
                    bool shiftFlag, bool rotateFlag);

    double operator()(const std::vector<double>& x, std::vector<double>& aux);
    double cheat(const std::vector<double>& x, std::vector<double>& aux);
    int missCount() const
    {
        return functionCalls;
    }
    int hitCount() const
    {
        return cacheHits;
    }
    double getMinimum() const
    {
        return minimum;
    }
    double getEpsilon() const
    {
        return epsilon;
    }
    std::string getFunctionName() const
    {
        return functionName;
    }

    static int rebalance;

  private:
    double callFunction(const std::vector<double>& x, std::vector<double>& aux);
    double callFunctionAndUpdateCache(const std::vector<double>& x,
                                      std::vector<double>& aux);

    // TODO: Reorder
    const std::string functionName;
    const int maxFes = 200'000;
    double epsilon = 0.1;
    double minimum = std::numeric_limits<double>::infinity();
    int functionCalls = 0;
    int cacheHits = 0;

    std::function<double(const std::vector<double>&, std::vector<double>&)>
        function;
    cache_layer::KDTreeCache cache;
};

} // namespace function_layer
