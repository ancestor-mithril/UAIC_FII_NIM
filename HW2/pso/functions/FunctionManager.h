#pragma once

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
    FunctionManager(std::string_view functionName, int dimensions,
                    bool shiftFlag, bool rotateFlag);

    double operator()(std::vector<double>& x, std::vector<double>& aux);
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

  private:
    double callFunction(std::vector<double>& x, std::vector<double>& aux);
    double callFunctionAndUpdateCache(
        std::vector<double>& x, std::vector<double>& aux,
        std::map<std::vector<double>, double>::const_iterator it);

    // TODO: Reorder
    const std::string functionName;
    const int maxFes;
    const double epsilon; // Maybe this should be dinamically adjusted
    double minimum = std::numeric_limits<double>::infinity();
    int functionCalls = 0;
    int cacheHits = 0;

    std::function<double(std::vector<double>&, std::vector<double>&)> function;
};

} // namespace function_layer
