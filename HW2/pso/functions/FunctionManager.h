#pragma once

#include <functional>
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
    int count() const
    {
        return functionCalls;
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
    int functionCalls = 0;

    std::function<double(std::vector<double>&, std::vector<double>&)> function;
};

} // namespace function_layer
