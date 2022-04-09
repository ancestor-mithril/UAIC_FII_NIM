#pragma once

#include <functional>
#include <string>
#include <string_view>
#include <vector>

namespace function_layer {

class FunctionManager
{
  public:
    FunctionManager(std::string_view functionName, int dimensions,
                    bool shiftFlag, bool rotateFlag);

  private:
    std::function<double(std::vector<double>&, std::vector<double>&)> function;

    const std::string functionName;
    const int dimensions; // TODO: remove if not needed
    const int maxFes;
    int functionCalls = 0;

    // TODO: Add caching layer
};

} // namespace function_layer
