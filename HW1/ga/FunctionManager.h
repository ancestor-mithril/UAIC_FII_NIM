#pragma once
#include <functional>
#include <string>
#include <vector>

namespace ga {

class FunctionManager
{
  public:
    FunctionManager(const std::string& functionName, int dimensions,
                    bool shiftFlag, bool rotateFlag);

    double operator()(std::vector<double>& x, std::vector<double>& aux);

    std::string toString() const;

  private:
    std::function<double(std::vector<double>&, std::vector<double>&)>
    initFunction(int dimensions, bool shiftFlag, bool rotateFlag);

    std::string functionName;
    std::function<double(std::vector<double>&, std::vector<double>&)> function;

    int functionCalls;
    std::vector<double> values;
};

} // namespace ga
