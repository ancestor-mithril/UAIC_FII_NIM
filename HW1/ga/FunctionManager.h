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
    double f(std::vector<double>& x, std::vector<double>& aux) const;

    std::string toString() const;
    int count() const;

  private:
    std::function<double(std::vector<double>&, std::vector<double>&)>
    initFunction(int dimensions, bool shiftFlag, bool rotateFlag);

    std::string functionName;
    std::function<double(std::vector<double>&, std::vector<double>&)> function;

    int functionCalls = 0;
    int maxFes;
    std::vector<double> values;
};

} // namespace ga
