#pragma once
#include <functional>
#include <string>
#include <vector>

namespace ga {

class FunctionManager
{
  public:
    FunctionManager(std::string&& functionName, int dimensions, bool shiftFlag,
                    bool rotateFlag);
    double operator()(std::vector<double>& x, std::vector<double>& aux) const;

  private:
    std::function<double(std::vector<double>&, std::vector<double>&)>
    initFunction(int dimensions, bool shiftFlag, bool rotateFlag);

    std::string functionName;
    std::function<double(std::vector<double>&, std::vector<double>&)> function;
};

} // namespace ga
