#include "cec22/Cec22.h"
#include "pso/PSO.h"

#include <iostream>

void runDefault();
int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[])
{
    std::cout << cec22::sanity_check() << '\n';
    runDefault();
    return 0;
}

void runFunction(std::string_view functionName, int dimensions)
{
    std::cout << std::endl;
    auto pso = pso::getDefault(functionName, dimensions);
    auto result = pso.run();
    std::cout << functionName << ' ' << result << std::endl;
}

void runDefault()
{
    runFunction("zakharov_func", 10);
    runFunction("rosenbrock_func", 10);
    runFunction("schaffer_F7_func", 10);
    runFunction("rastrigin_func", 10);
    runFunction("levy_func", 10);
    runFunction("hf01", 10);
    runFunction("hf02", 10);
    runFunction("hf03", 10);
    runFunction("cf01", 10);
    runFunction("cf02", 10);
    runFunction("cf03", 10);
    runFunction("cf04", 10);
}
