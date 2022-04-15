#include "cec22/Cec22.h"
#include "pso/PSO.h"

#include <iostream>

void runDefault();
int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[])
{
    // std::cout << cec22::sanity_check() << '\n';
    runDefault();
    return 0;
}

void runDefault()
{
    auto pso = pso::getDefault("zakharov_func", 10);
    auto result = pso.run();
    std::cout << "Result: " << result << std::endl;
}
