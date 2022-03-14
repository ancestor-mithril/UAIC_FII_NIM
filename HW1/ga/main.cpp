
#include "Cec22.h"
#include "GeneticAlgorithm.h"

#include <iostream>

int main()
{
    ga::functions::sanity_check();
    auto ga = ga::getDefault("rastrigin_func");
    ga.sanityCheck();
    ga.run();
    return 0;
}
