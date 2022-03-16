
#include "Cec22.h"
#include "GeneticAlgorithm.h"

#include <iostream>

void test();
int main()
{
    test();
    return 0;
}

void test()
{
    ga::functions::sanity_check();
    auto ga = ga::getDefault("rastrigin_func");
    ga.sanityCheck();
    ga.run();
}
