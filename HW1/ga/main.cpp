
#include "Cec22.h"
#include "GeneticAlgorithm.h"

#include <fstream>
#include <iostream>

void test();
void runExperiments();
int main()
{
    runExperiments();
    return 0;
}

void test()
{
    ga::functions::sanity_check();
    auto ga = ga::getDefault("rastrigin_func");
    ga.sanityCheck();
    std::cout << ga.run() << '\n';
}

void runExperiments()
{
    auto crossovers = {0.3, 0.5, 0.7};
    auto mutations = {0.0001, 0.001, 0.005, 0.01, 0.05};
    auto hyperMutations = {0.005, 0.01, 0.05, 0.1, 0.15};
    auto elites = {0.0, 0.02, 0.04, 0.08, 0.1};
    auto selectionPressure = {0.80, 1.0, 1.2, 2.0, 5.0, 8.0, 10.0};
    auto crossoverTypes = {
        ga::CrossoverType::Chaotic,
        ga::CrossoverType::Classic,
        ga::CrossoverType::Sorted,
    };
    auto hillclimbings = {
        ga::HillclimbingType::BestImprovement,
        ga::HillclimbingType::FirstImprovement,
        ga::HillclimbingType::FirstImprovementRandom,
    };
    auto hyperMutationSteps = {5, 10, 20};
    auto encodingChanges = {1, 2, 5, 10, 15, 20};

    auto i = 0;
    for (auto c : crossovers) {
        for (auto m : mutations) {
            for (auto h : hyperMutations) {
                for (auto e : elites) {
                    for (auto s : selectionPressure) {
                        for (auto ct : crossoverTypes) {
                            for (auto hc : hillclimbings) {
                                for (auto hs : hyperMutationSteps) {
                                    for (auto ec : encodingChanges) {
                                        std::cout << i << '\n';
                                        std::ofstream f{std::to_string(i) +
                                                        ".txt"};
                                        for (auto exp = 0; exp < 5; ++exp) {
                                            ga::GeneticAlgorithm ga{
                                                c,  m,    h,      e,    s,
                                                ct, hc,   100,    10,   hs,
                                                ec, 2000, "cf01", true, true};
                                            f << ga.run() << '\n';
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
