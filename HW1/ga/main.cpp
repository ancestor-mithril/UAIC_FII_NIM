
#include "Cec22.h"
#include "GeneticAlgorithm.h"

#include <execution>
#include <fstream>
#include <iostream>

void test();
void runExperiments();
int main()
{
    // runExperiments();
    test();
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
    auto mutations = {// 0.0001,
                      0.001, 0.005, 0.01, 0.05};
    auto elites = {
        // 0.0,
        0.02,
        0.04,
        0.08,
    };
    auto selectionPressure = {0.80, 1.0, 1.2, 5.0, 10.0};
    auto crossoverTypes = {
        // ga::CrossoverType::Chaotic,
        ga::CrossoverType::Classic,
        // ga::CrossoverType::Sorted,
    };
    auto hillclimbings = {
        ga::HillclimbingType::BestImprovement,
        // ga::HillclimbingType::FirstImprovement,
        // ga::HillclimbingType::FirstImprovementRandom,
    };
    auto hyperMutationSteps = {
        5,
        10,
        20,
    };
    auto encodingChanges = {// 1,
                            2, 5, 10, 20};

    auto i = 0;
    for (auto c : crossovers) {
        for (auto m : mutations) {
            for (auto e : elites) {
                for (auto s : selectionPressure) {
                    for (auto ct : crossoverTypes) {
                        for (auto hc : hillclimbings) {
                            for (auto hs : hyperMutationSteps) {
                                // why doesn't parallelization work
                                std::for_each(
                                    std::execution::unseq,
                                    encodingChanges.begin(),
                                    encodingChanges.end(), [&](auto ec) {
                                        auto copy = i++;
                                        std::cout << copy << '\n';
                                        std::ofstream f{"experiments/" +
                                                        std::to_string(copy) +
                                                        ".txt"};
                                        f << "Crossover : " << c << '\n';
                                        f << "Mutation : " << m << '\n';
                                        f << "hyperMutations : " << m * 10
                                          << '\n';
                                        f << "elites : " << e << '\n';
                                        f << "selectionPressure : " << s
                                          << '\n';
                                        f << "crossoverTypes : " << (int)ct
                                          << '\n';
                                        f << "hillclimbings : " << (int)hc
                                          << '\n';
                                        f << "hyperMutationSteps : " << hs
                                          << '\n';
                                        f << "encodingChanges : " << ec << '\n';
                                        for (auto exp = 0; exp < 3; ++exp) {
                                            ga::GeneticAlgorithm ga{
                                                c,
                                                m,
                                                m * 10,
                                                e,
                                                s,
                                                ct,
                                                hc,
                                                100,
                                                10,
                                                hs,
                                                ec,
                                                2000,
                                                "rastrigin_func",
                                                true,
                                                true};
                                            auto rez = ga.run();
                                            f << rez << '\n';
                                        }
                                    });
                            }
                        }
                    }
                }
            }
        }
    }
}
