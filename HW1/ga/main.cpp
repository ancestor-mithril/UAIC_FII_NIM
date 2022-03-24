
#include "Cec22.h"
#include "GeneticAlgorithm.h"

#include <execution>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

void test();
void runExperiments();
void runExperiments1();
void runExperiments2();
int main(int argc, char** argv)
{

    if (argc < 2) {
        // runExperiments();
        test();
    } else {
        if (argv[1] == std::string{"1"}) {
            runExperiments1();
            return 0;
        } else if (argv[1] == std::string{"2"}) {
            runExperiments2();
            return 0;
        }

        std::ofstream fout{"experiments/" + std::string{argv[1]}};
        for (auto i = 1; i < 30; ++i) {
            auto ga = ga::getDefault(std::string{argv[1]});
            fout << ga.run() << '\n';
        }
    }

    return 0;
}

void runExperiments1()
{
    // Fine tunning crossovers, mutations and selection pressures
    std::map<double, double> crossovers = {
        {0.1, 0.0}, //
        {0.3, 0.0}, //
        {0.5, 0.0}, //
        {0.7, 0.0}, //
        {0.9, 0.0}, //
    };
    std::map<double, double> mutations = {
        {0.0005, 0.0}, //
        {0.001, 0.0},  //
        {0.005, 0.0},  //
        {0.01, 0.0},   //
        {0.05, 0.0},   //
        {0.1, 0.0},
    };

    std::map<double, double> selectionPressure = {
        {0.8, 0.0},  //
        {1.0, 0.0},  //
        {1.2, 0.0},  //
        {2.0, 0.0},  //
        {5.0, 0.0},  //
        {7.0, 0.0},  //
        {10.0, 0.0}, //
        {15.0, 0.0},
    };

    const std::vector<std::string> func = {
        "zakharov_func",    //
        "rosenbrock_func",  //
        "schaffer_F7_func", //
        "rastrigin_func",   //
        "levy_func",
        // "hf01", // TODO: make this faster
        // "cf01", // TODO: check this for problems + make faster
        // "cf02", // TODO: check this for problems + maka faster
    };

    auto minVal = 2e100; // TODO: intialize with first
    auto minCrossover = 0.1;
    auto minMutation = 0.0005;
    auto minSelection = 0.8;

    const auto hypermutationRate = 10.0; // * mutation
    const auto elitesPercentage = 0.04;
    const auto crossoverType = ga::CrossoverType::Classic;
    const auto hillclimbingType = ga::HillclimbingType::BestImprovement;
    const auto populationSize = 100;
    const auto dimensions = 10;
    const auto stepsToHypermutation = 10;
    const auto encodingChangeRate = 5;
    const auto maxNoImprovementSteps = 1900;
    const auto applyShift = true;
    const auto applyRotation = true;
    auto k = 0;
    for (auto& [c, c_value] : crossovers) {
        for (auto& [m, m_value] : mutations) {
            for (auto& [s, s_value] : selectionPressure) {
                std::cout << k++ << '\n';
                for (int i = 0; i < 3; ++i) {

                    auto allF = 0.0;
                    for (const auto& f : func) {
                        ga::GeneticAlgorithm ga{c,
                                                m,
                                                m * hypermutationRate,
                                                elitesPercentage,
                                                s,
                                                crossoverType,
                                                hillclimbingType,
                                                populationSize,
                                                dimensions,
                                                stepsToHypermutation,
                                                encodingChangeRate,
                                                maxNoImprovementSteps,
                                                f,
                                                applyShift,
                                                applyRotation};
                        auto rez = ga.run();
                        
                        allF += rez;
                    }
                    
                    // separate experiment for cf02, because it has big values
                    ga::GeneticAlgorithm ga{c,
                                            m,
                                            m * hypermutationRate,
                                            elitesPercentage,
                                            s,
                                            crossoverType,
                                            hillclimbingType,
                                            populationSize,
                                            dimensions,
                                            stepsToHypermutation,
                                            encodingChangeRate,
                                            maxNoImprovementSteps,
                                            "cf02",
                                            applyShift,
                                            applyRotation};
                    auto rez = ga.run() / 30.000; // maybe another value is better
                    allF += rez;
                    c_value += allF;
                    m_value += allF;
                    s_value += allF;

                    if (allF < minVal) {
                        minVal = allF;
                        minCrossover = c;
                        minMutation = m;
                        minSelection = s;
                    }
                }
            }
        }
    }

    std::ofstream fout{"experiments/runExperiments1.txt"};
    fout << "minVal : " << minVal << '\n';
    fout << "minCrossover : " << minCrossover << '\n';
    fout << "minMutation : " << minMutation << '\n';
    fout << "minSelection : " << minSelection << '\n';
    fout << "minSelection : " << minSelection << '\n';

    fout << "\ncrossovers\n";
    for (const auto [key, val] : crossovers) {
        fout << key << " -> " << val << '\n';
    }
    fout << "\nmutations\n";
    for (const auto [key, val] : mutations) {
        fout << key << " -> " << val << '\n';
    }
    fout << "\nselectionPressure\n";
    for (const auto [key, val] : selectionPressure) {
        fout << key << " -> " << val << '\n';
    }
}

void runExperiments2()
{
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
