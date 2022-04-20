#include "cec22/Cec22.h"
#include "functions/FunctionManager.h"
#include "pso/PSO.h"

#include <execution>
#include <fstream>
#include <iostream>

void runDefault();
void runTest();
void runExperiment(int dimensions, double inertia, double cognition,
                   double social, double chaosCoef, bool augment);

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[])
{
    // std::cout << cec22::sanity_check() << '\n';
    // runDefault();
    // runTest();
    runExperiment(10, 0.3, 1.0, 3.0, 0.001, true);
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

void testVector(const std::vector<double>& v, std::string_view func,
                int dimensions)
{
    auto x = v;
    auto aux = v;
    auto f = function_layer::FunctionManager{func, dimensions, true, true};
    std::cout << f(x, aux);
}

void runTest()
{
    auto func = "levy_func";
    auto v = std::vector<double>{-24.189693, -1.571958, 24.609199, 54.102313,
                                 -22.873517, -4.937637, 15.888914, -1.903588,
                                 -54.722542, 40.666694};
    testVector(v, func, 10);
}

double runOnce(std::string_view functionName, int dimensions, double inertia,
               double cognition, double social, double chaosCoef, bool augment)
{
    auto pso = pso::PSO(functionName, dimensions, 100, inertia, cognition,
                        social, chaosCoef, augment, true, true);
    return pso.run();
}

std::vector<double>
run30Times(std::string_view functionName, int dimensions, double inertia,
           double cognition, double social, double chaosCoef, bool augment)
{
    auto ret = std::vector<double>(30, -100.0);
    std::transform(std::execution::par_unseq, ret.begin(), ret.end(),
                   ret.begin(), [&]([[maybe_unused]] const auto& x) {
                       return runOnce(functionName, dimensions, inertia,
                                      cognition, social, chaosCoef, augment);
                   });
    return ret;
}

void runExperiment(int dimensions, double inertia, double cognition,
                   double social, double chaosCoef, bool augment)
{
    auto functions = std::vector<std::string>{"zakharov_func",
                                              "rosenbrock_func",
                                              "schaffer_F7_func",
                                              "rastrigin_func",
                                              "levy_func",
                                              "hf01",
                                              "hf02",
                                              "hf03",
                                              "cf01",
                                              "cf02",
                                              "cf03",
                                              "cf04"};
    for (auto& f : functions) {
        auto rez = run30Times(f, dimensions, inertia, cognition, social,
                              chaosCoef, augment);
        const auto mean =
            std::accumulate(rez.begin(), rez.end(), 0.0) / rez.size();
        std::cout << f << " " << mean << std::endl;

        auto fileName =
            "experiments/" + f + '_' + std::to_string(dimensions) + '_' +
            std::to_string(inertia) + '_' + std::to_string(cognition) + '_' +
            std::to_string(social) + '_' + std::to_string(chaosCoef) + '_' +
            std::to_string(augment);
        std::ofstream file{fileName};
        for (auto x : rez) {
            file << x << ' ';
        }
        file << '\n';
    }
}
