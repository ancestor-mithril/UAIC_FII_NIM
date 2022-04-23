#include "cec22/Cec22.h"
#include "functions/FunctionManager.h"
#include "pso/PSO.h"

#include <execution>
#include <fstream>
#include <iostream>
#include <thread>

void runDefault();
void runTest();
void runExperiment(int dimensions, double inertia, double cognition,
                   double social, double chaosCoef, bool augment);

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[])
{
    // std::cout << cec22::sanity_check() << '\n';
    // runDefault();
    // runTest();
    runExperiment(10, 0.3, 1.0, 3.0, 0.01, true);
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

std::pair<double, int>
runOnce(std::string_view functionName, int dimensions, double inertia,
        double cognition, double social, double chaosCoef, bool augment)
{
    auto pso = pso::PSO(functionName, dimensions, 500, inertia, cognition,
                        social, chaosCoef, augment, true, true);
    // TODO: Add time measurements and write them to file
    auto value = pso.run();

    auto cacheHits = pso.getCacheHits();
    return {value, cacheHits};
}

std::vector<std::pair<double, int>>
run30Times(std::string_view functionName, int dimensions, double inertia,
           double cognition, double social, double chaosCoef, bool augment)
{
    auto ret = std::vector<std::pair<double, int>>(30, {-100.0, 0});
    std::transform(std::execution::par_unseq, ret.begin(), ret.end(),
                   ret.begin(), [&]([[maybe_unused]] const auto& x) {
                       return runOnce(functionName, dimensions, inertia,
                                      cognition, social, chaosCoef, augment);
                   });
    return ret;
}

double
runForFunction(std::string_view f, int dimensions, double inertia,
               double cognition, double social, double chaosCoef, bool augment)
{
    auto rez = run30Times(f, dimensions, inertia, cognition, social, chaosCoef,
                          augment);
    const auto mean = std::accumulate(rez.begin(), rez.end(), 0.0,
                                      [](auto f, auto elem) {
                                          return std::move(f) + elem.first;
                                      }) /
                      rez.size();
    std::cout << f << " " << mean << std::endl;

    auto fileName = "experiments/" + std::string{f} + '_' +
                    std::to_string(dimensions) + '_' + std::to_string(inertia) +
                    '_' + std::to_string(cognition) + '_' +
                    std::to_string(social) + '_' + std::to_string(chaosCoef) +
                    '_' + std::to_string(augment) + "2";
    std::ofstream file{fileName};
    for (auto [x, _] : rez) {
        file << x << ' ';
    }
    file << '\n';
    for (auto [_, x] : rez) {
        file << x << ' ';
    }
    file << '\n';

    return mean;
}

void runExperiment(int dimensions, double inertia, double cognition,
                   double social, double chaosCoef, bool augment)
{
    auto functions = std::vector<std::string>{
        "zakharov_func",
        "rosenbrock_func",
        "schaffer_F7_func",
        "rastrigin_func", // blocked
        "levy_func",      // blocked
        "hf01",           // semi-blocked
        "hf02",           // blocked
        "hf03",           // blocked
        "cf01",           // blocked
        "cf02",           // blocked
        "cf03",           // blocked
        "cf04"            // blocked
    };

    std::vector<std::jthread> futures;
    futures.reserve(12);
    for (auto& f : functions) {
        std::cout << f << std::endl;
        // futures.push_back(std::jthread{runForFunction, f, dimensions,
        // inertia,
        //                                cognition, social, chaosCoef,
        //                                augment});
        runForFunction(f, dimensions, inertia, cognition, social, chaosCoef,
                       augment);
    }
    for (auto& f : futures) {
        f.join();
    }
    // std::cout << meanSum << '\n';
}
