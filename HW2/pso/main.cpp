#include "cec22/Cec22.h"
#include "functions/FunctionManager.h"
#include "pso/PSO.h"
#include "utils/Timer.h"

#include <chrono>
#include <execution>
#include <fstream>
#include <iostream>
#include <thread>

using cacheStrategy =
    function_layer::cache_layer::KDTreeCache::CacheRetrievalStrategy;

void runDefault();
void runTest();
void runExperiment(int dimensions, int resetThreshold, double inertia,
                   double cognition, double social, double chaosCoef,
                   cacheStrategy cacheRetrievalStrategy, pso::topology topology,
                   bool augment);
void timeTest();

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[])
{
    // std::cout << cec22::sanity_check() << '\n';
    // runDefault();
    // runTest();

    runExperiment(10, 100, 0.3, 1.0, 3.0, 0.001, cacheStrategy::FirstNeighbor,
                  pso::topology::StaticRing, true);

    // timeTest();
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
    auto f = function_layer::FunctionManager{
        func, dimensions, cacheStrategy::FirstNeighbor, true, true};
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
runOnce(std::string_view functionName, int dimensions, int resetThreshold,
        double inertia, double cognition, double social, double chaosCoef,
        cacheStrategy cacheRetrievalStrategy, pso::topology topology,
        bool augment)
{
    auto pso = pso::PSO(functionName, dimensions, 100, resetThreshold, inertia,
                        cognition, social, chaosCoef, cacheRetrievalStrategy,
                        topology, augment, true, true);
    // TODO: Add time measurements and write them to file
    auto value = pso.run();

    auto cacheHits = pso.getCacheHits();
    return {value, cacheHits};
}

std::vector<std::pair<double, int>>
run30Times(std::string_view functionName, int dimensions, int resetThreshold,
           double inertia, double cognition, double social, double chaosCoef,
           cacheStrategy cacheRetrievalStrategy, pso::topology topology,
           bool augment, int runs)
{
    auto ret = std::vector<std::pair<double, int>>(runs, {-100.0, 0});
    std::transform(std::execution::par_unseq, ret.begin(), ret.end(),
                   ret.begin(), [=]([[maybe_unused]] const auto& x) {
                       return runOnce(functionName, dimensions, resetThreshold,
                                      inertia, cognition, social, chaosCoef,
                                      cacheRetrievalStrategy, topology,
                                      augment);
                   });
    return ret;
}

void timeTest()
{
    const auto functions = std::vector<std::string>{
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
    for (auto i = 1; i < 10; ++i) {
        std::cout << i << std::endl;
        function_layer::FunctionManager::rebalance = i;

        for (auto& func : functions) {
            run30Times(func, 10, i, 0.3, 1.0, 3.0, 0.001,
                       cacheStrategy::FirstNeighbor, pso::topology::Star, true,
                       10);
        }
        std::cout << utils::timer::Timer::getStatistics() << std::endl;
        utils::timer::Timer::clean();
    }
    for (auto i = 11; i < 50; i += 5) {
        std::cout << i << std::endl;
        function_layer::FunctionManager::rebalance = i;

        for (auto& func : functions) {
            run30Times(func, 10, i, 0.3, 1.0, 3.0, 0.001,
                       cacheStrategy::FirstNeighbor, pso::topology::Star, true,
                       10);
        }
        std::cout << utils::timer::Timer::getStatistics() << std::endl;
        utils::timer::Timer::clean();
    }
    for (auto i = 51; i < 100; i += 10) {
        std::cout << i << std::endl;
        function_layer::FunctionManager::rebalance = i;

        for (auto& func : functions) {
            run30Times(func, 10, i, 0.3, 1.0, 3.0, 0.001,
                       cacheStrategy::FirstNeighbor, pso::topology::Star, true,
                       10);
        }
        std::cout << utils::timer::Timer::getStatistics() << std::endl;
        utils::timer::Timer::clean();
    }
}

double runForFunction(std::string_view f, int dimensions, int resetThreshold,
                      double inertia, double cognition, double social,
                      double chaosCoef, cacheStrategy cacheRetrievalStrategy,
                      pso::topology topology, bool augment)
{
    const auto start = std::chrono::high_resolution_clock::now();
    const auto rez =
        run30Times(f, dimensions, resetThreshold, inertia, cognition, social,
                   chaosCoef, cacheRetrievalStrategy, topology, augment, 30);
    const auto end = std::chrono::high_resolution_clock::now();
    const auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    const auto miliseconds = duration.count();

    const auto mean = std::accumulate(rez.begin(), rez.end(), 0.0,
                                      [](auto f, auto elem) {
                                          return std::move(f) + elem.first;
                                      }) /
                      rez.size();

    std::cout << f << " " << mean << " and took " << miliseconds / 1000.0
              << " seconds " << std::endl;
    std::cout << utils::timer::Timer::getStatistics() << std::endl;
    utils::timer::Timer::clean();

    const auto fileName =
        "experiments/" + std::string{f} + '_' + std::to_string(dimensions) +
        '_' + std::to_string(inertia) + '_' + std::to_string(resetThreshold) +
        '_' + std::to_string(cognition) + '_' + std::to_string(social) + '_' +
        std::to_string(chaosCoef) + '_' +
        std::to_string((int)cacheRetrievalStrategy) + '_' +
        std::to_string(augment) + "2";
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

void runExperiment(int dimensions, int resetThreshold, double inertia,
                   double cognition, double social, double chaosCoef,
                   cacheStrategy cacheRetrievalStrategy, pso::topology topology,
                   bool augment)
{
    const auto functions = std::vector<std::string>{
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
        // futures.push_back(std::jthread{
        //     runForFunction, f, dimensions, resetThreshold, inertia,
        //     cognition, social, chaosCoef, cacheRetrievalStrategy, augment});
        runForFunction(f, dimensions, resetThreshold, inertia, cognition,
                       social, chaosCoef, cacheRetrievalStrategy, topology,
                       augment);
    }
    for (auto& f : futures) {
        f.join();
    }
    // std::cout << meanSum << '\n';
}
