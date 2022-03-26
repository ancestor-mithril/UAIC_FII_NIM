#include "FunctionManager.h"

#include "Cec22.h"

#include <exception>
// #include <format> // ;( format is not available yet
#include <filesystem>
#include <fstream>
#include <iostream>
#include <ranges>
#include <set>
#include <sstream>
#include <unordered_map>

namespace fs = std::filesystem;

namespace ga {

namespace {

std::string getInputDir()
{
    if (auto p = std::getenv("GA_ROOT")) {
        return std::string{p} + "/input_data/";
    }
    return fs::current_path().string() + "/input_data/";
    // TODO: implement recursive search
}

std::vector<double> readShift(std::size_t dimensions, int index, bool shiftFlag)
{
    if (not shiftFlag) {
        return {};
    }

    // std format is not available yet
    const auto file =
        getInputDir() + "shift_data_" + std::to_string(index) + ".txt";
    if (not fs::exists(file)) {
        throw std::runtime_error{"File " + file + " does not exist"};
    }

    // std::cout << "Reading file " << file << '\n';
    std::ifstream in{file};
    std::vector<double> x(std::istream_iterator<double>{in},
                          std::istream_iterator<double>{});

    if (x.size() < dimensions) {
        std::cerr << "Read " << x.size() << " doubles, expected " << dimensions
                  << "\n";
        throw std::runtime_error{"Read error"};
    }
    x.resize(dimensions);
    return x;
}

std::vector<std::vector<double>>
readRotate(std::size_t rows, std::size_t columns, int index, bool rotateFlag)
{
    if (not rotateFlag) {
        return {{}};
    }

    const auto file = getInputDir() + "M_" + std::to_string(index) + "_D" +
                      std::to_string(columns) + ".txt";
    if (not fs::exists(file)) {
        throw std::runtime_error{"File " + file + " does not exist"};
    }

    // std::cout << "Reading file " << file << '\n';
    std::ifstream in{file};
    std::vector<std::vector<double>> rotate;
    for (std::string str; std::getline(in, str);) {
        std::istringstream ss{str};
        rotate.push_back({std::istream_iterator<double>{ss},
                          std::istream_iterator<double>{}});
    }

    if (rotate.size() < rows) {
        std::cerr << "Rotate has " << rotate.size() << " rows, expected "
                  << rows << "\n";
        throw std::runtime_error{"Read error"};
    }
    rotate.resize(rows);
    for (auto& row : rotate) {
        if (row.size() < columns) {
            std::cerr << "Row has " << row.size() << " columns, expected "
                      << columns << "\n";
            throw std::runtime_error{"Read error"};
        }
        row.resize(columns);
    }
    return rotate;
}

std::vector<std::size_t> readShuffle(std::size_t dimensions, int index)
{
    const auto file = getInputDir() + "shuffle_data_" + std::to_string(index) +
                      "_D" + std::to_string(dimensions) + ".txt";
    if (not fs::exists(file)) {
        throw std::runtime_error{"File " + file + " does not exist"};
    }

    // std::cout << "Reading file " << file << '\n';
    std::ifstream in{file};
    std::vector<std::size_t> x(std::istream_iterator<std::size_t>{in},
                               std::istream_iterator<std::size_t>{});

    // TODO: Debug
    // Validation
    std::set<std::size_t> set;
    for (auto& i : x) {
        --i;
        if (i > dimensions) {
            // i can't be smaller than 0 because it's size_t
            std::cerr << "I is " << i << ", max accepted value is "
                      << dimensions - 1 << '\n';
            throw std::runtime_error{"Read error"};
        }
        auto [it, inserted] = set.insert(i);
        if (inserted) {
            std::cerr << "Index is repeating: " << i << '\n';
            throw std::runtime_error{"Read error"};
        }
    }
    if (set.size() != dimensions) {
        std::cerr << "Not valid indices from 0 to " << dimensions - 1;
        throw std::runtime_error{"Read error"};
    }

    return x;
}

} // namespace

FunctionManager::FunctionManager(const std::string& functionName,
                                 int dimensions, bool shiftFlag,
                                 bool rotateFlag)
    : functionName{functionName}
{

    if ((shiftFlag or rotateFlag) and dimensions != 10 and dimensions != 20) {
        // not very flexible...
        throw std::runtime_error{"This Function Manager accepts only 10 or 20 "
                                 "dimensions when shifting or rotating"};
    }
    function = initFunction(dimensions, shiftFlag, rotateFlag);
}

std::string FunctionManager::toString() const
{
    std::ostringstream ss;
    ss << functionCalls << '\n';
    for (const auto v : values) {
        ss << v << ' ';
    }
    ss << '\n';
    return ss.str();
}

double
FunctionManager::operator()(std::vector<double>& x, std::vector<double>& aux)
{
    ++functionCalls;
    auto ret = function(x, aux);
    if (functionCalls % 1000) {
        values.push_back(ret);
    }
    return ret;
}

std::function<double(std::vector<double>&, std::vector<double>&)>
FunctionManager::initFunction(int dimensions, bool shiftFlag, bool rotateFlag)
{
    using namespace std::string_literals;
    using namespace ga::functions;
    const std::unordered_map<
        std::string,
        std::tuple<int,
                   std::function<double(
                       std::vector<double>&, std::vector<double>&,
                       const std::vector<double>&,
                       const std::vector<std::vector<double>>&, bool, bool)>,
                   double>>
        basicFunctions = {
            // TODO: find what 300.0 and the likes are and if we even need to
            // add them
            {"zakharov_func"s, {1, zakharov_func, 300.0}},
            {"rosenbrock_func"s, {2, rosenbrock_func, 400.0}},
            {"schaffer_F7_func"s, {3, schaffer_F7_func, 600.0}},
            {"rastrigin_func"s, {4, rastrigin_func, 800.0}},
            {"levy_func"s, {5, levy_func, 900.0}},
        };
    const std::unordered_map<
        std::string,
        std::tuple<
            int,
            std::function<double(std::vector<double>&, std::vector<double>&,
                                 const std::vector<double>&,
                                 const std::vector<std::vector<double>>&,
                                 const std::vector<std::size_t>&, bool, bool)>,
            double>>
        hybridFunctions = {
            {"hf01"s, {6, hf01, 1800.0}},
            {"hf02"s, {7, hf02, 2000.0}},
            {"hf03"s, {8, hf03, 2200.0}},
        };
    const std::unordered_map<
        std::string,
        std::tuple<int,
                   std::function<double(
                       std::vector<double>&, std::vector<double>&,
                       const std::vector<double>&,
                       const std::vector<std::vector<double>>&, bool)>,
                   double, int>>
        compositionFunctions = {
            {"cf01"s, {9, cf01, 2300.0, 5}},
            {"cf02"s, {10, cf02, 2400.0, 3}},
            {"cf03"s, {11, cf03, 2600.0, 5}},
            {"cf04"s, {12, cf04, 2700.0, 6}},
        };

    if (basicFunctions.find(functionName) != basicFunctions.end()) {
        const auto [index, f, fStar] = basicFunctions.at(functionName);
        return [=, f = std::move(f),
                shift = readShift(dimensions, index, shiftFlag),
                rotate = readRotate(dimensions, dimensions, index, rotateFlag)](
                   std::vector<double>& x, std::vector<double>& aux) {
            return f(x, aux, shift, rotate, shiftFlag, rotateFlag);
            // + fStar;
        };
    }

    if (hybridFunctions.find(functionName) != hybridFunctions.end()) {
        const auto [index, f, fStar] = hybridFunctions.at(functionName);
        return [=, f = std::move(f),
                shift = readShift(dimensions, index, shiftFlag),
                rotate = readRotate(dimensions, dimensions, index, rotateFlag),
                indices = readShuffle(dimensions, index)](
                   std::vector<double>& x, std::vector<double>& aux) {
            return f(x, aux, shift, rotate, indices, shiftFlag, rotateFlag);
            // + fStar;
        };
    }

    if (compositionFunctions.find(functionName) != compositionFunctions.end()) {
        const auto [index, f, fStar, n] = compositionFunctions.at(functionName);
        return [=, f = std::move(f),
                shift = readShift(dimensions * n, index, true),
                rotate = readRotate(dimensions * n, dimensions, index, true)](
                   std::vector<double>& x, std::vector<double>& aux) {
            return f(x, aux, shift, rotate, true); // always rotate
            // + fStar;
        };
    }

    // No function found
    std::cerr << "Available basic functions: \n";
    for (const auto& f : std::views::keys(basicFunctions)) {
        std::cerr << '\t' << f << '\n';
    }
    std::cerr << "\nAvailable hybrid functions: \n";
    for (const auto& f : std::views::keys(hybridFunctions)) {
        std::cerr << '\t' << f << '\n';
    }
    std::cerr << "\nAvailable composition functions: \n";
    for (const auto& f : std::views::keys(compositionFunctions)) {
        std::cerr << '\t' << f << '\n';
    }
    throw std::runtime_error{"Function "s + functionName + " not found."s};
}

} // namespace ga
