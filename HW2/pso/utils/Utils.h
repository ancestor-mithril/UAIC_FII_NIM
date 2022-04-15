#pragma once

#include <cmath>
#include <functional>
#include <numeric>
#include <vector>

namespace utils {

double l2distance(const std::vector<double>& a, const std::vector<double>& b)
{
    // assuming a.size() == b.size()
    return std::sqrt(std::inner_product(
        a.begin(), a.end(), b.begin(), 0.0, std::plus<double>(),
        [](auto x, auto y) { return (x - y) * (x - y); }));
}

struct L2Distance {
    double
    operator()(const std::vector<double>& a, const std::vector<double>& b) const
    {
        return l2distance(a, b);
    }
};

struct L2Norm {
    double operator()(const std::vector<double>& a) const
    {
        // Accumulate vs reduce
        return std::sqrt(std::accumulate(
            a.begin(), a.end(), 0.0, [](auto x, auto y) { return x + y * y; }));
    }
};

} // namespace utils
