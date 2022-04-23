#pragma once

#include <cmath>
#include <functional>
#include <numeric>
#include <vector>

namespace utils {

double l2dSquared(const std::vector<double>& a, const std::vector<double>& b);

double l2d(const std::vector<double>& a, const std::vector<double>& b);

struct L2Distance {
    double
    operator()(const std::vector<double>& a, const std::vector<double>& b) const
    {
        return l2d(a, b);
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
