#include "Utils.h"

namespace utils {

double l2dSquared(const std::vector<double>& a, const std::vector<double>& b)
{
    // assuming a.size() == b.size()
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0,
                              std::plus<double>(),
                              [](auto x, auto y) { return (x - y) * (x - y); });
}

double l2d(const std::vector<double>& a, const std::vector<double>& b)
{
    return std::sqrt(l2dSquared(a, b));
}
} // namespace utils
