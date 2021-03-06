#pragma once

#include <concepts>
#include <type_traits>

namespace ga::utils {
/// https://stackoverflow.com/a/27270730/18441695
template <typename T> T consteval pow_c(T base, std::integral auto exponent)
{
    return exponent == 0 ? 1 : base * pow_c(base, exponent - 1);
}

/// https://stackoverflow.com/a/23782939/18441695
std::size_t consteval floor_log2(std::size_t x)
{
    return x == 1 ? 0 : 1 + floor_log2(x >> 1);
}

// floor_log2 is at most unsigned
unsigned consteval ceil_log2(std::size_t x)
{
    return x == 1 ? 0 : floor_log2(x - 1) + 1;
}
} // namespace ga::utils
