#pragma once

#include "Constexpr.h"

namespace ga::constants {

// changeable
inline constexpr auto precision = 8;
inline constexpr auto populationSize = 100;

// problem specific
inline constexpr auto minimum = -100.0;
inline constexpr auto maximum = 100.0;

// constants
inline constexpr auto valuesRange = maximum - minimum;
inline constexpr auto count =
    static_cast<std::size_t>(utils::pow_c(10, precision) * valuesRange);
inline constexpr auto bitsPerVariable =
    static_cast<int>(utils::ceil_log2(count));
inline constexpr auto discriminator = (1LL << bitsPerVariable) - 1.0;

} // namespace ga::constants
