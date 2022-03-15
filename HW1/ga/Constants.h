#pragma once

#include "Constexpr.h"

namespace ga::constants {

inline constexpr auto minimum = -100.0;
inline constexpr auto maximum = 100.0;
inline constexpr auto valuesRange = maximum - minimum;
inline constexpr auto precision = 8;
inline constexpr auto count = pow_c(10, precision) * valuesRange;
inline constexpr auto bitsPerVariable = ceil_log2(count);
inline constexpr auto discriminator = (1ll << bitsPerVariable) - 1.0;

} // ga::constants