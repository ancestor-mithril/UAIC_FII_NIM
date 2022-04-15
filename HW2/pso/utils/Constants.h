#pragma once

namespace utils::constants {

// changeable
inline constexpr auto precision = 8;
inline constexpr auto populationSize = 100;

// problem specific
inline constexpr auto minimum = -100.0;
inline constexpr auto maximum = 100.0;
inline constexpr auto best = 1.0e-8;

// constants
inline constexpr auto valuesRange = maximum - minimum;

} // namespace utils::constants
