#include "Timer.h"

namespace utils::timer {
std::map<std::string, int> Timer::timers = {};

Timer::Timer(const std::string& name) : name{name}
{
    start = std::chrono::high_resolution_clock::now();
}

Timer::~Timer()
{
    const auto stop = std::chrono::high_resolution_clock::now();
    const auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    timers[name] += duration.count();
}

// static
void Timer::clean()
{
    timers = {};
}

// static
std::string Timer::getStatistics()
{
    // TODO: operator << might be nicer
    std::stringstream ss;
    ss << "Timer statistics:\n";
    for (const auto& [name, time] : timers) {
        ss << name << ": " << time / 1e6 << " s\n";
    }
    return ss.str();
}

} // namespace utils::timer
