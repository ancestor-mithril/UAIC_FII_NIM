#pragma once
#include <chrono>
#include <map>
#include <sstream>
#include <string>

namespace utils::timer {

// TODO: prepare for multithreading timer
class Timer
{
  public:
    Timer() = delete;
    Timer(const std::string& name);
    ~Timer();

    static void clean();
    static std::string getStatistics();

  private:
    std::chrono::high_resolution_clock::time_point start;
    std::string name;

    static std::map<std::string, int> timers;
};

} // namespace utils::timer
