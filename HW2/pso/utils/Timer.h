#pragma once
#include <chrono>
#include <string>
#include <sstream>
#include <map>

namespace utils::timer {

// TODO: Create .cpp file
class Timer {
public:
    Timer() = delete;
    Timer(const std::string& name) : name{name} {
        start = std::chrono::high_resolution_clock::now();
    }
    ~Timer() {
        const auto stop = std::chrono::high_resolution_clock::now();
        const auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        timers[name] += duration.count();
    }

    static void clean() {
        timers = {};
    }

    static std::string getStatistics() {
        // TODO: operator << might be nicer
        std::stringstream ss;
        ss << "Timer statistics:\n";
        for (const auto& [name, time] : timers) {
            ss << name << ": " << time / 1e6 << " s " << time % (int) 1e6 / 1e3 << " ms\n";
        }
        return ss.str();
    }
    
    private:
     std::chrono::high_resolution_clock::time_point start;
     std::string name;

     static std::map<std::string, int> timers;

};

} // namespace utils::timer
