// timing_harness.hpp

#ifndef TIMING_HARNESS_HPP
#define TIMING_HARNESS_HPP

#include <chrono>
#include <string>
#include <unordered_map>

// Starts timing a function
void time_fun_start(const std::string& function_name);

// Ends timing a function and records the elapsed time
void time_fun_end(const std::string& function_name);

// Writes the timing data to a CSV file
void write_timing_data(const std::string& filename);

// Macro for convenient timing of function calls
#define TIMED_CALL(func, ...)               \
    do {                                    \
        time_fun_start(#func);              \
        func(__VA_ARGS__);                  \
        time_fun_end(#func);                \
    } while(0)

#endif // TIMING_HARNESS_HPP

