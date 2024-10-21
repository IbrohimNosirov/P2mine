// timing_harness.cpp

#include "timing_harness.hpp"
#include <fstream>
#include <iostream>
#include <mutex>
#include <unordered_map>

static std::unordered_map<std::string, double> function_times;
static std::unordered_map<std::string, std::chrono::high_resolution_clock::time_point> start_times;
static std::mutex timing_mutex;

void time_fun_start(const std::string& function_name) {
    std::lock_guard<std::mutex> lock(timing_mutex);
    start_times[function_name] = std::chrono::high_resolution_clock::now();
}

void time_fun_end(const std::string& function_name) {
    auto end_time = std::chrono::high_resolution_clock::now();
    std::lock_guard<std::mutex> lock(timing_mutex);
    auto start_time_it = start_times.find(function_name);
    if (start_time_it != start_times.end()) {
        double elapsed = std::chrono::duration<double>(end_time - start_time_it->second).count();
        function_times[function_name] += elapsed;
        start_times.erase(start_time_it); // Remove the start time entry
    } else {
        std::cerr << "Warning: time_fun_end called without matching time_fun_start for function "
                  << function_name << std::endl;
    }
}

void write_timing_data(const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << "Function,Time(s)\n";
        for (const auto& entry : function_times) {
            file << entry.first << "," << entry.second << "\n";
        }
        file.close();
    } else {
        std::cerr << "Error opening timing data file for writing: " << filename << std::endl;
    }
}

