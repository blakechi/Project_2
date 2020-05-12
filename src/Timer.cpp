#include <iostream>
#include <chrono>
#include <string>

#include "Timer.hpp"


Timer::Timer()
{
    std::cout << "[Timer] Start...\n";
    _start = std::chrono::high_resolution_clock::now();
}

Timer::Timer(const std::string& message)
{
    std::cout << "[Timer] " << message << '\n';
    std::cout << "[Timer]     Start...\n";
    _start = std::chrono::high_resolution_clock::now();
}

Timer::~Timer()
{
    stop();
    duration();
}

void Timer::stop()
{
    std::cout << "[Timer]     Stopped\n";
    _end = std::chrono::high_resolution_clock::now();
}

void Timer::duration() const
{
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(_end - _start);
    std::cout << "[Timer]     Execution Time : " 
        << duration.count()/static_cast<float>(1000*1000) << " (s), " 
        << duration.count()/static_cast<float>(1000) << " (ms)\n";
}