#pragma once


class Timer
{
private:
    std::chrono::steady_clock::time_point _start;
    std::chrono::steady_clock::time_point _end;

public:
    Timer();
    Timer(const std::string& message);
    ~Timer();

    void stop();
    void duration() const;
};