#ifndef TIME_TRACKER_H
#define TIME_TRACKER_H

#include <chrono>
#include <iostream>
#include <string>
#include <unordered_map>

class TimeTracker
{
private:
  std::unordered_map<std::string, std::chrono::steady_clock::time_point>
    start_times_;
  std::unordered_map<std::string, std::chrono::milliseconds> durations_;

public:
  void start(const std::string &task_name);
  void stop(const std::string &task_name);
  void print_summary_report() const;
};

#endif  // TIME_TRACKER_H
