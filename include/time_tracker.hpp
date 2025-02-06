#ifndef TIME_TRACKER_H
#define TIME_TRACKER_H

#include <chrono>
#include <string>
#include <unordered_map>

class TimeTracker
{
private:
  std::unordered_map<std::string, std::chrono::steady_clock::time_point>
    start_times_;
  std::unordered_map<std::string, std::chrono::milliseconds> durations_;
  std::string name_;

public:
  void set_name(std::string);
  void start(const std::string &task_name);
  void stop(const std::string &task_name);
  void swap(const std::string &t1, const std::string &t2);
  void print_summary_report() const;

  // Find the duration of a particular task
  std::chrono::milliseconds duration(const std::string &task_name) const;
};

#endif  // TIME_TRACKER_H
