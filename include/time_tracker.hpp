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

  std::chrono::steady_clock::time_point program_start_;

public:
  TimeTracker();
  explicit TimeTracker(std::string name);

  void set_name(std::string);
  void start(const std::string &task_name);
  void stop(const std::string &task_name);
  void swap(const std::string &t1, const std::string &t2);
  void print_summary_report() const;

  // Find the duration of a particular task
  std::chrono::milliseconds duration(const std::string &task_name) const;

  // Total elapsed time from the CTOR of the object till now in seconds
  double total_elapsed_time_in_seconds() const;
};

#endif  // TIME_TRACKER_H
