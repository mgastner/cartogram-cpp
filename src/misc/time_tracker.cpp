#include "time_tracker.hpp"
#include <iostream>

void TimeTracker::start(const std::string &task_name)
{
  start_times_[task_name] = std::chrono::steady_clock::now();
}

void TimeTracker::stop(const std::string &task_name)
{
  auto iter = start_times_.find(task_name);
  if (iter != start_times_.end()) {
    auto now = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      now - iter->second);
    durations_[task_name] += duration;
    start_times_.erase(iter);
  } else {
    std::cerr << "Error: Task " << task_name << " was not started."
              << std::endl;
  }
}

void TimeTracker::print_summary_report() const
{
  std::cerr << "\n********** Time Report **********" << std::endl;
  for (const auto &pair : durations_) {
    std::cerr << pair.first << ": " << pair.second.count() << " ms"
              << std::endl;
  }
  std::cerr << "*********************************" << std::endl;
}
