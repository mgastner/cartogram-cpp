#include "time_tracker.h"

void TimeTracker::start(const std::string &task_name)
{
  start_times[task_name] = std::chrono::steady_clock::now();
}

void TimeTracker::stop(const std::string &task_name)
{
  auto iter = start_times.find(task_name);
  if (iter != start_times.end()) {
    auto now = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      now - iter->second);
    durations[task_name] += duration;
    start_times.erase(iter);
  } else {
    std::cerr << "Error: Task " << task_name << " was not started."
              << std::endl;
  }
}

void TimeTracker::print_summary_report() const
{
  std::cout << "\n********** Time Report **********" << std::endl;
  for (const auto &pair : durations) {
    std::cout << pair.first << ": " << pair.second.count() << " ms"
              << std::endl;
  }
  std::cout << "*********************************" << std::endl;
}
