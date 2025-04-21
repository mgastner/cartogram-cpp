#include "time_tracker.hpp"
#include <iostream>
#include <vector>

void TimeTracker::set_name(std::string name)
{
  name_ = name;
}

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

void TimeTracker::swap(const std::string &t1, const std::string &t2)
{
  stop(t1);
  start(t2);
}

void TimeTracker::print_summary_report() const
{
  std::cerr << "\n********** Time Report **********" << std::endl;
  std::cerr << "(" << name_ << ")\n" << std::endl;
  std::vector<std::pair<std::string, std::chrono::milliseconds>>
    sorted_durations;
  for (const auto &[task, time_taken] : durations_) {
    sorted_durations.push_back({task, time_taken});
  }
  std::sort(
    sorted_durations.begin(),
    sorted_durations.end(),
    [](const auto &a, const auto &b) {
      return a.second > b.second;
    });

  for (const auto &[task, time_taken] : sorted_durations) {
    std::cerr << task << ": " << time_taken.count() << " ms" << std::endl;
  }
  std::cerr << "*********************************" << std::endl;
}

std::chrono::milliseconds TimeTracker::duration(
  const std::string &task_name) const
{
  try {
    return durations_.at(task_name);
  } catch (const std::out_of_range &e) {
    std::cerr << "ERROR: Key '" << task_name
              << "' not found in target_areas_. "
              << "Exception: " << e.what() << std::endl;
    // Re-throw, or return a default value
    throw;
  }
}