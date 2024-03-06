#ifndef PROGRESS_TRACKER_H
#define PROGRESS_TRACKER_H

#include "inset_state.hpp"

class ProgressTracker
{
public:
  // Constructor
  explicit ProgressTracker(double total_geo_divs);

  // Destructor
  ~ProgressTracker();

  // Method to update the progress
  void print_progress_mid_integration(const InsetState &inset_state);
  void update_and_print_progress(const InsetState &inset_state);

  // Method to print the current progress
  void print_progress(const double progress);

private:
  double total_geo_divs_;  // Total number of GeoDivs to monitor progress
  double progress_;  // Progress measured on a scale from 0 (start) to 1 (end)
};

#endif  // PROGRESS_TRACKER_H
