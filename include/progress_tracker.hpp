#ifndef PROGRESS_TRACKER_H
#define PROGRESS_TRACKER_H

#include "indicators/block_progress_bar.hpp"
#include "indicators/color.hpp"
#include "indicators/cursor_control.hpp"
#include "indicators/cursor_movement.hpp"
#include "indicators/details/stream_helper.hpp"
#include "indicators/display_width.hpp"
#include "indicators/dynamic_progress.hpp"
#include "indicators/font_style.hpp"
#include "indicators/indeterminate_progress_bar.hpp"
#include "indicators/multi_progress.hpp"
#include "indicators/progress_bar.hpp"
#include "indicators/progress_spinner.hpp"
#include "indicators/progress_type.hpp"
#include "indicators/setting.hpp"
#include "indicators/termcolor.hpp"
#include "indicators/terminal_size.hpp"
#include <cstdint>

class ProgressTracker
{
public:
  // Constructor
  explicit ProgressTracker(double);

  // Destructor
  ~ProgressTracker();

  // Method to update the progress
  void print_progress_mid_integration(
    double max_area_error,
    unsigned int n_geo_div_in_inset,
    unsigned int n_finished_integrations);
  void update_and_print_progress_end_integration(
    const unsigned int n_geo_divs_in_inset);

  // Method to print the current progress
  void print_progress(const double);

  // Method to print the progress bar
  void print_progress_bar(const double);

private:
  double total_geo_divs_;  // Total number of GeoDivs to monitor progress
  double progress_;  // Progress measured on a scale from 0 (start) to 1 (end)
  double max_progress_;  // Maximum progress value ever reached
  indicators::ProgressBar bar_;
};

#endif  // PROGRESS_TRACKER_H
