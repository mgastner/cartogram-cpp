#include "progress_tracker.hpp"
#include "constants.hpp"

// Constructor
ProgressTracker::ProgressTracker(double total_geo_divs)
    : total_geo_divs_(total_geo_divs), progress_(0.0)
{
}

// Destructor
ProgressTracker::~ProgressTracker() = default;

// Method to print the current progress mid integration
void ProgressTracker::print_progress_mid_integration(
  const InsetState &inset_state)
{
  // Calculate progress percentage. We assume that the maximum area
  // error is typically reduced to 1/5 of the previous value.
  const double ratio_actual_to_permitted_max_area_error =
    inset_state.max_area_error().value / max_permitted_area_error;
  const double n_predicted_integrations =
    std::max((log(ratio_actual_to_permitted_max_area_error) / log(5)), 1.0);

  // We make the approximation that the progress towards generating the
  // cartogram is proportional to the number of GeoDivs that are in the
  // finished insets
  const double inset_max_frac = inset_state.n_geo_divs() / total_geo_divs_;
  double progress = progress_ + (inset_max_frac / n_predicted_integrations);
  print_progress(progress);
  print_progress_bar(progress);
}

// Method to update the progress and print progress at the end of the
// integrations of the inset
void ProgressTracker::update_and_print_progress_end_integration(
  const InsetState &inset_state)
{
  const double inset_max_frac = inset_state.n_geo_divs() / total_geo_divs_;
  progress_ += inset_max_frac;
  print_progress(progress_);
  print_progress_bar(progress_);
}

// Method to print the current progress
void ProgressTracker::print_progress(double progress)
{
  std::cerr << "Progress: " << progress << std::endl;
}

// Method to print the progress bar
void ProgressTracker::print_progress_bar(double progress)
{
  bar_.set_progress(round(progress * 100));
  std::cerr << std::endl;
}