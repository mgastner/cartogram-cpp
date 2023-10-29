#include "progress_tracker.h"
#include "constants.h"

// Constructor
ProgressTracker::ProgressTracker(double total_geo_divs)
    : total_geo_divs(total_geo_divs), progress(0.0)
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
  const double inset_max_frac = inset_state.n_geo_divs() / total_geo_divs;
  print_progress(progress + (inset_max_frac / n_predicted_integrations));
}

// Method to update the progress and print it
void ProgressTracker::update_and_print_progress(const InsetState &inset_state)
{
  const double inset_max_frac = inset_state.n_geo_divs() / total_geo_divs;
  progress += inset_max_frac;
  print_progress(progress);
}

// Method to print the current progress
void ProgressTracker::print_progress(double progress)
{
  std::cerr << "Progress: " << progress << std::endl;
}