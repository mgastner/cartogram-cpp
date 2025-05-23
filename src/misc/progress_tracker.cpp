#include "progress_tracker.hpp"
#include "constants.hpp"

// Constructor
ProgressTracker::ProgressTracker(double total_geo_divs)
    : total_geo_divs_(total_geo_divs), progress_(0),
      bar_(
        indicators::ProgressBar(
          indicators::option::BarWidth{75},
          indicators::option::Start{"["},
          indicators::option::Fill{"■"},
          indicators::option::Lead{"■"},
          indicators::option::Remainder{"-"},
          indicators::option::End{"]"},
          indicators::option::FontStyles(
            std::vector<indicators::FontStyle>{indicators::FontStyle::bold}),
          indicators::option::ShowPercentage{true},
          indicators::option::ShowElapsedTime{true},
          indicators::option::Stream{std::cerr}))
{
  bar_.set_progress(0);  // Initialize progress to 0 at the start
}

// Destructor
ProgressTracker::~ProgressTracker() = default;

// Method to print the current progress mid integration
void ProgressTracker::print_progress_mid_integration(
  double max_area_error,
  unsigned int n_geo_div_in_inset,
  unsigned int n_finished_integrations)
{
  // Calculate progress percentage. We assume that the maximum area
  // error is typically reduced to 1/5 of the previous value.
  const double ratio_actual_to_permitted_max_area_error =
    max_area_error / max_permitted_area_error;
  const double n_predicted_integrations =
    std::max((log(ratio_actual_to_permitted_max_area_error) / log(5)), 1.0);

  // We make the approximation that the progress towards generating the
  // cartogram is proportional to the number of GeoDivs that are in the
  // finished insets
  const double inset_max_frac = n_geo_div_in_inset / total_geo_divs_;
  double progress = progress_ + (inset_max_frac / n_predicted_integrations);

  // Change how much progress increases by, so it never reaches 100 here
  double remaining_progress = 1.0 - max_progress_;
  double dynamic_increment = remaining_progress * 0.1;

  // Leave buffer at end so that we don't reach 100% prematurely
  progress = std::min(progress, 0.75);

  // Our assumption above causes the progress bar to start at 36%.
  // Thus, we temper it down for the first few integrations.
  if (n_finished_integrations < 4) {
    progress = std::min(progress, max_progress_);
  }

  // Increase max_progress by dynamic increment that gets smaller
  // as we get closer to 100%.
  progress = std::max(progress, max_progress_ + dynamic_increment);

  max_progress_ = progress;
  print_progress(progress);
  print_progress_bar(progress);
}

// Method to update the progress and print progress at the end of the
// integrations of the inset
void ProgressTracker::update_and_print_progress_end_integration(
  const unsigned int n_geo_divs_in_inset)
{
  max_progress_ = 0;
  const double inset_max_frac = n_geo_divs_in_inset / total_geo_divs_;
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
  bar_.set_progress(static_cast<size_t>(round(progress * 100)));
  std::cerr << std::endl;  // Print a new line after the progress bar
}