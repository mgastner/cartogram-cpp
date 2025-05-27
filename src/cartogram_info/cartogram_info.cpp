#include "cartogram_info.hpp"
#include "constants.hpp"
#include "csv.hpp"
#include "round_point.hpp"
#include <cstdlib>
#include <iostream>

CartogramInfo::CartogramInfo(const Arguments args) : args_(args)
{
  is_world_map_ = args_.world;
  timer.start("Total time");
  crs_ = "+proj=longlat";

  if (!args.visual_file_name.empty()) {

    // Determine name of input map based on the CSV file and store it
    set_map_name(args_.visual_file_name);
  } else {

    // User wants to likely make CSV or output equal area map
    set_map_name(args_.geo_file_name);
  }
  timer.set_name(map_name_);
}

double CartogramInfo::cart_initial_total_target_area() const
{
  double target_area = 0.0;

  // Iterate over inset states. Syntax from:
  // https://stackoverflow.com/questions/13087028/can-i-easily-iterate-over-
  // the-values-of-a-map-using-a-range-based-for-loop
  for (const InsetState &inset_state : inset_states_) {
    target_area += inset_state.initial_target_area();
  }
  return target_area;
}

bool CartogramInfo::converged() const
{
  for (const InsetState &inset_state : inset_states_) {
    if (!inset_state.converged()) {
      return false;
    }
  }
  return true;
}

double CartogramInfo::area() const
{
  double area = 0.0;

  // Iterate over inset states. Syntax from:
  // https://stackoverflow.com/questions/13087028/can-i-easily-iterate-over-
  // the-values-of-a-map-using-a-range-based-for-loop
  for (const InsetState &inset_state : inset_states_) {
    area += inset_state.total_inset_area();
  }
  return area;
}

bool CartogramInfo::is_world_map() const
{
  return is_world_map_;
}

size_t CartogramInfo::n_geo_divs() const
{
  size_t n_geo_divs = 0;
  for (const InsetState &inset_state : inset_states_) {
    n_geo_divs += inset_state.n_geo_divs();
  }
  return n_geo_divs;
}

size_t CartogramInfo::n_insets() const
{
  return inset_states_.size();
}

void CartogramInfo::plot_input()
{

  // Color if colors are not provided
  for (InsetState &inset_state : inset_states_) {
    inset_state.auto_color();
  }

  // Create copy of cart_info
  CartogramInfo tmp_ci = *this;

  for (InsetState &inset_state : tmp_ci.ref_to_inset_states()) {
    inset_state.normalize_inset_area(
      tmp_ci.cart_initial_total_target_area(),
      true);
  }

  // Shift insets so that they do not overlap
  tmp_ci.reposition_insets();
  tmp_ci.write_svg("input");
}

void CartogramInfo::preprocess()
{

  // Replace missing and zero target areas with positive values
  replace_missing_and_zero_target_areas();

  for (InsetState &inset_state : inset_states_) {

    // Determine the name of the inset
    std::string inset_name = map_name_;
    if (n_insets() > 1) {
      inset_name += "_" + inset_state.pos();
    }
    inset_state.set_inset_name(inset_name);

    // Preprocess inset
    inset_state.preprocess();
  }

  if (args_.export_preprocessed) {
    // Output rescaled GeoJSON
    write_geojson("input_processed");
    // processed = simplified + rescaled
    // and potentially projected + small polygons removed

    // Output preprocessed CSV file
    write_csv(map_name_ + "_input_processed");
  }
}

void CartogramInfo::project_to_equal_area()
{

  // Project map and ensure that all holes are inside polygons
  for (InsetState &inset_state : inset_states_) {

    // Check for errors in the input topology
    inset_state.check_topology();

    // Can the coordinates be interpreted as longitude and latitude?
    // TODO: The "crs" field for GeoJSON files seems to be deprecated.
    //       However, in earlier specifications, the coordinate reference
    //       system used to be written in the format specified here:
    //       https://geojson.org/geojson-spec.html#coordinate-reference-system-objects.
    //       It may be a good idea to make a list of possible entries
    //       corresponding to longitude and lattitude projection.
    //       "urn:ogc:def:crs:OGC:1.3:CRS84" is one such entry.
    if (!args_.skip_projection || args_.output_equal_area_map) {
      // TODO: Potentially check for the CRS field we output?

      // If yes, transform the coordinates with the Albers projection if the
      // input map is not a world map. Otherwise, use the Smyth-Craster
      // projection.
      if (args_.world) {
        inset_state.apply_smyth_craster_projection();
      } else {
        inset_state.apply_albers_projection();
      }
    }
  }

  if (args_.output_equal_area_map) {

    if (args_.skip_projection) {
      std::cerr << "WARNING: --skip_projection flag ignored as it "
                << "contradicts --output_equal_area_map flag. " << std::endl;
    }

    // If a visual file has also been provided, it may have insets defined.
    if (!args_.visual_file_name.empty()) {

      // If the user does not want the projection to be changed according to
      // defined insets, then, they should not provide an input visual file.
      // With such a file provided, we project once per inset, with each
      // projection having different input parameters.
      // That is, the reference longitude and latitude for the Albers
      // projection along with the standard parallels, will be calculated
      // according to each inset, rather than the entire map.
      reposition_insets();
    }
    write_geojson("equal_area");
    std::exit(EXIT_SUCCESS);
  }
}

std::vector<InsetState> &CartogramInfo::ref_to_inset_states()
{
  return inset_states_;
}

void CartogramInfo::replace_missing_and_zero_target_areas()
{
  // Get total current area and total target area
  double total_start_area_with_data = 0.0;
  double total_target_area_with_data = 0.0;
  for (const InsetState &inset_state : inset_states_) {
    for (const auto &gd : inset_state.geo_divs()) {
      if (!inset_state.target_area_is_missing(gd.id())) {
        total_start_area_with_data += gd.area();
        total_target_area_with_data += inset_state.target_area_at(gd.id());
      }
    }
  }

  const double mean_density =
    total_target_area_with_data / total_start_area_with_data;

  std::cerr << "Mean density: " << mean_density << std::endl;

  // Calculate threshold for small target areas. For GeoDivs below this
  // threshold, the target area is scaled up for easier calculation.
  const double small_target_area_threshold =
    total_target_area_with_data * small_area_threshold_frac;

  std::cerr << "Using Small target area threshold: "
            << small_target_area_threshold << std::endl;

  // Check whether target areas exist that are missing or very small
  bool small_target_area_exists = false;
  bool missing_target_area_exists = false;
  for (InsetState &inset_state : inset_states_) {
    for (const auto &gd : inset_state.geo_divs()) {
      const double target_area = inset_state.target_area_at(gd.id());
      inset_state.insert_whether_input_target_area_is_missing(
        gd.id(),
        target_area < 0.0);
      if (target_area < 0.0) {
        missing_target_area_exists = true;
      } else if (target_area <= small_target_area_threshold) {
        small_target_area_exists = true;
      }
    }
  }

  // Deal with target areas that are below the threshold
  if (small_target_area_exists) {
    double replacement_target_area;

    // We replace the zero and small areas, if any, with
    // small_target_area_threshold if not all target areas are initially
    // missing or zero
    if (small_target_area_threshold > 0.0) {
      std::cerr << "Replacing small target areas..." << std::endl;
      replacement_target_area = small_target_area_threshold;
    } else {

      // If all target areas are zero or missing, we assign the minimum GeoDiv
      // area (instead of the minimum target area) as replacement_target_area.
      std::cerr << "No non-zero target area.\n"
                << "Setting zero target areas to the minimum positive area."
                << std::endl;
      double min_positive_area = dbl_inf;
      for (const InsetState &inset_state : inset_states_) {
        for (const auto &gd : inset_state.geo_divs()) {
          min_positive_area = std::min(min_positive_area, gd.area());
        }
      }
      replacement_target_area = min_positive_area;
    }

    // Replace the small target areas
    for (InsetState &inset_state : inset_states_) {
      for (const auto &gd : inset_state.geo_divs()) {

        // Current target area
        const double target_area = inset_state.target_area_at(gd.id());
        if (
          (target_area >= 0.0) &&
          (target_area <= small_target_area_threshold)) {

          // Do not allow the replacement target area to be smaller than the
          // GeoDiv's target area
          double gd_specific_replacement_target_area = std::max(
            std::min(replacement_target_area, gd.area() * mean_density),
            target_area);
          inset_state.replace_target_area(
            gd.id(),
            gd_specific_replacement_target_area);
          std::cerr << gd.id() << ": " << target_area << " to "
                    << gd_specific_replacement_target_area
                    << " Area: " << gd.area() << "\n";

          // Update total target area
          total_target_area_with_data +=
            (gd_specific_replacement_target_area - target_area);
        }
      }
    }
  }

  // Deal with missing target areas
  if (missing_target_area_exists) {

    // Assign new target areas to GeoDivs
    for (InsetState &inset_state : inset_states_) {
      for (const auto &gd : inset_state.geo_divs()) {
        if (inset_state.target_area_is_missing(gd.id())) {
          double new_target_area;

          // If all target areas are missing, make all GeoDivs equal to their
          // geographic area
          if (almost_equal(total_target_area_with_data, 0.0)) {
            new_target_area = gd.area();
          } else {

            // Replace target_area
            const double adjusted_mean_density =
              total_target_area_with_data / total_start_area_with_data;
            new_target_area = adjusted_mean_density * gd.area();
          }
          inset_state.replace_target_area(gd.id(), new_target_area);
        }
      }
    }
  }
}

void CartogramInfo::rescale_insets()
{

  // Iterate over insets and normalize areas
  for (InsetState &inset_state : inset_states_) {
    inset_state.normalize_inset_area(cart_initial_total_target_area());

    // Rescale copy of original map too
    if (args_.redirect_exports_to_stdout) {
      inset_state.normalize_inset_area(
        cart_initial_total_target_area(),
        false,
        true);
    }
  }
}

std::string CartogramInfo::set_map_name(const std::string &map_name)
{
  map_name_ = map_name;
  if (map_name_.find_last_of("/\\") != std::string::npos) {
    map_name_ = map_name_.substr(map_name_.find_last_of("/\\") + 1);
  }
  if (map_name_.find('.') != std::string::npos) {
    map_name_ = map_name_.substr(0, map_name_.find('.'));
  }
  return map_name_;
}

void CartogramInfo::print_time_report()
{
  // Stop main timer
  timer.stop("Total time");

  // Print time report for each inset
  for (const InsetState &inset_state : inset_states_) {
    inset_state.print_time_report();
    if (args_.export_time_report) {
      inset_state.export_time_report();
    }
  }

  // Print Total time
  timer.print_summary_report();
}

void CartogramInfo::write_csv(const std::string &csv_file_name)
{
  // Write a csv file with the current target areas
  std::ofstream out_file_csv;
  out_file_csv.open(csv_file_name + ".csv");
  if (!out_file_csv) {
    std::cerr << "ERROR writing CSV: failed to open " << csv_file_name
              << ".csv" << std::endl;
  }

  // Each vector of strings will represent one row, starting with column names
  std::vector<std::vector<std::string> > csv_rows(1);

  csv_rows[0].push_back(id_header_);
  csv_rows[0].push_back("Target Area");

  // Fill up the rows with the IDs and target areas
  for (const InsetState &inset_state : inset_states_) {
    for (const GeoDiv &gd : inset_state.geo_divs()) {
      csv_rows.push_back(
        {gd.id(), std::to_string(inset_state.target_area_at(gd.id()))});
    }
  }

  // Write to CSV object
  auto writer = csv::make_csv_writer(out_file_csv);
  for (const auto &row : csv_rows) {
    writer << row;
  }

  // Close out_file and exit
  out_file_csv.close();
}

InsetState CartogramInfo::convert_to_inset_state()
{

  InsetState new_inset_state("", args_);

  for (const InsetState &inset_state : inset_states_) {
    for (const auto &geo_div : inset_state.geo_divs()) {
      new_inset_state.push_back(geo_div);
      new_inset_state.insert_color(
        geo_div.id(),
        inset_state.color_at(geo_div.id()));
    }
  }
  return new_inset_state;
}

void CartogramInfo::write_shifted_insets()
{

  // Normalize areas
  for (InsetState &inset_state : inset_states_) {

    // The following condition is not possible because
    // project_to_equal_area should take care of it.
    // if (!args_.output_equal_area_map)
    inset_state.adjust_for_dual_hemisphere();
    inset_state.normalize_inset_area(cart_initial_total_target_area(), true);
  }
  // Shift insets so that they do not overlap
  reposition_insets();

  // Output to GeoJSON
  write_geojson("insets_shifted");
  std::exit(EXIT_SUCCESS);
}

void CartogramInfo::write_svg(const std::string &suffix)
{
  InsetState insets_combined = convert_to_inset_state();
  insets_combined.rescale_map();

  // TODO: Figure out how to add a grid
  // scale_factor = sqrt(scale_factor);
  // double scale_factor = cart_initial_total_target_area() /
  // insets_combined.total_inset_area();

  // Figure out combined name
  std::string inset_names = "";
  for (const InsetState &inset_state : inset_states_) {
    inset_names += inset_state.pos();
  }
  // Only attach inset names to filename if there's more than inset
  if (n_insets() > 1) {
    insets_combined.write_map(
      map_name_ + "_" + inset_names + "_" + suffix,
      false);
  } else {
    insets_combined.write_map(map_name_ + "_" + suffix, false);
  }
}