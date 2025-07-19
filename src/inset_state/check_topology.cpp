#include "inset_state.hpp"

// Returns error if there are holes not inside their respective polygons
void InsetState::holes_inside_polygons() const
{
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto &ext_ring = pwh.outer_boundary();
      for (auto &h : pwh.holes()) {
        for (auto &p : h) {

          // TODO: In the future, a better method would be to only check
          // whether one point in each hole is on the bounded side of the
          // exterior ring. Next, check whether the hole intersects the
          // polygon at any point. For this, the function
          // "do_intersect(Polygon, Polygon)" may help.
          if (ext_ring.bounded_side(p) == CGAL::ON_UNBOUNDED_SIDE) {
            CGAL::set_pretty_mode(std::cerr);
            std::cerr << "ERROR: Hole detected outside polygon!";
            std::cerr << " Hole: " << h;
            std::cerr << ". Polygon: " << ext_ring;
            std::cerr << ". GeoDiv: " << gd.id() << std::endl;
            std::exit(20);
          }
        }
      }
    }
  }
}

void InsetState::is_simple(const char *caller_func) const
{
  if (!args_.simplify)
    return;

  // Only check topology if simplification and densification is enabled.
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      bool okay = true;
      if (!pwh.outer_boundary().is_simple() && okay) {
        std::cerr << "ERROR: Outer boundary is not simple for GeoDiv "
                  << gd.id();
        std::cerr << ". is_simple() called from " << caller_func << std::endl;
        write_map(
          inset_name_ + "_" + std::to_string(n_finished_integrations_) +
            "_not_simple_after_" + caller_func,
          false);
        okay = false;
      }
      for (const auto &h : pwh.holes()) {
        if (!h.is_simple() && okay) {
          std::cerr << "ERROR: Hole is not simple for GeoDiv " << gd.id();
          std::cerr << ". is_simple() called from " << caller_func
                    << std::endl;
          write_map(
            inset_name_ + "_" + std::to_string(n_finished_integrations_) +
              "_not_simple_after_" + caller_func,
            false);
          okay = false;
        }
      }
      if (!okay) {
        std::vector<std::pair<size_t, size_t>> changes = densification_changes;
        double max_ratio = 0.0;
        double avg_ratio = 0.0;
        std::pair<size_t, size_t> worst_result = {0, 0};
        for (size_t i = 0; i < changes.size(); ++i) {
          double ratio =
            static_cast<double>(changes[i].second) / changes[i].first;
          std::cerr << "Iteration: " << i + 1 << ", "
                    << "From: " << changes[i].first << ", "
                    << "To: " << changes[i].second << ", "
                    << "Ratio: " << ratio << std::endl;
          avg_ratio += ratio;
          if (ratio > max_ratio) {
            max_ratio = ratio;
            worst_result = changes[i];
          }
        }
        avg_ratio /= changes.size();
        std::cerr << "Maximum ratio: " << max_ratio << std::endl;
        std::cerr << "Average ratio: " << avg_ratio << std::endl;
        std::string csv_file_name = "densification_changes_";
        if (args_.use_munching_method_for_densification) {
          csv_file_name +=
            "munching_" + std::to_string(args_.munching_segment_length) + "_";
        } else {
          csv_file_name += "delaunay_";
        }
        csv_file_name += ".csv";
        if (!std::filesystem::exists(csv_file_name)) {
          std::ofstream out_file_csv(csv_file_name);
          if (!out_file_csv) {
            std::cerr << "ERROR writing CSV: failed to open " << csv_file_name
                      << std::endl;
          } else {
            out_file_csv << "csv_file,from,to,max_ratio,avg_ratio\n";
          }
        }

        // Append to the CSV file
        std::ofstream out_file_csv(csv_file_name, std::ios_base::app);
        out_file_csv << inset_name() << "," << worst_result.first << ","
                     << worst_result.second << ","
                     << static_cast<double>(worst_result.second) /
                          worst_result.first
                     << "," << avg_ratio << "\n";
        out_file_csv.close();
        exit(1);
      }
    }
  }
}

void InsetState::check_topology() const
{
  holes_inside_polygons();
  is_simple(__func__);
}
