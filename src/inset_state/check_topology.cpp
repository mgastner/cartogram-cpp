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
  if (!args_.simplify) return;

  // Only check topology if simplification and densification is enabled.
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      if (!pwh.outer_boundary().is_simple()) {
        std::cerr << "ERROR: Outer boundary is not simple for GeoDiv "
                  << gd.id();
        std::cerr << ". is_simple() called from " << caller_func << std::endl;
        write_map(
          inset_name_ + "_" + std::to_string(n_finished_integrations_) +
            "_not_simple_after_" + caller_func,
          false);
        exit(1);
      }
      for (const auto &h : pwh.holes()) {
        if (!h.is_simple()) {
          std::cerr << "ERROR: Hole is not simple for GeoDiv " << gd.id();
          std::cerr << ". is_simple() called from " << caller_func
                    << std::endl;
          write_map(
            inset_name_ + "_" + std::to_string(n_finished_integrations_) +
              "_not_simple_after_" + caller_func,
            false);
          exit(1);
        }
      }
    }
  }
}

void InsetState::check_topology() const
{
  holes_inside_polygons();
  is_simple(__func__);
}
