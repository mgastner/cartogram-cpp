#include "constants.h"
#include "inset_state.h"

void InsetState::rescale_map(
  unsigned int max_n_grid_rows_or_cols,
  bool is_world_map)
{
  // padding 1 means no padding
  double padding = (is_world_map ? 1 : padding_unless_world);
  Bbox bb;
  if (is_world_map) {

    // Bbox for Smyth-Craster projection. Equivalent to applying Smyth-Craster
    // projection to -180, -90, 90, 180.
    // TODO: It would be more self-documenting to replace (-2.50663, -1.25331)
    // with point_after_smyth_craster_projection(Point(-180.0, -90.0)) and
    // similarly for the other two bounding box coordinate
    std::cerr << "Rescaling world map..." << std::endl;
    bb = Bbox(-2.50663, -1.25331, 2.50663, 1.25331);
  } else {
    bb = bbox();
  }

  // Expand bounding box to guarantee a minimum padding
  double new_xmin =
    0.5 * ((1.0 - padding) * bb.xmax() + (1.0 + padding) * bb.xmin());
  double new_xmax =
    0.5 * ((1.0 + padding) * bb.xmax() + (1.0 - padding) * bb.xmin());
  double new_ymin =
    0.5 * ((1.0 - padding) * bb.ymax() + (1.0 + padding) * bb.ymin());
  double new_ymax =
    0.5 * ((1.0 + padding) * bb.ymax() + (1.0 - padding) * bb.ymin());

  // Ensure that the grid dimensions lx and ly are integer powers of 2
  unsigned int lx, ly;
  if (
    (max_n_grid_rows_or_cols <= 0) ||
    ((max_n_grid_rows_or_cols & (~max_n_grid_rows_or_cols + 1)) !=
     max_n_grid_rows_or_cols)) {
    std::cerr
      << "ERROR: max_n_grid_rows_or_cols must be an integer power of 2."
      << std::endl;
    _Exit(15);
  }
  double latt_const;
  if (bb.xmax() - bb.xmin() > bb.ymax() - bb.ymin()) {
    lx = max_n_grid_rows_or_cols;
    latt_const = (new_xmax - new_xmin) / lx;
    ly = 1 << static_cast<int>(ceil(log2((new_ymax - new_ymin) / latt_const)));
    new_ymax = 0.5 * (bb.ymax() + bb.ymin()) + 0.5 * ly * latt_const;
    new_ymin = 0.5 * (bb.ymax() + bb.ymin()) - 0.5 * ly * latt_const;
  } else {
    ly = max_n_grid_rows_or_cols;
    latt_const = (new_ymax - new_ymin) / ly;
    lx = 1 << static_cast<int>(ceil(log2((new_xmax - new_xmin) / latt_const)));
    new_xmax = 0.5 * (bb.xmax() + bb.xmin()) + 0.5 * lx * latt_const;
    new_xmin = 0.5 * (bb.xmax() + bb.xmin()) - 0.5 * lx * latt_const;
  }
  std::cerr << "Rescaling to " << lx << "-by-" << ly
            << " grid with bounding box\n\t(" << new_xmin << ", " << new_ymin
            << ", " << new_xmax << ", " << new_ymax << ")" << std::endl;
  set_grid_dimensions(lx, ly);

  // Rescale and translate all GeoDiv coordinates
  const Transformation translate(
    CGAL::TRANSLATION,
    CGAL::Vector_2<Scd>(-new_xmin, -new_ymin));
  const Transformation scale(CGAL::SCALING, (1.0 / latt_const));
  for (auto &gd : geo_divs_) {
    for (auto &pwh : gd.ref_to_polygons_with_holes()) {
      auto *ext_ring = &pwh.outer_boundary();
      *ext_ring = transform(translate, *ext_ring);
      *ext_ring = transform(scale, *ext_ring);
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        *h = transform(translate, *h);
        *h = transform(scale, *h);
      }
    }
  }

  // Transformed Bounding Box:
  std::cerr << "New bounding box: " << bbox() << std::endl;
}

void InsetState::normalize_inset_area(
  double total_cart_target_area,
  bool equal_area)
{
  const auto bb = bbox();

  // Calculate scale_factor that makes inset areas proportional to their
  // target areas on the cartogram
  const double inset_area_prop = total_target_area() / total_cart_target_area;
  const double scale_factor =
    equal_area ? 1.0 : sqrt(inset_area_prop / total_inset_area());

  // Rescale and translate all GeoDiv coordinates
  const Transformation translate(
    CGAL::TRANSLATION,
    CGAL::Vector_2<Scd>(
      -(bb.xmin() + bb.xmax()) / 2,
      -(bb.ymin() + bb.ymax()) / 2));
  const Transformation scale(CGAL::SCALING, scale_factor);
  for (auto &gd : geo_divs_) {
    for (auto &pwh : gd.ref_to_polygons_with_holes()) {
      auto *ext_ring = &pwh.outer_boundary();
      *ext_ring = transform(translate, *ext_ring);
      *ext_ring = transform(scale, *ext_ring);
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        *h = transform(translate, *h);
        *h = transform(scale, *h);
      }
    }
  }
}
