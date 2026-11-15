#include "constants.hpp"
#include "inset_state.hpp"

void InsetState::rescale_map()
{
  // padding 1 means no padding
  double padding = (args_.world ? 1 : padding_unless_world);
  Bbox bb;
  if (args_.world) {

    // Bbox for Smyth-Craster projection. Equivalent to applying Smyth-Craster
    // projection to -180, -90, 90, 180.
    // TODO: It would be more self-documenting to replace (-2.50663, -1.25331)
    // with point_after_smyth_craster_projection(Point(-180.0, -90.0)) and
    // similarly for the other two bounding box coordinate
    std::cerr << "Rescaling world map..." << std::endl;
    bb = Bbox(-2.50663, -1.25331, 2.50663, 1.25331);

    // With --periodic_x, the sides of the frame must be the antimeridian
    // exactly, so that the vertices where the input cuts regions at 180
    // degrees lie on the seam of the periodic flow
    if (args_.periodic_x) {
      bb = Bbox(
        point_after_smyth_craster_projection(Point(-180.0, 0.0)).x(),
        bb.ymin(),
        point_after_smyth_craster_projection(Point(180.0, 0.0)).x(),
        bb.ymax());
    }
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
    (args_.n_grid_rows_or_cols <= 0) ||
    ((args_.n_grid_rows_or_cols & (~args_.n_grid_rows_or_cols + 1)) !=
     args_.n_grid_rows_or_cols)) {
    std::cerr
      << "ERROR: args_.n_grid_rows_or_cols must be an integer power of 2."
      << std::endl;
    std::exit(15);
  }
  double latt_const;
  if (bb.xmax() - bb.xmin() > bb.ymax() - bb.ymin()) {
    lx = args_.n_grid_rows_or_cols;
    latt_const = (new_xmax - new_xmin) / lx;
    ly = 1 << static_cast<int>(ceil(log2((new_ymax - new_ymin) / latt_const)));
    new_ymax = 0.5 * (bb.ymax() + bb.ymin()) + 0.5 * ly * latt_const;
    new_ymin = 0.5 * (bb.ymax() + bb.ymin()) - 0.5 * ly * latt_const;
  } else {
    ly = args_.n_grid_rows_or_cols;
    latt_const = (new_ymax - new_ymin) / ly;
    lx = 1 << static_cast<int>(ceil(log2((new_xmax - new_xmin) / latt_const)));
    new_xmax = 0.5 * (bb.xmax() + bb.xmin()) + 0.5 * lx * latt_const;
    new_xmin = 0.5 * (bb.xmax() + bb.xmin()) - 0.5 * lx * latt_const;
  }
  latt_const_ = latt_const;
  std::cerr << "Rescaling to " << lx << "-by-" << ly
            << " grid with bounding box\n\t(" << new_xmin << ", " << new_ymin
            << ", " << new_xmax << ", " << new_ymax << ")" << std::endl;
  set_grid_dimensions(lx, ly);

  // With --periodic_x, x is mapped by (x - new_xmin) / latt_const. Because lx
  // is a power of two, latt_const = (new_xmax - new_xmin) / lx is exact, so
  // the antimeridian lands exactly on x = 0 and x = lx. Its vertices become
  // the cut vertices.
  if (args_.periodic_x) {
    const double dlx = lx;
    transform_points([&](const Point &p) {
      return Point(
        on_periodic_grid((p.x() - new_xmin) / latt_const, dlx),
        (p.y() - new_ymin) / latt_const);
    });
    cut_vertices_.clear();
    auto collect = [&](const Polygon &ring) {
      for (const Point &p : ring) {
        const bool on_west_side = !(p.x() < 0.0) && !(p.x() > 0.0);
        const bool on_east_side = !(p.x() < dlx) && !(p.x() > dlx);
        if (on_west_side || on_east_side)
          cut_vertices_.insert(p);
      }
    };
    for (const auto &gd : geo_divs_) {
      for (const auto &pwh : gd.polygons_with_holes()) {
        collect(pwh.outer_boundary());
        for (const auto &h : pwh.holes())
          collect(h);
      }
    }
    std::cerr << "New bounding box: " << bbox() << std::endl;
    std::cerr << cut_vertices_.size() << " vertices on the antimeridian"
              << std::endl;
    return;
  }

  // Rescale and translate all GeoDiv coordinates
  const Transformation translate(
    CGAL::TRANSLATION,
    CGAL::Vector_2<Scd>(-new_xmin, -new_ymin));
  const Transformation scale(CGAL::SCALING, (1.0 / latt_const));
  transform_points(translate);
  transform_points(scale);

  // Transformed Bounding Box:
  std::cerr << "New bounding box: " << bbox() << std::endl;
}

void InsetState::normalize_inset_area(
  double total_cart_target_area,
  bool equal_area,
  bool normalize_original)
{
  const auto bb = bbox(normalize_original);

  // Calculate scale_factor that makes inset areas proportional to their
  // target areas on the cartogram
  const double inset_area_prop =
    equal_area ? 1.0 : initial_target_area() / total_cart_target_area;
  const double scale_factor =
    equal_area ? 1.0
               : sqrt(inset_area_prop / total_inset_area(normalize_original));

  // Rescale and translate all GeoDiv coordinates
  const Transformation translate(
    CGAL::TRANSLATION,
    CGAL::Vector_2<Scd>(
      -(bb.xmin() + bb.xmax()) / 2,
      -(bb.ymin() + bb.ymax()) / 2));
  const Transformation scale(CGAL::SCALING, scale_factor);
  transform_points(translate, normalize_original);
  transform_points(scale, normalize_original);
}
