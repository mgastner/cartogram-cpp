#include "constants.h"
#include "cartogram_info.h"
#include "inset_state.h"

void rescale_map(unsigned int long_grid_side_length,
                 InsetState *inset_state,
                 bool is_world_map)
{
  double padding = (is_world_map ?  1.0 : padding_unless_world);

  // Initialize bounding box of map with bounding box of 0-th
  // Polygon_with_holes in 0-th GeoDiv
  GeoDiv gd0 = inset_state->geo_divs()[0];
  std::vector<Polygon_with_holes> pwhs = gd0.polygons_with_holes();
  CGAL::Bbox_2 bb0 = pwhs[0].bbox();
  double map_xmin = bb0.xmin();
  double map_xmax = bb0.xmax();
  double map_ymin = bb0.ymin();
  double map_ymax = bb0.ymax();

  // Expand bounding box to enclose all GeoDivs
  for (auto gd : inset_state->geo_divs()) {
    for (auto pwh : gd.polygons_with_holes()) {
      CGAL::Bbox_2 bb = pwh.bbox();
      map_xmin = std::min(map_xmin, bb.xmin());
      map_ymin = std::min(map_ymin, bb.ymin());
      map_xmax = std::max(map_xmax, bb.xmax());
      map_ymax = std::max(map_ymax, bb.ymax());
    }
  }

  // Expand bounding box to guarantee a minimum padding
  double new_xmin = 0.5 * ((1.0-padding)*map_xmax + (1.0+padding)*map_xmin);
  double new_xmax = 0.5 * ((1.0+padding)*map_xmax + (1.0-padding)*map_xmin);
  double new_ymin = 0.5 * ((1.0-padding)*map_ymax + (1.0+padding)*map_ymin);
  double new_ymax = 0.5 * ((1.0+padding)*map_ymax + (1.0-padding)*map_ymin);

  // Ensure that the grid dimensions lx and ly are integer powers of 2
  if ((long_grid_side_length <= 0) ||
      ((long_grid_side_length &
        (~long_grid_side_length + 1)) != long_grid_side_length)) {
    std::cerr << "ERROR: long_grid_side_length must be an integer"
              << "power of 2." << std::endl;
    _Exit(15);
  }
  unsigned int lx, ly;
  double latt_const;
  if (map_xmax-map_xmin > map_ymax-map_ymin) {
    lx = long_grid_side_length;
    latt_const = (new_xmax-new_xmin) / lx;
    ly = 1 << ((int) ceil(log2((new_ymax-new_ymin) / latt_const)));
    new_ymax = 0.5*(map_ymax+map_ymin) + 0.5*ly*latt_const;
    new_ymin = 0.5*(map_ymax+map_ymin) - 0.5*ly*latt_const;
  } else {
    ly = long_grid_side_length;
    latt_const = (new_ymax-new_ymin) / ly;
    lx = 1 << ((int) ceil(log2((new_xmax-new_xmin) / latt_const)));
    new_xmax = 0.5*(map_xmax+map_xmin) + 0.5*lx*latt_const;
    new_xmin = 0.5*(map_xmax+map_xmin) - 0.5*lx*latt_const;
  }
  std::cerr << "Rescaling to " << lx << "-by-" << ly
            << " grid with bounding box" << std::endl;
  std::cerr << "\t("
            << new_xmin << ", " << new_ymin << ", "
            << new_xmax << ", " << new_ymax << ")"
            << std::endl;
  inset_state->set_grid_dimensions(lx, ly);

  // Rescale all GeoDiv coordinates
  Transformation translate(CGAL::TRANSLATION,
                           CGAL::Vector_2<Epick>(-new_xmin, -new_ymin));
  Transformation scale(CGAL::SCALING, (1.0/latt_const));
  for (auto &gd : *inset_state->ref_to_geo_divs()) {
    for (auto &pwh : *gd.ref_to_polygons_with_holes()) {
      Polygon *ext_ring = &pwh.outer_boundary();
      *ext_ring = transform(translate, *ext_ring);
      *ext_ring = transform(scale, *ext_ring);
      for (auto hi = pwh.holes_begin(); hi != pwh.holes_end(); ++hi) {
        *hi = transform(translate, *hi);
        *hi = transform(scale, *hi);
      }
    }
  }

  // Storing coordinates to rescale in future
  inset_state->set_xmin(new_xmin);
  inset_state->set_ymin(new_ymin);
  inset_state->set_map_scale(latt_const);
  return;
}

void normalize_inset_area(InsetState *inset_state,
                          double total_target_area, 
                          bool equal_area)
{

  // Initialize bounding box of map with bounding box of 0-th
  // Polygon_with_holes in 0-th GeoDiv
  GeoDiv gd0 = inset_state->geo_divs()[0];
  std::vector<Polygon_with_holes> pwhs = gd0.polygons_with_holes();
  CGAL::Bbox_2 bb0 = pwhs[0].bbox();
  double map_xmin = bb0.xmin();
  double map_xmax = bb0.xmax();
  double map_ymin = bb0.ymin();
  double map_ymax = bb0.ymax();

  // Expand bounding box to enclose all GeoDivs
  for (auto gd : inset_state->geo_divs()) {
    for (auto pwh : gd.polygons_with_holes()) {
      CGAL::Bbox_2 bb = pwh.bbox();
      map_xmin = std::min(map_xmin, bb.xmin());
      map_ymin = std::min(map_ymin, bb.ymin());
      map_xmax = std::max(map_xmax, bb.xmax());
      map_ymax = std::max(map_ymax, bb.ymax());
    }
  }

  // Calculates scale_factor value to make insets proportionate to each other
  double inset_size_proportion =
    inset_state->total_target_area() / total_target_area;
  double scale_factor = sqrt(1.0/inset_state->cart_area()
                        * inset_size_proportion);
  if (equal_area) {
    scale_factor = 1;
  }

  // Rescale all GeoDiv coordinates
  Transformation translate(CGAL::TRANSLATION,
                           CGAL::Vector_2<Epick>(-(map_xmin + map_xmax) / 2,
                                                 -(map_ymin + map_ymax) / 2));
  Transformation scale(CGAL::SCALING, scale_factor);
  for (auto &gd : *inset_state->ref_to_geo_divs()) {
    for (auto &pwh : *gd.ref_to_polygons_with_holes()) {
      Polygon *ext_ring = &pwh.outer_boundary();
      *ext_ring = transform(translate, *ext_ring);
      *ext_ring = transform(scale, *ext_ring);
      for (auto hi = pwh.holes_begin(); hi != pwh.holes_end(); ++hi) {
        *hi = transform(translate, *hi);
        *hi = transform(scale, *hi);
      }
    }
  }
  return;
}

void shift_insets_to_target_position(CartogramInfo *cart_info)
{
  // For simplicity's sake, let us formally insert bounding boxes for
  // all conceivable inset positions
  std::map<std::string, CGAL::Bbox_2> bboxes;
  std::string possible_inset_positions[5] = {"C", "B", "L", "T", "R"};
  for (auto pos : possible_inset_positions) {
    bboxes.insert(std::pair<std::string, CGAL::Bbox_2>(pos, {0, 0, 0, 0}));
  }

  // If the inset actually exists, we get its current bounding box
  for (auto &inset_state : *cart_info->ref_to_inset_states()) {
    bboxes.at(inset_state.pos()) = inset_state.bbox();
  }
  double x = 0;
  double y = 0;
  for (auto &inset_state : *cart_info->ref_to_inset_states()) {
    const std::string pos = inset_state.pos();
    if (pos == "R") {
      x = std::max(bboxes.at("C").xmax(), bboxes.at("B").xmax());
      x = std::max(x, bboxes.at("T").xmax());
      x += bboxes.at("R").xmax();
      y = 0;
    } else if (pos == "L") {
      x = std::min(bboxes.at("C").xmin(), bboxes.at("B").xmin());
      x = std::min(x, bboxes.at("T").xmin());
      x += bboxes.at("L").xmin();
      y = 0;
    } else if (pos == "T") {
      x = 0;
      y = std::max(bboxes.at("C").ymax(), bboxes.at("R").ymax());
      y = std::max(y, bboxes.at("L").ymax());
      y += bboxes.at("T").ymax();
    } else if (pos == "B") {
      x = 0;
      y = std::min(bboxes.at("C").ymin(), bboxes.at("R").ymin());
      y = std::min(y, bboxes.at("L").ymin());
      y += bboxes.at("B").ymin();
    }
    Transformation translate(CGAL::TRANSLATION,
                             CGAL::Vector_2<Epick>(x, y));
    for (auto &gd : *inset_state.ref_to_geo_divs()) {
      for (auto &pwh : *gd.ref_to_polygons_with_holes()) {
        Polygon *ext_ring = &pwh.outer_boundary();
        *ext_ring = transform(translate, *ext_ring);
        for (auto hi = pwh.holes_begin(); hi != pwh.holes_end(); ++hi) {
          *hi = transform(translate, *hi);
        }
      }
    }
  }
  return;
}
