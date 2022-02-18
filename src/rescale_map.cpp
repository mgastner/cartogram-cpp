#include "constants.h"
#include "cartogram_info.h"
#include "inset_state.h"
#include "smyth_projection.h"

void rescale_map(unsigned int max_n_graticule_rows_or_cols,
                 InsetState *inset_state,
                 bool is_world_map)
{
  double padding = (is_world_map ?  1.0 : padding_unless_world);
  Bbox bb;
  if (is_world_map) {
    bb = Bbox(project_x_to_smyth(-180.0),
              project_y_to_smyth(-90.0),
              project_x_to_smyth(180.0),
              project_y_to_smyth(90.0));
  } else {
    bb = inset_state->bbox();
  }

  // Expand bounding box to guarantee a minimum padding
  double new_xmin =
    0.5 * ((1.0-padding)*bb.xmax() + (1.0+padding)*bb.xmin());
  double new_xmax =
    0.5 * ((1.0+padding)*bb.xmax() + (1.0-padding)*bb.xmin());
  double new_ymin =
    0.5 * ((1.0-padding)*bb.ymax() + (1.0+padding)*bb.ymin());
  double new_ymax =
    0.5 * ((1.0+padding)*bb.ymax() + (1.0-padding)*bb.ymin());

  // Ensure that the grid dimensions lx and ly are integer powers of 2
  if ((max_n_graticule_rows_or_cols <= 0) ||
      ((max_n_graticule_rows_or_cols &
        (~max_n_graticule_rows_or_cols + 1)) != max_n_graticule_rows_or_cols)) {
    std::cerr << "ERROR: max_n_graticule_rows_or_cols must be an integer"
              << "power of 2."
              << std::endl;
    _Exit(15);
  }
  unsigned int lx, ly;
  double latt_const;
  if (bb.xmax()-bb.xmin() > bb.ymax()-bb.ymin()) {
    lx = max_n_graticule_rows_or_cols;
    latt_const = (new_xmax-new_xmin) / lx;
    ly = 1 << static_cast<int>(ceil(log2((new_ymax-new_ymin) / latt_const)));
    new_ymax = 0.5*(bb.ymax()+bb.ymin()) + 0.5*ly*latt_const;
    new_ymin = 0.5*(bb.ymax()+bb.ymin()) - 0.5*ly*latt_const;
  } else {
    ly = max_n_graticule_rows_or_cols;
    latt_const = (new_ymax-new_ymin) / ly;
    lx = 1 << static_cast<int>(ceil(log2((new_xmax-new_xmin) / latt_const)));
    new_xmax = 0.5*(bb.xmax()+bb.xmin()) + 0.5*lx*latt_const;
    new_xmin = 0.5*(bb.xmax()+bb.xmin()) - 0.5*lx*latt_const;
  }
  std::cerr << "Rescaling to " << lx << "-by-" << ly
            << " grid with bounding box" << std::endl;
  std::cerr << "\t("
            << new_xmin << ", " << new_ymin << ", "
            << new_xmax << ", " << new_ymax << ")"
            << std::endl;
  inset_state->set_grid_dimensions(lx, ly);

  // Rescale and translate all GeoDiv coordinates
  const Transformation translate(CGAL::TRANSLATION,
                           CGAL::Vector_2<Scd>(-new_xmin, -new_ymin));
  const Transformation scale(CGAL::SCALING, (1.0/latt_const));
  for (auto &gd : *inset_state->ref_to_geo_divs()) {
    for (auto &pwh : *gd.ref_to_polygons_with_holes()) {
      auto *ext_ring = &pwh.outer_boundary();
      *ext_ring = transform(translate, *ext_ring);
      *ext_ring = transform(scale, *ext_ring);
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        *h = transform(translate, *h);
        *h = transform(scale, *h);
      }
    }
  }
  return;
}

void normalize_inset_area(InsetState *inset_state,
                          double total_cart_target_area,
                          bool equal_area)
{
  const auto bb = inset_state->bbox();

  // Calculate scale_factor value to make insets proportional to each other
  const double inset_area_prop =
    inset_state->total_target_area() / total_cart_target_area;
  const double scale_factor =
    equal_area ?
    1.0 :
    sqrt(inset_area_prop / inset_state->total_inset_area());

  // Rescale and translate all GeoDiv coordinates
  const Transformation translate(
    CGAL::TRANSLATION,
    CGAL::Vector_2<Scd>(-(bb.xmin() + bb.xmax()) / 2,
                          -(bb.ymin() + bb.ymax()) / 2));
  const Transformation scale(CGAL::SCALING, scale_factor);
  for (auto &gd : *inset_state->ref_to_geo_divs()) {
    for (auto &pwh : *gd.ref_to_polygons_with_holes()) {
      auto *ext_ring = &pwh.outer_boundary();
      *ext_ring = transform(translate, *ext_ring);
      *ext_ring = transform(scale, *ext_ring);
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        *h = transform(translate, *h);
        *h = transform(scale, *h);
      }
    }
  }
  return;
}

void shift_insets_to_target_position(CartogramInfo *cart_info)
{
  // For simplicity's sake, let us formally insert bounding boxes for
  // all conceivable inset positions
  std::map<std::string, Bbox> bboxes;
  std::string possible_inset_positions[5] = {"C", "B", "L", "T", "R"};
  for (const auto &pos : possible_inset_positions) {
    bboxes.insert(std::pair<std::string, Bbox>(pos, {0, 0, 0, 0}));
  }

  // If the inset actually exists, we get its current bounding box
  for (const auto &[inset_pos, inset_state] :
       *cart_info->ref_to_inset_states()) {
    bboxes.at(inset_pos) = inset_state.bbox();
  }

  // Calculate the width and height of all positioned insets without spacing
  double width = bboxes.at("C").xmax() - bboxes.at("C").xmin() +
                 bboxes.at("L").xmax() - bboxes.at("L").xmin() +
                 bboxes.at("R").xmax() - bboxes.at("R").xmin();

  // Considering edge cases where the width of the insets named "T" or "B"
  // might be greater than the width of "C", "L", "R" insets combined
  width = std::max({
    bboxes.at("T").xmax() - bboxes.at("T").xmin(),  // width of inset T
    bboxes.at("B").xmax() - bboxes.at("B").xmin(),  // width of inset B
    width  // width of inset C + L + R
  });

  // Similarly for height instead of width
  double height = bboxes.at("C").ymax() - bboxes.at("C").ymin() +
                  bboxes.at("T").ymax() - bboxes.at("T").ymin() +
                  bboxes.at("B").ymax() - bboxes.at("B").ymin();
  height = std::max({
    bboxes.at("R").ymax() - bboxes.at("R").ymin(),  // height of inset R
    bboxes.at("L").ymax() - bboxes.at("L").ymin(),  // height of inset L
    height  // height of inset C + T + B
  });

  // Spacing between insets
  const double inset_spacing = std::max(width, height) * inset_spacing_factor;
  for (auto &[inset_pos, inset_state] : *cart_info->ref_to_inset_states()) {

    // Assuming X and Y value of translation vector to be 0 to begin with
    double x = 0;
    double y = 0;
    const std::string pos = inset_pos;

    // We only need to modify either x-coordinates or y-coordinates, depening
    // on the inset position
    if (pos == "R") {
      x = std::max({bboxes.at("C").xmax(),
                    bboxes.at("B").xmax(),
                    bboxes.at("T").xmax()});
      x += bboxes.at("R").xmax();
      x += inset_spacing;
    } else if (pos == "L") {
      x = std::min({bboxes.at("C").xmin(),
                    bboxes.at("B").xmin(),
                    bboxes.at("T").xmin()});

      // At "L", xmin is negative and lies in the 2nd and 3rd quadrant
      x += bboxes.at("L").xmin();
      x -= inset_spacing;
    } else if (pos == "T") {
      y = std::max({bboxes.at("C").ymax(),
                    bboxes.at("R").ymax(),
                    bboxes.at("L").ymax()});
      y += bboxes.at("T").ymax();
      y += inset_spacing;
    } else if (pos == "B") {
      y = std::min({bboxes.at("C").ymin(),
                    bboxes.at("R").ymin(),
                    bboxes.at("L").ymin()});

      // At "B", ymin is negative and lies in the 3rd and 4th quadrant
      y += bboxes.at("B").ymin();
      y -= inset_spacing;
    }

    // Translating inset according to translation vector calculated above
    const Transformation translate(CGAL::TRANSLATION,
                                   CGAL::Vector_2<Scd>(x, y));

    // TODO: WE USE THE for-LOOP BELOW SEVERAL TIMES IN THIS FILE. CAN IT BE
    // TURNED INTO A FUNCTION TO REDUCE CODE DUPLICATION?
    for (auto &gd : *inset_state.ref_to_geo_divs()) {
      for (auto &pwh : *gd.ref_to_polygons_with_holes()) {
        auto *ext_ring = &pwh.outer_boundary();
        *ext_ring = transform(translate, *ext_ring);
        for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
          *h = transform(translate, *h);
        }
      }
    }
  }
  return;
}
