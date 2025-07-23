#include "cartogram_info.hpp"
#include "constants.hpp"

void CartogramInfo::reposition_insets()
{

  // Warn user about repositoning insets with `--skip_projection` flag
  if (args_.skip_projection && n_insets() > 1) {
    std::cerr << "WARNING: Trying to repostion insets with ";
    if (crs_ == custom_crs) {
      std::cerr << "custom coordinate reference system " << custom_crs << ". ";
    } else {
      std::cerr << "`--skip_projection` flag present. ";
    }
    std::cerr << "This implies that map has already been projected with "
              << "standard parallels based on original, unprojected map. "
              << "Insets may appear skewed. " << std::endl;
  }

  // For simplicity's sake, let us formally insert bounding boxes for
  // all conceivable inset positions
  std::map<std::string, Bbox> bboxes;
  std::string possible_inset_positions[5] = {"C", "B", "L", "T", "R"};
  for (const auto &pos : possible_inset_positions) {
    bboxes.insert(std::pair<std::string, Bbox>(pos, {0, 0, 0, 0}));
  }

  // If the inset actually exists, we get its current bounding box
  for (const InsetState &inset_state : inset_states_) {
    std::string inset_pos = inset_state.pos();
    bboxes.at(inset_pos) = inset_state.bbox(args_.redirect_exports_to_stdout);
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
  for (InsetState &inset_state : inset_states_) {
    std::string inset_pos = inset_state.pos();

    // Assuming X and Y value of translation vector to be 0 to begin with
    double x = 0;
    double y = 0;
    const std::string pos = inset_pos;

    // We only need to modify either x-coordinates or y-coordinates, depending
    // on the inset position
    if (pos == "R") {
      x = std::max(
        {bboxes.at("C").xmax(), bboxes.at("B").xmax(), bboxes.at("T").xmax()});
      x += bboxes.at("R").xmax();
      x += inset_spacing;
    } else if (pos == "L") {
      x = std::min(
        {bboxes.at("C").xmin(), bboxes.at("B").xmin(), bboxes.at("T").xmin()});

      // At "L", xmin is negative and lies in the 2nd and 3rd quadrant
      x += bboxes.at("L").xmin();
      x -= inset_spacing;
    } else if (pos == "T") {
      y = std::max(
        {bboxes.at("C").ymax(), bboxes.at("R").ymax(), bboxes.at("L").ymax()});
      y += bboxes.at("T").ymax();
      y += inset_spacing;
    } else if (pos == "B") {
      y = std::min(
        {bboxes.at("C").ymin(), bboxes.at("R").ymin(), bboxes.at("L").ymin()});

      // At "B", ymin is negative and lies in the 3rd and 4th quadrant
      y += bboxes.at("B").ymin();
      y -= inset_spacing;
    }

    // Translating inset according to translation vector calculated above
    const Transformation translate(
      CGAL::TRANSLATION,
      CGAL::Vector_2<Scd>(x, y));

    // Apply translation to all points
    inset_state.transform_points(translate, false);

    if (args_.redirect_exports_to_stdout) {
      inset_state.transform_points(translate, true);
    }
  }
}
