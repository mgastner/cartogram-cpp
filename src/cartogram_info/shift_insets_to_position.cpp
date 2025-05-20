#include "cartogram_info.hpp"
#include "constants.hpp"

void CartogramInfo::reposition_insets(bool output_to_stdout)
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

    try {
      bboxes.at(inset_pos) = inset_state.bbox(output_to_stdout);
    } catch (const std::out_of_range &e) {
      std::cerr << "ERROR: Key '" << inset_pos << "' not found in bboxes. "
                << "Exception: " << e.what() << std::endl;
      // Re-throw, or return a default value
      throw;
    }
  }

  Bbox bboxC;
  Bbox bboxB;
  Bbox bboxL;
  Bbox bboxT;
  Bbox bboxR;

  try {
    bboxC = bboxes.at("C");
  } catch (const std::out_of_range &e) {
    std::cerr << "ERROR: Key '" << "C" << "' not found in bboxes. "
              << "Exception: " << e.what() << std::endl;
    // Re-throw, or return a default value
    throw;
  }

  try {
    bboxB = bboxes.at("B");
  } catch (const std::out_of_range &e) {
    std::cerr << "ERROR: Key '" << "B" << "' not found in bboxes. "
              << "Exception: " << e.what() << std::endl;
    // Re-throw, or return a default value
    throw;
  }

  try {
    bboxL = bboxes.at("L");
  } catch (const std::out_of_range &e) {
    std::cerr << "ERROR: Key '" << "L" << "' not found in bboxes. "
              << "Exception: " << e.what() << std::endl;
    // Re-throw, or return a default value
    throw;
  }

  try {
    bboxT = bboxes.at("T");
  } catch (const std::out_of_range &e) {
    std::cerr << "ERROR: Key '" << "T" << "' not found in bboxes. "
              << "Exception: " << e.what() << std::endl;
    // Re-throw, or return a default value
    throw;
  }

  try {
    bboxR = bboxes.at("R");
  } catch (const std::out_of_range &e) {
    std::cerr << "ERROR: Key '" << "R" << "' not found in bboxes. "
              << "Exception: " << e.what() << std::endl;
    // Re-throw, or return a default value
    throw;
  }

  // Calculate the width and height of all positioned insets without spacing
  double width = bboxC.xmax() - bboxC.xmin() +
                 bboxL.xmax() - bboxL.xmin() +
                 bboxR.xmax() - bboxR.xmin();

  // Considering edge cases where the width of the insets named "T" or "B"
  // might be greater than the width of "C", "L", "R" insets combined
  width = std::max({
    bboxT.xmax() - bboxT.xmin(),  // width of inset T
    bboxB.xmax() - bboxB.xmin(),  // width of inset B
    width  // width of inset C + L + R
  });

  // Similarly for height instead of width
  double height = bboxC.ymax() - bboxC.ymin() +
                  bboxT.ymax() - bboxT.ymin() +
                  bboxB.ymax() - bboxB.ymin();
  height = std::max({
    bboxR.ymax() - bboxR.ymin(),  // height of inset R
    bboxL.ymax() - bboxL.ymin(),  // height of inset L
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
        {bboxC.xmax(), bboxB.xmax(), bboxT.xmax()});
      x += bboxR.xmax();
      x += inset_spacing;
    } else if (pos == "L") {
      x = std::min(
        {bboxC.xmin(), bboxB.xmin(), bboxT.xmin()});

      // At "L", xmin is negative and lies in the 2nd and 3rd quadrant
      x += bboxL.xmin();
      x -= inset_spacing;
    } else if (pos == "T") {
      y = std::max(
        {bboxC.ymax(), bboxR.ymax(), bboxL.ymax()});
      y += bboxT.ymax();
      y += inset_spacing;
    } else if (pos == "B") {
      y = std::min(
        {bboxC.ymin(), bboxR.ymin(), bboxL.ymin()});

      // At "B", ymin is negative and lies in the 3rd and 4th quadrant
      y += bboxB.ymin();
      y -= inset_spacing;
    }

    // Translating inset according to translation vector calculated above
    const Transformation translate(
      CGAL::TRANSLATION,
      CGAL::Vector_2<Scd>(x, y));

    // Apply translation to all points
    inset_state.transform_points(translate, false);

    if (output_to_stdout) {
      inset_state.transform_points(translate, true);
    }
  }
}
