#include "inset_state.hpp"
#include "interpolate_bilinearly.hpp"
#include "matrix.hpp"
#include "round_point.hpp"

void InsetState::project()
{

  if (args_.qtdt_method) {

    // Update triangulation adding shorter diagonal as constraint for better
    // shape similarity
    update_delaunay_t();
    if (args_.simplify) {
      densify_geo_divs_using_delaunay_t();
    }

    // Plot if requested
    if (args_.plot_quadtree) {
      write_delaunay_triangles(
        file_prefix_ + "c_updated_delaunay_t_after_flatten",
        false);
    }

    // Project using the updated Delaunay triangulation and plot
    project_with_delaunay_t(args_.redirect_exports_to_stdout);

    if (args_.plot_quadtree) {
      write_delaunay_triangles(
        file_prefix_ + "d_projected_with_updated_delaunay_t",
        true);
    }

  } else if (args_.triangulation) {

    // Only densify if we will also simplify later.
    if (args_.simplify) {

      // Choose diagonals that are inside grid cells, then densify.
      fill_grid_diagonals();
      densify_geo_divs();
    }

    // Project with triangulation
    project_with_triangulation();

  } else {

    // Project using bilinear interpolation
    project_with_bilinear_interpolation();
  }

  if (args_.simplify) {

    simplify(args_.target_points_per_inset);
  }
  if (args_.plot_intersections) {
    write_intersections_image();
  }
}

Point interpolate_point_bilinearly(
  const Point p1,
  const boost::multi_array<double, 2> &xdisp,
  const boost::multi_array<double, 2> &ydisp,
  const unsigned int lx,
  const unsigned int ly)
{
  const double intp_x =
    interpolate_bilinearly(p1.x(), p1.y(), xdisp, 'x', lx, ly);
  const double intp_y =
    interpolate_bilinearly(p1.x(), p1.y(), ydisp, 'y', lx, ly);
  return {p1.x() + intp_x, p1.y() + intp_y};
}

void InsetState::project_with_bilinear_interpolation()
{
  timer.start("Project (Bilinear Interpolation)");
  // Calculate displacement from proj array
  boost::multi_array<double, 2> xdisp(boost::extents[lx_][ly_]);
  boost::multi_array<double, 2> ydisp(boost::extents[lx_][ly_]);

#pragma omp parallel for default(none) shared(xdisp, ydisp)
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      xdisp[i][j] = proj_[i][j].x() - i - 0.5;
      ydisp[i][j] = proj_[i][j].y() - j - 0.5;
    }
  }

// Cumulative projection
#pragma omp parallel for default(none) shared(xdisp, ydisp)
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {

      // TODO: Should the interpolation be made on the basis of triangulation?
      // Calculate displacement for cumulative grid coordinates
      const double grid_intp_x = interpolate_bilinearly(
        cum_proj_[i][j].x(),
        cum_proj_[i][j].y(),
        xdisp,
        'x',
        lx_,
        ly_);
      const double grid_intp_y = interpolate_bilinearly(
        cum_proj_[i][j].x(),
        cum_proj_[i][j].y(),
        ydisp,
        'y',
        lx_,
        ly_);

      // Update cumulative grid coordinates
      cum_proj_[i][j] = Point(
        cum_proj_[i][j].x() + grid_intp_x,
        cum_proj_[i][j].y() + grid_intp_y);
    }
  }

  // Specialize (i.e., curry) interpolate_point_bilinearly() such that it only
  // requires one argument (Point p1).
  std::function<Point(Point)> lambda =
    [&xdisp, &ydisp, lx = lx_, ly = ly_](Point p1) {
      return interpolate_point_bilinearly(p1, xdisp, ydisp, lx, ly);
    };

  // Apply "lambda" to all points
  transform_points(lambda);
  is_simple(__func__);
  timer.stop("Project (Bilinear Interpolation)");
}

Point interpolate_point_with_barycentric_coordinates(
  const Point &p,
  const Delaunay &dt,
  const std::unordered_map<Point, Point> &proj_map)
{
  // Find the triangle containing the point
  const Face_handle fh = dt.locate(p);

  // Get the three vertices
  const Point v1 = fh->vertex(0)->point();
  const Point v2 = fh->vertex(1)->point();
  const Point v3 = fh->vertex(2)->point();

  // Calculate barycentric coordinates
  const std::tuple<Scd::FT, Scd::FT, Scd::FT> bary_coor =
    CGAL::Barycentric_coordinates::triangle_coordinates_in_tuple_2<Point>(
      v1,
      v2,
      v3,
      p);

  // Get the barycentric coordinates
  const double bary_x = std::get<0>(bary_coor);
  const double bary_y = std::get<1>(bary_coor);
  const double bary_z = std::get<2>(bary_coor);

  // Get projected vertices
  const Point v1_proj = proj_map.at(v1);
  const Point v2_proj = proj_map.at(v2);
  const Point v3_proj = proj_map.at(v3);

  // Calculate projected point of p
  return {
    bary_x * v1_proj.x() + bary_y * v2_proj.x() + bary_z * v3_proj.x(),
    bary_x * v1_proj.y() + bary_y * v2_proj.y() + bary_z * v3_proj.y()};
}

void InsetState::project_with_delaunay_t(bool output_to_stdout)
{
  timer.start("Project (Delanuay Triangulation)");
  std::function<Point(Point)> lambda_bary =
    [&dt = proj_qd_.dt,
     &proj_map = proj_qd_.triangle_transformation](Point p1) {
      return interpolate_point_with_barycentric_coordinates(p1, dt, proj_map);
    };
  transform_points(lambda_bary);

  if (output_to_stdout) {
    transform_points(lambda_bary, true);
  }
  is_simple(__func__);
  timer.stop("Project (Delanuay Triangulation)");
}

// In chosen_diag() and transformed_triangle(), the input x-coordinates can
// only be 0, lx, or 0.5, 1.5, ..., lx-0.5. A similar rule applies to the
// y-coordinates.
void InsetState::exit_if_not_on_grid_or_edge(const Point p1) const
{
  const double frac_x = p1.x() - std::floor(p1.x());
  const double frac_y = p1.y() - std::floor(p1.y());

  const bool bad_x = !almost_equal(p1.x(), 0.0) &&
                     !almost_equal(p1.x(), lx_) && !almost_equal(frac_x, 0.5);

  const bool bad_y = !almost_equal(p1.y(), 0.0) &&
                     !almost_equal(p1.y(), ly_) && !almost_equal(frac_y, 0.5);

  if (bad_x || bad_y) {
    std::cerr << "ERROR: Invalid input coordinate in triangulation. "
              << "\tpt = (" << p1.x() << ", " << p1.y() << ")" << std::endl;
    std::exit(1);
  }
}

Point InsetState::projected_point(const Point &p1, const bool project_original)
  const
{
  auto &proj = project_original ? cum_proj_ : proj_;

  exit_if_not_on_grid_or_edge(p1);
  const unsigned int proj_x = std::min(
    static_cast<unsigned int>(lx_) - 1,
    static_cast<unsigned int>(p1.x()));
  const unsigned int proj_y = std::min(
    static_cast<unsigned int>(ly_) - 1,
    static_cast<unsigned int>(p1.y()));
  return {
    (almost_equal(p1.x(), 0.0) || almost_equal(p1.x(), lx_))
      ? p1.x()
      : proj[proj_x][proj_y].x(),
    (almost_equal(p1.y(), 0.0) || almost_equal(p1.y(), ly_))
      ? p1.y()
      : proj[proj_x][proj_y].y()};
}

// Apply projection to all points in set
void InsetState::project_point_set(std::unordered_set<Point> &unprojected)
{
  std::function<Point(Point)> lambda_bary =
    [&dt = proj_qd_.dt,
     &proj_map = proj_qd_.triangle_transformation](Point p1) {
      return interpolate_point_with_barycentric_coordinates(p1, dt, proj_map);
    };
  std::unordered_set<Point> projected;
  for (const Point &pt : unprojected) {
    Point pp = lambda_bary(pt);
    projected.insert(pp);
  }
  unprojected = std::move(projected);
}

// TODO: chosen_diag() seems to be more naturally thought of as a boolean
//       than an integer.

// For a grid cell with corners stored in the Point array v, determine
// whether the diagonal from v[0] to v[2] is inside the grid cell. If
// yes, return 0. Otherwise, if the diagonal from v[1] to v[3] is inside the
// grid cell, return 1. If neither of the two diagonals is inside the
// grid cell, then the cell's topology is invalid; thus, we exit with an
// error message.
int InsetState::chosen_diag(
  const Point v[4],
  unsigned int &num_concave,
  const bool project_original) const
{
  // The input v[i].x can only be 0, lx, or 0.5, 1.5, ..., lx-0.5. A similar
  // rule applies to the y-coordinates.
  for (unsigned int i = 0; i < 4; ++i) {
    exit_if_not_on_grid_or_edge(v[i]);
  }

  // Transform the coordinates in v to the corresponding coordinates on the
  // projected grid. If the x-coordinate is 0 or lx, we keep the input. The
  // input v[i].x can only be 0, lx, or 0.5, 1.5, ..., lx-0.5. A similar rule
  // applies to the y-coordinates. This condition is checked in
  // projected_point().
  Point tv[4];
  for (unsigned int i = 0; i < 4; ++i) {
    tv[i] = projected_point(v[i], project_original);
  }

  // Get the two possible midpoints
  const Point midpoint_diag_0(
    (tv[0].x() + tv[2].x()) / 2,
    (tv[0].y() + tv[2].y()) / 2);
  const Point midpoint_diag_1(
    (tv[1].x() + tv[3].x()) / 2,
    (tv[1].y() + tv[3].y()) / 2);

  // Get the transformed grid cell as a polygon
  Polygon trans_grid;
  for (auto &i : tv) {
    trans_grid.push_back(i);
  }

  // Check if grid cell is concave
  if (!trans_grid.is_convex()) {
    num_concave += 1;
  }
  if (trans_grid.bounded_side(midpoint_diag_0) == CGAL::ON_BOUNDED_SIDE) {
    return 0;
  }
  if (trans_grid.bounded_side(midpoint_diag_1) == CGAL::ON_BOUNDED_SIDE) {
    return 1;
  }
  std::cerr << "Invalid grid cell! At\n";
  std::cerr << "(" << tv[0].x() << ", " << tv[0].y() << ")\n";
  std::cerr << "(" << tv[1].x() << ", " << tv[1].y() << ")\n";
  std::cerr << "(" << tv[2].x() << ", " << tv[2].y() << ")\n";
  std::cerr << "(" << tv[3].x() << ", " << tv[3].y() << ")\n";
  std::cerr << "Original: \n";
  std::cerr << "(" << v[0].x() << ", " << v[0].y() << ")\n";
  std::cerr << "(" << v[1].x() << ", " << v[1].y() << ")\n";
  std::cerr << "(" << v[2].x() << ", " << v[2].y() << ")\n";
  std::cerr << "(" << v[3].x() << ", " << v[3].y() << ")\n";
  std::cerr << "i: " << static_cast<unsigned int>(v[0].x())
            << ", j: " << static_cast<unsigned int>(v[0].y()) << std::endl;
  exit(1);
}

void InsetState::fill_grid_diagonals(const bool project_original)
{
  timer.start("Densification (using Grid Diagonals)");
  // Initialize array if running for the first time
  if (grid_diagonals_.shape()[0] != lx_ || grid_diagonals_.shape()[1] != ly_) {
    grid_diagonals_.resize(boost::extents[lx_ - 1][ly_ - 1]);
  }
  unsigned int n_concave = 0;  // Count concave grid cells

#pragma omp parallel for default(none) shared(n_concave, project_original)
  for (unsigned int i = 0; i < lx_ - 1; ++i) {
    for (unsigned int j = 0; j < ly_ - 1; ++j) {
      Point v[4];
      v[0] = Point(double(i) + 0.5, double(j) + 0.5);
      v[1] = Point(double(i) + 1.5, double(j) + 0.5);
      v[2] = Point(double(i) + 1.5, double(j) + 1.5);
      v[3] = Point(double(i) + 0.5, double(j) + 1.5);
      grid_diagonals_[i][j] = chosen_diag(v, n_concave, project_original);
    }
  }
  std::cerr << "Number of concave grid cells: " << n_concave << std::endl;
  timer.stop("Densification (using Grid Diagonals)");
}

std::array<Point, 3> InsetState::transformed_triangle(
  const std::array<Point, 3> &tri,
  const bool project_original) const
{
  std::array<Point, 3> transf_tri;
  for (unsigned int i = 0; i < 3; ++i) {
    exit_if_not_on_grid_or_edge(tri[i]);
    const auto transf_pt = projected_point(tri[i], project_original);
    transf_tri[i] = transf_pt;
  }
  return transf_tri;
}

// Determine if a point `pt` is on the boundary of a triangle by using cross
// products to find areas spanned by pt and each triangle edge. Idea from:
// https://stackoverflow.com/questions/7050186/find-if-point-lies-on-line-segment
// This function is needed because, sometimes,
// `triangle.bounded_side(Point(x, y)) == CGAL::ON_BOUNDARY` does not return
// `true` even if the point is on the boundary.
bool is_on_triangle_boundary(const Point &pt, const Polygon &triangle)
{
  for (unsigned int i = 0; i < triangle.size(); ++i) {
    const auto t1 = triangle[i];
    const auto t2 = triangle[(i == triangle.size() - 1) ? 0 : i + 1];
    const double area = (t1.x() - pt.x()) * (t2.y() - pt.y()) -
                        (t2.x() - pt.x()) * (t1.y() - pt.y());
    if (almost_equal(area, 0.0)) {
      return true;
    }
  }
  return false;
}

// Get the untransformed coordinates of the triangle in which the point `pt`
// is located. After transformation, this triangle must be entirely inside
// the transformed grid cell.
std::array<Point, 3> InsetState::untransformed_triangle(
  const Point &pt,
  const bool project_original) const
{
  if (pt.x() < 0 || pt.x() > lx_ || pt.y() < 0 || pt.y() > ly_) {
    CGAL::set_pretty_mode(std::cerr);
    std::cerr << "ERROR: coordinate outside bounding box in " << __func__
              << "(). pt = " << pt << std::endl;
    exit(1);
  }

  // Get original grid coordinates
  Point v[4];
  v[0] = Point(
    std::max(0.0, floor(pt.x() + 0.5) - 0.5),
    std::max(0.0, floor(pt.y() + 0.5) - 0.5));
  v[1] = Point(
    std::min(static_cast<double>(lx_), floor(pt.x() + 0.5) + 0.5),
    v[0].y());
  v[2] = Point(
    v[1].x(),
    std::min(static_cast<double>(ly_), floor(pt.y() + 0.5) + 0.5));
  v[3] = Point(v[0].x(), v[2].y());

  // TODO: diag SEEMS TO BE MORE NATURALLY THOUGHT OF AS bool INSTEAD OF int.
  // Assuming that the transformed grid does not have self-intersections,
  // at least one of the diagonals must be completely inside the grid.
  // We use that diagonal to split the grid into two triangles.
  int diag;
  if (
    almost_equal(v[0].x(), 0.0) || almost_equal(v[0].y(), 0.0) ||
    almost_equal(v[2].x(), lx_) || almost_equal(v[2].y(), ly_)) {

    // Case when the grid is on the edge of the grid.
    // We calculate the chosen diagonal because grid_diagonals_ does not
    // store the diagonals for edge grid cells.
    unsigned int n_concave = 0;
    diag = chosen_diag(v, n_concave, project_original);
  } else {

    // Case when the grid is not on the edge of the grid. We can find the
    // already computed chosen diagonal in grid_diagonals_.
    const auto x = static_cast<unsigned int>(v[0].x());
    const auto y = static_cast<unsigned int>(v[0].y());
    diag = grid_diagonals_[x][y];
  }

  // Get the two possible triangles
  Polygon triangle1;
  Polygon triangle2;
  if (diag == 0) {
    triangle1.push_back(v[0]);
    triangle1.push_back(v[1]);
    triangle1.push_back(v[2]);
    triangle2.push_back(v[0]);
    triangle2.push_back(v[2]);
    triangle2.push_back(v[3]);
  } else {
    triangle1.push_back(v[0]);
    triangle1.push_back(v[1]);
    triangle1.push_back(v[3]);
    triangle2.push_back(v[1]);
    triangle2.push_back(v[2]);
    triangle2.push_back(v[3]);
  }

  // Determine which untransformed triangle the given point is in. If the
  // point is in neither, an error is raised.
  std::array<Point, 3> triangle_coordinates;
  if (
    (triangle1.bounded_side(pt) == CGAL::ON_BOUNDED_SIDE) ||
    is_on_triangle_boundary(pt, triangle1)) {
    for (unsigned int i = 0; i < triangle1.size(); ++i) {
      triangle_coordinates[i] = triangle1[i];
    }
  } else if (
    (triangle2.bounded_side(pt) == CGAL::ON_BOUNDED_SIDE) ||
    is_on_triangle_boundary(pt, triangle2)) {
    for (unsigned int i = 0; i < triangle2.size(); ++i) {
      triangle_coordinates[i] = triangle2[i];
    }
  } else {
    std::cerr << "Point not in grid cell!\n";
    std::cerr << "Point coordinates:\n";
    std::cerr << "(" << pt.x() << ", " << pt.y() << ")\n";
    std::cerr << "Original grid cell:\n";
    std::cerr << "(" << v[0].x() << ", " << v[0].y() << ")\n";
    std::cerr << "(" << v[1].x() << ", " << v[1].y() << ")\n";
    std::cerr << "(" << v[2].x() << ", " << v[2].y() << ")\n";
    std::cerr << "(" << v[3].x() << ", " << v[3].y() << ")\n";
    std::cerr << "Chosen diagonal: " << diag << "\n";
    exit(1);
  }
  return triangle_coordinates;
}

Point affine_trans(
  const std::array<Point, 3> &tri,
  const std::array<Point, 3> &org_tri,
  const Point &pt)
{
  // For each point, we make the following transformation. Suppose we find
  // that, before the cartogram transformation, a point (x, y) is in the
  // triangle (a, b, c). We want to find its position in the projected
  // triangle (p, q, r). We locally approximate the cartogram transformation
  // by an affine transformation T such that T(a) = p, T(b) = q and T(c) = r.
  // We can think of T as a 3x3 matrix
  //    -----------
  //   |t11 t12 t13|
  //   |t21 t22 t23|  such that
  //   | 0   0   1 |
  //    -----------
  //    -----------   ----------     ----------
  //   |t11 t12 t13| | a1 b1 c1 |   | p1 q1 r1 |
  //   |t21 t22 t23| | a2 b2 c2 | = | p2 q2 r2 | or TA = P.
  //   | 0   0   1 | | 1  1  1  |   |  1  1  1 |
  //    -----------   ----------     ----------
  // Hence, T = PA^{-1}.
  //                              -----------------------
  //                             |b2-c2 c1-b1 b1*c2-b2*c1|
  // We have A^{-1} = (1/det(A)) |c2-a2 a1-c1 a2*c1-a1*c2|. By multiplying
  //                             |a2-b2 b1-a1 a1*b2-a2*b1|
  //                              -----------------------
  // PA^{-1} we obtain t11, t12, t13, t21, t22, and t23. If the original
  // coordinates are (x, y) on the unprojected map, then the transformed
  // coordinates are:
  // post.x = t11*x + t12*y + t13, post.y = t21*x + t22*y + t23.
  const Point pre(pt.x(), pt.y());

  // Old triangle (a, b, c) expressed as matrix A
  const Matrix abc_mA(org_tri[0], org_tri[1], org_tri[2]);

  // New triangle (p, q, r) expressed as matrix P
  const Matrix pqr_mP(tri[0], tri[1], tri[2]);

  // Transformation matrix T
  const auto mT = pqr_mP.multiplied_with(abc_mA.inverse());

  // Transformed point
  return mT.transformed_point(pre);
}

Point InsetState::projected_point_with_triangulation(
  const Point &pt,
  const bool project_original) const
{
  // Get the untransformed triangle the point pt is in
  const auto old_triangle = untransformed_triangle(pt, project_original);

  // Get the coordinates of the transformed triangle
  const auto new_triangle =
    transformed_triangle(old_triangle, project_original);

  // Get the transformed point and return it
  const auto transformed_pt = affine_trans(new_triangle, old_triangle, pt);
  return rounded_point(transformed_pt, lx_, ly_);
}

void InsetState::project_with_triangulation()
{
  timer.start("Project (Triangulation)");
  // Store reference to current object and call member function
  // projected_point_with_triangulation
  // https://www.nextptr.com/tutorial/ta1430524603/
  // capture-this-in-lambda-expression-timeline-of-change
  std::function<Point(Point)> lambda = [&](Point p1) {
    return projected_point_with_triangulation(p1);
  };

  // Transforming all points based on triangulation
  transform_points(lambda);

  // Cumulative projection
#pragma omp parallel for default(none)
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      cum_proj_[i][j] = projected_point_with_triangulation(cum_proj_[i][j]);
    }
  }
  is_simple(__func__);
  timer.stop("Project (Triangulation)");
}

void InsetState::project_with_cum_proj()
{
  std::function<Point(Point)> lambda = [&](Point p1) {
    return projected_point_with_triangulation(p1, true);
  };

  // Transforming all points based on triangulation
  transform_points(lambda, true);
  is_simple(__func__);
}
