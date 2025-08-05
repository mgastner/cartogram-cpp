#include "inset_state.hpp"
#include "interpolate_bilinearly.hpp"
#include "round_point.hpp"

void InsetState::project()
{

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

  if (args_.simplify) {

    simplify(args_.target_points_per_inset);
  }
  if (args_.plot_intersections) {
    write_intersections_image();
  }
}

static Point interpolate_point_with_barycentric_coordinates(
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
  timer.start("Project");
  auto lambda_bary = [&](Point p1) {
    return interpolate_point_with_barycentric_coordinates(
      p1,
      proj_qd_.dt,
      proj_qd_.triangle_transformation);
  };
  transform_points(lambda_bary);

  if (output_to_stdout) {
    transform_points(lambda_bary, true);
  }
  is_simple(__func__);
  timer.stop("Project");
}

// Apply projection to all points in set
void InsetState::project_point_set(std::unordered_set<Point> &unprojected)
{
  auto lambda_bary = [&dt = proj_qd_.dt,
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
