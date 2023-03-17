#include "constants.h"
#include "inset_state.h"

void InsetState::min_ellipses()
{
  for (auto &gd : geo_divs_) {
    gd.clear_min_ellipses();
    for (auto &pwh : gd.polygons_with_holes()) {
      auto ext_ring_actual = pwh.outer_boundary();
      std::vector<double> x_coords;
      std::vector<double> y_coords;
      for (const auto &pt : ext_ring_actual) {
        x_coords.push_back(pt.x());
        y_coords.push_back(pt.y());
      }
      std::nth_element(
        x_coords.begin(),
        x_coords.begin() + x_coords.size() / 2,
        x_coords.end());
      std::nth_element(
        y_coords.begin(),
        y_coords.begin() + y_coords.size() / 2,
        y_coords.end());
        
      double x_median = x_coords[x_coords.size() / 2];
      double y_median = y_coords[y_coords.size() / 2];
      
      // Shift the polygon so that its median is at the origin
      // This is necessary to avoid floating point errors
      Polygon ext_ring;
      for (const auto &v : ext_ring_actual) {
        ext_ring.push_back(Point(v.x() - x_median, v.y() - y_median));
      }
      Ellipse ell;

      // The minimum ellipse is not uniquely defined if there are fewer than
      // 6 points, see:
      // https://math.stackexchange.com/questions/3063610/how-many-points-are-needed-to-uniquely-define-an-ellipse
      // In that case, we use the minimum circle instead of the minimum
      // ellipse.
      if (ext_ring.size() < 6) {
        Min_circle mc(ext_ring.vertices_begin(), ext_ring.vertices_end());
        ell.center = mc.circle().center();
        
        // Shift the center back to the original coordinate system
        ell.center =
          Point(ell.center.x() + x_median, ell.center.y() + y_median);
        ell.semimajor = sqrt(mc.circle().squared_radius());
        ell.semiminor = ell.semimajor;
        ell.cos_theta = 1.0;
        ell.sin_theta = 0.0;
      } else {
        Min_ellipse me(
          ext_ring.vertices_begin(),
          ext_ring.vertices_end(),
          true);
        double a, b, c, d, e, f;

        // Following example at
        // https://doc.cgal.org/latest/Bounding_volumes/classCGAL_1_1Min__ellipse__2.html
        // The order "a, c, b" is deliberate so that it is easier to match the
        // coefficients with those at
        // https://en.wikipedia.org/wiki/Ellipse#General_ellipse
        me.ellipse().double_coefficients(a, c, b, d, e, f);

        // If a < 0, we flip the signs of all coefficients so that we identify
        // correctly which axis is the semimajor axis.
        if (a < 0) {
          a *= -1;
          b *= -1;
          c *= -1;
          d *= -1;
          e *= -1;
          f *= -1;
        }
        double denom = (b * b) - (4 * a * c);
        double fac1 =
          (a * e * e) + (c * d * d) - (b * d * e) + ((b * b - 4 * a * c) * f);
        double inner_sqrt = sqrt(((a - c) * (a - c)) + (b * b));
        ell.semimajor = -sqrt(2 * fac1 * (a + c + inner_sqrt)) / denom;
        ell.semiminor = -sqrt(2 * fac1 * (a + c - inner_sqrt)) / denom;
        ell.center = Point(
          ((2 * c * d) - (b * e)) / denom,
          ((2 * a * e) - (b * d)) / denom);
          
        // Shift the center back to the original coordinate system
        ell.center =
          Point(ell.center.x() + x_median, ell.center.y() + y_median);

        double theta = (a < c) ? 0.0 : pi / 2;
        if (b != 0.0) {
          theta = atan((c - a - inner_sqrt) / b);
        }
        ell.theta = theta;
        ell.cos_theta = cos(theta);
        ell.sin_theta = sin(theta);
      }
      gd.push_back_ellipse(ell);
    }
  }
}

double ellipse_density_prefactor(
  double rho_p,
  double rho_mean,
  double pwh_area,
  double nu)
{
  return nu * pwh_area * (rho_p - rho_mean) / pi;
}

double ellipse_density_polynomial(double r_tilde_sq)
{
  if (r_tilde_sq >= 4 * xi_sq)
    return 0.0;
  return -(
    (r_tilde_sq - xi_sq) * (r_tilde_sq - 4 * xi_sq) *
    (r_tilde_sq - 4 * xi_sq) / (16 * xi_sq * xi_sq * xi_sq));
}

double ellipse_flux_prefactor(
  double r_tilde_sq,
  double rho_p,
  double rho_mean,
  double pwh_area,
  double nu)
{
  if (r_tilde_sq >= 4 * xi_sq)
    return 0.0;
  double xi_to_6 = xi_sq * xi_sq * xi_sq;
  return nu * pwh_area * (rho_p - rho_mean) * (4 * xi_sq - r_tilde_sq) *
         (4 * xi_sq - r_tilde_sq) * (4 * xi_sq - r_tilde_sq) /
         (128 * pi * xi_to_6);
}

bool all_map_points_are_in_domain(
  double delta_t,
  std::unordered_map<Point, Point> &proj_map,
  std::unordered_map<Point, Point> &v_intp,
  const unsigned int lx,
  const unsigned int ly)
{
  // Return false if and only if there exists a point that would be outside
  // [0, lx] x [0, ly]
  for (auto &[key, val] : proj_map) {
    double x = val.x() + 0.5 * delta_t * v_intp[key].x();
    double y = val.y() + 0.5 * delta_t * v_intp[key].y();
    if (x < 0.0 || x > lx || y < 0.0 || y > ly) {
      return false;
    }
  }
  return true;
}

void calculate_velocity(
  std::unordered_map<Point, double> &rho_mp,
  std::unordered_map<Point, Point> &flux_mp,
  std::unordered_map<Point, Point> &triangle_transformation,
  std::unordered_map<Point, Point> &velocity)
{
  for (const auto &[key, val] : triangle_transformation) {
    velocity[key] =
      Point(flux_mp[key].x() / rho_mp[key], flux_mp[key].y() / rho_mp[key]);
  }
}

Point interpolate(Point p, Delaunay &dt, std::unordered_map<Point, Point> &mp)
{
  // Find the triangle containing the point
  Face_handle fh = dt.locate(p);

  // Get three vertices
  Point v1 = fh->vertex(0)->point();
  Point v2 = fh->vertex(1)->point();
  Point v3 = fh->vertex(2)->point();

  // Calculate barycentric coordinates
  std::tuple<Scd::FT, Scd::FT, Scd::FT> bary_coor;
  bary_coor =
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
  auto v1_proj = mp[v1];
  auto v2_proj = mp[v2];
  auto v3_proj = mp[v3];

  // Calculate projected point of p
  const Point p_proj = Point(
    bary_x * v1_proj.x() + bary_y * v2_proj.x() + bary_z * v3_proj.x(),
    bary_x * v1_proj.y() + bary_y * v2_proj.y() + bary_z * v3_proj.y());

  return p_proj;
}

void InsetState::flatten_ellipse_density()
{
  std::cerr << "In flatten_ellipse_density()" << std::endl;

  std::unordered_map<Point, double> rho_mp;
  std::unordered_map<Point, Point> flux_mp;

  // eul[i][j] will be the new position of proj_[i][j] proposed by a simple
  // Euler step: move a full time interval delta_t with the velocity at time t
  // and position (proj_[i][j].x, proj_[i][j].y)
  std::unordered_map<Point, Point> eul;

  // mid[i][j] will be the new displacement proposed by the midpoint
  // method (see comment below for the formula)
  std::unordered_map<Point, Point> mid;

  // (vx_intp, vy_intp) will be the velocity at position (proj_.x, proj_.y) at
  // time t
  std::unordered_map<Point, Point> v_intp;

  // (vx_intp_half, vy_intp_half) will be the velocity at the midpoint
  // (proj_.x + 0.5*delta_t*vx_intp, proj_.y + 0.5*delta_t*vy_intp) at time
  // t + 0.5*delta_t
  std::unordered_map<Point, Point> v_intp_half;

  std::unordered_map<Point, Point> velocity;

  // Constants for the numerical integrator
  const double inc_after_acc = 1.1;
  const double dec_after_not_acc = 0.75;
  const double abs_tol = (std::min(lx_, ly_) * 1e-6);

  // Clear previous triangle transformation data
  proj_qd_.triangle_transformation.clear();

  for (const Point &pt : unique_quadtree_corners_) {
    proj_qd_.triangle_transformation.insert_or_assign(pt, pt);
  }

  // TODO: Does the next line adjust rho_mean if the total inset area has
  //       changed in previous runs of flatten_density() or
  //       flatten_ellipse_density()?
  // double rho_mean = total_target_area() / total_inset_area();
  
  double rho_mean = (1.0 - (initial_area_ / (lx_ * ly_))) /
                    (1.0 - (total_inset_area() / (lx_ * ly_)));
  std::cerr << "rho_mean = " << rho_mean << std::endl;

  // Determine attenuation factor nu that keeps density changes caused by
  // any ellipse within a fraction f of the mean density.
  double nu = 0.3;

  std::vector<Ellipse> ells;
  std::vector<double> ell_density_prefactors;
  std::vector<double> pwh_rhos;
  std::vector<double> pwh_areas;

  for (auto &gd : geo_divs_) {
    double rho_p = (target_area_at(gd.id()) / gd.area());
    for (unsigned int pgon = 0; pgon < gd.n_polygons_with_holes(); ++pgon) {
      pwh_rhos.push_back(rho_p);
    }
  }

  for (auto &gd : geo_divs_) {
    double rho_p = (target_area_at(gd.id()) / gd.area());
    for (unsigned int pgon = 0; pgon < gd.n_polygons_with_holes(); ++pgon) {
      Polygon_with_holes pwh = gd.polygons_with_holes()[pgon];
      const auto ext_ring = pwh.outer_boundary();
      double pwh_area = ext_ring.area();
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        pwh_area += h->area();
      }
      pwh_areas.push_back(pwh_area);
      Ellipse ell = gd.min_ellipses()[pgon];
      ells.push_back(ell);
      ell_density_prefactors.push_back(
        ellipse_density_prefactor(rho_p, rho_mean, pwh_area, nu) /
        (ell.semimajor * ell.semiminor));
    }
  }

  double rho_min = *std::min_element(
    ell_density_prefactors.begin(),
    ell_density_prefactors.end());
  double rho_max = *std::max_element(
    ell_density_prefactors.begin(),
    ell_density_prefactors.end());

  std::cout << "rho_min = " << rho_min << ", rho_max = " << rho_max
            << std::endl;
  std::cout << "rho_mean = " << rho_mean << std::endl;

  // Update the prefactor densities
  double acceptable_max = nu * rho_mean;
  if (rho_max > acceptable_max) {
    double nu_max = acceptable_max / rho_max;
    if (nu_max < 1.0) {
      nu = nu_max;
    }
  }

  std::cerr << "nu = " << nu << std::endl;
  for (unsigned int pgn_index = 0; pgn_index < ell_density_prefactors.size();
       ++pgn_index) {
    ell_density_prefactors[pgn_index] *= nu;
  }

  // Calculate densities
  for (const auto &[start_pt, curr_pt] : proj_qd_.triangle_transformation) {
    double rho =
      rho_mean;  // to avoid division by zero while calculating velocity
    double flux_x = 0.0;
    double flux_y = 0.0;
    for (unsigned int pgn_index = 0; pgn_index < ells.size(); ++pgn_index) {
      auto ell = ells[pgn_index];
      auto pwh_area = pwh_areas[pgn_index];
      auto rho_p = pwh_rhos[pgn_index];
      for (int i = -2; i <= 2; ++i) {
        double x = ((i + abs(i) % 2) * static_cast<int>(lx_)) +
                   (curr_pt.x() * (i % 2 == 0 ? 1 : -1));
        for (int j = -2; j <= 2; ++j) {
          double y = ((j + abs(j) % 2) * static_cast<int>(ly_)) +
                     (curr_pt.y() * (j % 2 == 0 ? 1 : -1));
          double x_tilde = ((x - ell.center.x()) * ell.cos_theta +
                            (y - ell.center.y()) * ell.sin_theta) /
                           ell.semimajor;
          double y_tilde = ((-(x - ell.center.x()) * ell.sin_theta) +
                            (y - ell.center.y()) * ell.cos_theta) /
                           ell.semiminor;
          double r_tilde_sq = (x_tilde * x_tilde) + (y_tilde * y_tilde);
          rho += ell_density_prefactors[pgn_index] *
                 ellipse_density_polynomial(r_tilde_sq);
          double flux_prefac =
            (abs(i) % 2 == 0 ? 1 : -1) * (abs(j) % 2 == 0 ? 1 : -1) *
            ellipse_flux_prefactor(r_tilde_sq, rho_p, rho_mean, pwh_area, nu);
          double flux_tilde_x = flux_prefac * x_tilde * (1 / ell.semiminor);
          double flux_tilde_y = flux_prefac * y_tilde * (1 / ell.semimajor);
          flux_x +=
            flux_tilde_x * ell.cos_theta - flux_tilde_y * ell.sin_theta;
          flux_y +=
            flux_tilde_x * ell.sin_theta + flux_tilde_y * ell.cos_theta;
        }
      }
    }
    rho_mp[curr_pt] = rho;
    flux_mp[curr_pt] = Point(flux_x, flux_y);
  }

  // Initial time and step size
  double t = 0.0;
  double delta_t = 1e-2;  // Initial time step.
  unsigned int iter = 0;

  // write_cairo_map(
  //   inset_name() + "_" + std::to_string(n_finished_integrations()) +
  //   "_flux", false, flux_mp);

  // calculate_velocity(
  //   rho_mp,
  //   flux_mp,
  //   proj_qd_.triangle_transformation,
  //   velocity);
  // write_cairo_map(
  //   inset_name() + "_" + std::to_string(n_finished_integrations()) +
  //     "_velcity",
  //   false,
  //   velocity);

  while (t < 1.0) {
    calculate_velocity(
      rho_mp,
      flux_mp,
      proj_qd_.triangle_transformation,
      velocity);

    // calculating velocity at t by filling v_intp
    for (const auto &[key, val] : proj_qd_.triangle_transformation) {
      Point v_intp_val(interpolate(val, proj_qd_.dt, velocity));
      v_intp.insert_or_assign(key, v_intp_val);
    }
    bool accept = false;
    do {
      // Simple Euler step.
      for (const auto &[key, val] : proj_qd_.triangle_transformation) {
        Point eul_val(
          val.x() + v_intp[key].x() * delta_t,
          val.y() + v_intp[key].y() * delta_t);
        eul.insert_or_assign(key, eul_val);
      }

      calculate_velocity(
        rho_mp,
        flux_mp,
        proj_qd_.triangle_transformation,
        velocity);

      accept = all_map_points_are_in_domain(
        delta_t,
        proj_qd_.triangle_transformation,
        v_intp,
        lx_,
        ly_);

      if (accept) {
        for (const auto &[key, val] : proj_qd_.triangle_transformation) {
          Point n_val = Point(
            val.x() + 0.5 * delta_t * v_intp[key].x(),
            val.y() + 0.5 * delta_t * v_intp[key].y());
          Point v_intp_half_val(interpolate(n_val, proj_qd_.dt, velocity));
          v_intp_half.insert_or_assign(key, v_intp_half_val);
          Point mid_val(
            val.x() + v_intp_half[key].x() * delta_t,
            val.y() + v_intp_half[key].y() * delta_t);
          mid.insert_or_assign(key, mid_val);

          // Do not accept the integration step if the maximum squared
          // difference between the Euler and midpoint proposals exceeds
          // abs_tol. Neither should we accept the integration step if one
          // of the positions wandered out of the domain. If one of these
          // problems occurred, decrease the time step.
          const double sq_dist =
            (mid[key].x() - eul[key].x()) * (mid[key].x() - eul[key].x()) +
            (mid[key].y() - eul[key].y()) * (mid[key].y() - eul[key].y());

          if (
            sq_dist > abs_tol || mid[key].x() < 0.0 || mid[key].x() > lx_ ||
            mid[key].y() < 0.0 || mid[key].y() > ly_) {
            accept = false;
          }
        }
      }

      if (!accept) {
        delta_t *= dec_after_not_acc;
      }
    } while (!accept);
    // Control ouput
    if (iter % 10 == 0) {
      std::cerr << "iter = " << iter << ", t = " << t
                << ", delta_t = " << delta_t << "\n";
    }

    // When we get here, the integration step was accepted
    t += delta_t;
    ++iter;

    // Update the triangle transformation map
    proj_qd_.triangle_transformation = mid;
    delta_t *= inc_after_acc;  // Try a larger step next time
  }
}
