#include "constants.h"
#include "inset_state.h"

constexpr double xi_sq(4.0);

void InsetState::min_ellipses()
{
  for (auto &gd : geo_divs_) {
    // std::cout << "gd " << gd.id() << std::endl;
    for (const auto &pwh : gd.polygons_with_holes()) {
      auto ext_ring = pwh.outer_boundary();
      Ellipse ell;

      // The minimum ellipse is not uniquely defined if there are fewer than
      // 6 points, see:
      // https://math.stackexchange.com/questions/3063610/how-many-points-are-needed-to-uniquely-define-an-ellipse
      // In that case, we use the minimum circle instead of the minimum
      // ellipse.
      if (ext_ring.size() < 6) {
        Min_circle mc(ext_ring.vertices_begin(), ext_ring.vertices_end());
        ell.center = mc.circle().center();
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

// double delta_rho_of_polygon(
//   Ellipse ell,
//   double r_tilde_sq,
//   double rho_p,
//   double rho_mean,
//   double pwh_area)
// {
//   if (r_tilde_sq >= 4 * xi_sq) {
//     return 0.0;
//   }
//   double xi_to_6 = xi_sq * xi_sq * xi_sq;
//   double prefac = ((rho_p - rho_mean) * pwh_area) /
//                   (16 * pi * ell.semimajor * ell.semiminor * xi_to_6);
//   double postfac = r_tilde_sq - 4 * xi_sq;
//   double polynomial = -(r_tilde_sq - xi_sq) * postfac * postfac;
//   return prefac * polynomial;
// }

double ellipse_density_prefactor(
  Ellipse ell,
  double rho_p,
  double rho_mean,
  double pwh_area,
  double nu)
{
  // return nu *(pi * ell.semimajor * ell.semiminor)* (rho_p - rho_mean) /
  // pwh_area; // possibility
  return nu * pwh_area * (rho_p - rho_mean) /
         (pi * ell.semimajor * ell.semiminor);
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
  Ellipse ell,
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
         (128 * pi * ell.semimajor * ell.semiminor * xi_to_6);
}

void InsetState::calculate_den_prefactor_rho_p_ells()
{
  ell_density_prefactors_.clear();
  rho_p_vec_.clear();
  ells_.clear();
  pwh_areas_.clear();
  double rho_mean = total_target_area() / total_inset_area();
  for (auto &gd : geo_divs_) {
    double rho_p = (target_area_at(gd.id()) / gd.area());
    for (unsigned int pgon = 0; pgon < gd.n_polygons_with_holes(); ++pgon) {
      Polygon_with_holes pwh = gd.polygons_with_holes()[pgon];
      const auto ext_ring = pwh.outer_boundary();
      double pwh_area = ext_ring.area();
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        pwh_area += h->area();
      }
      pwh_areas_.push_back(pwh_area);
      Ellipse ell = gd.min_ellipses()[pgon];
      ells_.push_back(ell);
      ell_density_prefactors_.push_back(
        ellipse_density_prefactor(ell, rho_p, rho_mean, pwh_area, 1.0));
      rho_p_vec_.push_back(rho_p);
    }
  }
}

void InsetState::calculate_rho_flux(
  std::unordered_map<Point, double> &rho_mp,
  std::unordered_map<Point, double> &fluxx_mp,
  std::unordered_map<Point, double> &fluxy_mp,
  const double nu)
{
  for (const auto &[start_pt, curr_pt] : proj_qd_.triangle_transformation) {
    double rho_mean = total_target_area() / total_inset_area();
    double rho =
      rho_mean;  // to avoid division by zero while calculating velocity
    // double rho = 0.0;
    double flux_x = 0.0;
    double flux_y = 0.0;
    for (unsigned int pgn_index = 0; pgn_index < ells_.size(); ++pgn_index) {
      auto ell = ells_[pgn_index];
      auto pwh_area = pwh_areas_[pgn_index];
      auto rho_p = rho_p_vec_[pgn_index];
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
          rho += ell_density_prefactors_[pgn_index] *
                 ellipse_density_polynomial(r_tilde_sq);

          double prefac = ellipse_flux_prefactor(
            ell,
            r_tilde_sq,
            rho_p,
            rho_mean,
            pwh_area,
            nu);
          double flux_tilde_x = prefac * x_tilde;
          double flux_tilde_y = prefac * y_tilde;
          flux_x += ell.semimajor * flux_tilde_x * ell.cos_theta -
                    ell.semiminor * flux_tilde_y * ell.sin_theta;
          flux_y += ell.semimajor * flux_tilde_x * ell.sin_theta +
                    ell.semiminor * flux_tilde_y * ell.cos_theta;
        }
      }
    }
    rho_mp[curr_pt] = rho;
    fluxx_mp[curr_pt] = flux_x;
    fluxy_mp[curr_pt] = flux_y;
  }
}


double calculate_velocity_for_point(
  Point pt,
  char direction,
  double t,
  std::unordered_map<Point, double> &rho_mp,
  std::unordered_map<Point, double> &fluxx_mp,
  std::unordered_map<Point, double> &fluxy_mp)
{
  return (direction == 'x') ? (-fluxx_mp[pt] / rho_mp[pt])
                            : (-fluxy_mp[pt] / rho_mp[pt]);
}

void InsetState::flatten_ellipse_density2()
{
  std::cerr << "In flatten_ellipse_density2()" << std::endl;

  std::unordered_map<Point, double> rho_mp, fluxx_mp, fluxy_mp;
  
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

  // Constants for the numerical integrator
  // const double inc_after_acc = 1.1;
  // const double dec_after_not_acc = 0.75;
  // const double abs_tol = (std::min(lx_, ly_) * 1e-6);

  // Clear previous triangle transformation data
  proj_qd_.triangle_transformation.clear();

  for (const Point &pt : unique_quadtree_corners_) {
    proj_qd_.triangle_transformation.insert_or_assign(pt, pt);
  }

  // TODO: Does the next line adjust rho_mean if the total inset area has
  //       changed in previous runs of flatten_density() or
  //       flatten_ellipse_density()?
  double rho_mean = total_target_area() / total_inset_area();
  std::cout << "rho_mean = " << rho_mean << std::endl;

  // Determine attenuation factor nu that keeps density changes caused by
  // any ellipse within a fraction f of the mean density.
  double f = 0.1;

  calculate_den_prefactor_rho_p_ells();

  double rho_min = *std::min_element(
    ell_density_prefactors_.begin(),
    ell_density_prefactors_.end());
  double rho_max = *std::max_element(
    ell_density_prefactors_.begin(),
    ell_density_prefactors_.end());
  std::cout << "rho_min = " << rho_min << ", rho_max = " << rho_max
            << std::endl;

  // Update the prefactor densities
  double nu = 1.0;
  double acceptable_min = -f * rho_mean;
  double acceptable_max = f * rho_mean;
  if (rho_min < acceptable_min || rho_max > acceptable_max) {
    double nu_min = acceptable_min / rho_min;
    double nu_max = acceptable_max / rho_max;
    if (std::max(nu_min, nu_max) < 1.0) {
      nu = (nu_min < nu_max) ? nu_min : nu_max;
    }
  }
  std::cout << "nu = " << nu << std::endl;
  for (unsigned int pgn_index = 0; pgn_index < ell_density_prefactors_.size();
       ++pgn_index) {
    ell_density_prefactors_[pgn_index] *= nu;
  }

  // Calculate densities
  calculate_rho_flux(rho_mp, fluxx_mp, fluxy_mp, nu);

  // Initial time and step size
  double t = 0.0;
  double delta_t = 1e-2;  // Initial time step.
  unsigned int iter = 0;

  return;
}
