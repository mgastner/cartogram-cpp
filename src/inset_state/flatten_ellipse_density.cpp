#include "constants.h"
#include "inset_state.h"

#define XI (2)

void InsetState::min_ellipses()
{
  for (auto &gd : geo_divs_) {
    //std::cout << "gd " << gd.id() << std::endl;
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
        double theta = (a < c) ? 0.0 : pi;
        if (b != 0.0) {
          theta = atan((c - a - inner_sqrt) / b);
        }
        ell.cos_theta = cos(theta);
        ell.sin_theta = sin(theta);
      }
      gd.push_back_ellipse(ell);
    }
  }
}

double delta_rho_of_polygon(
  Ellipse ell,
  double r_tilde_sq,
  double rho_p,
  double rho_mean,
  double pwh_area)
{
  double xi_sq = XI * XI;
  if (r_tilde_sq >= 4 * xi_sq) {
    return 0.0;
  }
  double xi_to_6 = xi_sq * xi_sq * xi_sq;
  double prefac = ((rho_p - rho_mean) * pwh_area) /
                  (16 * pi * ell.semimajor * ell.semiminor * xi_to_6);
  double postfac = r_tilde_sq - 4 * xi_sq;
  double polynomial = -(r_tilde_sq - xi_sq) * postfac * postfac;
  return prefac * polynomial;
}

double ellipse_flux_prefactor(
  Ellipse ell,
  double r_tilde_sq,
  double rho_p,
  double rho_mean,
  double pwh_area,
  double nu)
{
  double xi_sq = XI * XI;
  if (r_tilde_sq >= 4 * xi_sq) {
    return 0.0;
  }
  double xi_to_6 = xi_sq * xi_sq * xi_sq;
  return nu * pwh_area * (rho_p - rho_mean) * (4 * xi_sq - r_tilde_sq) *
         (4 * xi_sq - r_tilde_sq) * (4 * xi_sq - r_tilde_sq) /
         (128 * pi * ell.semimajor * ell.semiminor * xi_to_6);
}

void InsetState::flatten_ellipse_density()
{
  std::cerr << "In flatten_ellipse_density()" << std::endl;

  // Constants for the numerical integrator
  const double inc_after_acc = 1.1;
  const double dec_after_not_acc = 0.75;
  const double abs_tol = (std::min(lx_, ly_) * 1e-6);

  // Clear previous triangle transformation data
  proj_qd_.triangle_transformation.clear();

  for (const Point &pt : unique_quadtree_corners_) {
    proj_qd_.triangle_transformation.insert_or_assign(pt, pt);
  }

  // eul[i][j] will be the new position of proj_[i][j] proposed by a simple
  // Euler step: move a full time interval delta_t with the velocity at
  // time t and position (proj_[i][j].x, proj_[i][j].y)
  std::unordered_map<Point, Point> eul;

  // mid[i][j] will be the new displacement proposed by the midpoint
  // method (see comment below for the formula)
  std::unordered_map<Point, Point> mid;

  // Initial time and step size
  double t = 0.0;
  double delta_t = 1e-2;  // Initial time step.
  unsigned int iter = 0;

  // Integrate
  while (t < 1.0) {
    for (const auto &[key, val] : proj_qd_.triangle_transformation) {
      // Calculate density, flux and velocity at (val.x(), val.y())
      return;
    }
  }
  //    calculate_velocity(
  //      t,
  //      grid_fluxx_init,
  //      grid_fluxy_init,
  //      rho_ft_,
  //      rho_init_,
  //      &grid_vx,
  //      &grid_vy,
  //      lx_,
  //      ly_);
  //
  //
  //      // We know, either because of the initialization or because of the
  //      // check at the end of the last iteration, that (proj_.x, proj_.y)
  //      // is inside the rectangle [0, lx_] x [0, ly_]. This fact guarantees
  //      // that interpolate_bilinearly() is given a point that cannot cause
  //      it
  //      // to fail.
  //      Point v_intp_val(
  //        interpolate_bilinearly(val.x(), val.y(), &grid_vx, 'x', lx_, ly_),
  //        interpolate_bilinearly(val.x(), val.y(), &grid_vy, 'y', lx_, ly_));
  //      v_intp.insert_or_assign(key, v_intp_val);
  //    }
  //
  //    bool accept = false;
  //    while (!accept) {
  //
  //      // Simple Euler step.
  //      for (const auto &[key, val] : proj_qd_.triangle_transformation) {
  //        Point eul_val(
  //          val.x() + v_intp[key].x() * delta_t,
  //          val.y() + v_intp[key].y() * delta_t);
  //        eul.insert_or_assign(key, eul_val);
  //      }
  //
  //      // Use "explicit midpoint method"
  //      // x <- x + delta_t * v_x(x + 0.5*delta_t*v_x(x,y,t),
  //      //                        y + 0.5*delta_t*v_y(x,y,t),
  //      //                        t + 0.5*delta_t)
  //      // and similarly for y.
  //      calculate_velocity(
  //        t + 0.5 * delta_t,
  //        grid_fluxx_init,
  //        grid_fluxy_init,
  //        rho_ft_,
  //        rho_init_,
  //        &grid_vx,
  //        &grid_vy,
  //        lx_,
  //        ly_);
  //
  //      // Make sure we do not pass a point outside [0, lx_] x [0, ly_] to
  //      // interpolate_bilinearly(). Otherwise decrease the time step below
  //      and
  //      // try again.
  //      accept = all_map_points_are_in_domain(
  //        delta_t,
  //        &proj_qd_.triangle_transformation,
  //        &v_intp,
  //        lx_,
  //        ly_);
  //      if (accept) {
  //
  //        // Okay, we can run interpolate_bilinearly()
  //        for (const auto &[key, val] : proj_qd_.triangle_transformation) {
  //          Point v_intp_half_val(
  //            interpolate_bilinearly(
  //              val.x() + 0.5 * delta_t * v_intp[key].x(),
  //              val.y() + 0.5 * delta_t * v_intp[key].y(),
  //              &grid_vx,
  //              'x',
  //              lx_,
  //              ly_),
  //            interpolate_bilinearly(
  //              val.x() + 0.5 * delta_t * v_intp[key].x(),
  //              val.y() + 0.5 * delta_t * v_intp[key].y(),
  //              &grid_vy,
  //              'y',
  //              lx_,
  //              ly_));
  //          v_intp_half.insert_or_assign(key, v_intp_half_val);
  //          Point mid_val(
  //            val.x() + v_intp_half[key].x() * delta_t,
  //            val.y() + v_intp_half[key].y() * delta_t);
  //          mid.insert_or_assign(key, mid_val);
  //
  //          // Do not accept the integration step if the maximum squared
  //          // difference between the Euler and midpoint proposals exceeds
  //          // abs_tol. Neither should we accept the integration step if one
  //          // of the positions wandered out of the domain. If one of these
  //          // problems occurred, decrease the time step.
  //          const double sq_dist =
  //            (mid[key].x() - eul[key].x()) * (mid[key].x() - eul[key].x()) +
  //            (mid[key].y() - eul[key].y()) * (mid[key].y() - eul[key].y());
  //          if (
  //            sq_dist > abs_tol || mid[key].x() < 0.0 || mid[key].x() > lx_
  //            || mid[key].y() < 0.0 || mid[key].y() > ly_) { accept = false;
  //          }
  //        }
  //      }
  //      if (!accept) {
  //        delta_t *= dec_after_not_acc;
  //      }
  //    }
  //
  //    // Control ouput
  //    if (iter % 10 == 0) {
  //      std::cerr << "iter = " << iter << ", t = " << t
  //                << ", delta_t = " << delta_t << "\n";
  //    }
  //
  //    // When we get here, the integration step was accepted
  //    t += delta_t;
  //    ++iter;
  //
  //    // Update the triangle transformation map
  //    proj_qd_.triangle_transformation = mid;
  //    delta_t *= inc_after_acc;  // Try a larger step next time
  //  }
  //
  //  // Add current proj to proj_sequence vector
  //  proj_sequence_.push_back(proj_qd_);
  //
  //  // Clean up
  //  grid_fluxx_init.destroy_fftw_plan();
  //  grid_fluxy_init.destroy_fftw_plan();
  //  grid_fluxx_init.free();
  //  grid_fluxy_init.free();

  return;
}
