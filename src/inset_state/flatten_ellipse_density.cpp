#include "inset_state.h"

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
      // Calculate density, velocity and flux at (val.x(), val.y())
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
