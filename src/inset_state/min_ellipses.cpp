#include "constants.h"
#include "inset_state.h"

#define XI (2)

void InsetState::min_ellipses()
{
  for (auto &gd : geo_divs_) {
    std::cout << "gd " << gd.id() << std::endl;
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
        ell.theta = 0.0;
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
        ell.theta = (a < c) ? 0.0 : pi;
        if (b != 0.0) {
          ell.theta = atan((c - a - inner_sqrt) / b);
        }
      }
      gd.push_back_ellipse(ell);
    }
  }
}

double InsetState::delta_rho(
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

double InsetState::ellipse_flux_prefactor(
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
  return nu * pwh_area * (rho_p - rho_mean) * (4 * xi_sq - r_tilde_sq) /
         (128 * pi * ell.semimajor * ell.semiminor * xi_to_6);
}

void InsetState::fill_with_ellipse_density_and_flux(
  bool plot_density,
  bool plot_flux)
{
  std::cout << "In fill_with_ellipse_density_and_flux()" << std::endl;
  for (unsigned int i = 0; i < lx_; i++) {
    for (unsigned int j = 0; j < ly_; j++) {
      rho_init_(i, j) = 0.0;
    }
  }
  double rho_mean = total_target_area() / total_inset_area();

  std::cout << "rho_mean = " << rho_mean << std::endl;

  for (auto gd : geo_divs_) {
    //    double rho_p = (target_area_at(gd.id()) / total_target_area()) *
    //                   (total_inset_area() / gd.area());

    double rho_p = (target_area_at(gd.id()) / gd.area());

    std::cout << "total_target_area = " << total_target_area() << std::endl;
    std::cout << "total_inset_area = " << total_inset_area() << std::endl;
    std::cout << "target_area = " << target_area_at(gd.id()) << std::endl;
    std::cout << "gd.area = " << gd.area() << std::endl;
    std::cout << "rho_p = " << rho_p << std::endl;

    for (unsigned int pgon = 0; pgon < gd.n_polygons_with_holes(); ++pgon) {
      Ellipse ell = gd.min_ellipses()[pgon];

      std::cout << "ellipse area = " << pi * ell.semimajor * ell.semimajor
                << std::endl;

      Polygon_with_holes pwh = gd.polygons_with_holes()[pgon];
      const auto ext_ring = pwh.outer_boundary();
      double pwh_area = ext_ring.area();
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        pwh_area += h->area();
      }
      double cos_theta = cos(ell.theta);
      double sin_theta = sin(ell.theta);
      for (double x_grid = 0.5; x_grid < lx_; ++x_grid) {
        for (double y_grid = 0.5; y_grid < ly_; ++y_grid) {
          for (int i = -2; i <= 2; ++i) {
            double x = ((i + abs(i) % 2) * static_cast<int>(lx_)) +
                       (x_grid * (i % 2 == 0 ? 1 : -1));
            for (int j = -2; j <= 2; ++j) {
              double y = ((j + abs(j) % 2) * static_cast<int>(ly_)) +
                         (y_grid * (j % 2 == 0 ? 1 : -1));
              double x_tilde = ((x - ell.center.x()) * cos_theta +
                                (y - ell.center.y()) * sin_theta) /
                               ell.semimajor;
              double y_tilde = ((-(x - ell.center.x()) * sin_theta) +
                                (y - ell.center.y()) * cos_theta) /
                               ell.semiminor;
              double r_tilde_sq = (x_tilde * x_tilde) + (y_tilde * y_tilde);
              rho_init_(
                static_cast<unsigned int>(x_grid),
                static_cast<unsigned int>(y_grid)) +=
                delta_rho(ell, r_tilde_sq, rho_p, rho_mean, pwh_area);
            }
          }
        }
      }
    }
  }

  // Scale rho_init_ so that the integration works. For example, negative
  // density must be avoided.
  double rho_min = dbl_inf;
  double rho_max = -dbl_inf;
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      rho_min = std::min(rho_min, rho_init_(i, j));
      rho_max = std::max(rho_max, rho_init_(i, j));
    }
  }

  std::cout << "pos A: rho_min = " << rho_min << ", rho_max = " << rho_max
            << std::endl;

  double acceptable_min = -0.2 * rho_mean;
  double acceptable_max = 0.2 * rho_mean;
  double nu = 1.0;
  if (rho_min < acceptable_min || rho_max > acceptable_max) {
    double nu_min = acceptable_min / rho_min;
    double nu_max = acceptable_max / rho_max;
    if (std::max(nu_min, nu_max) < 1.0) {
      nu = (nu_min < nu_max) ? nu_min : nu_max;
    }
  }
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      rho_init_(i, j) = nu * rho_init_(i, j) + rho_mean;
    }
  }
  if (plot_density) {
    std::string file_name = inset_name_ + "_ellipse_density_" +
                            std::to_string(n_finished_integrations_) + ".eps";
    std::cerr << "Writing " << file_name << std::endl;
    write_density_to_eps(file_name, rho_init_.as_1d_array());
  }

  // Flux
  // Initialize array if running for the first time
  if (ellipse_fluxx_.shape()[0] != lx_ || ellipse_fluxx_.shape()[1] != ly_) {
    ellipse_fluxx_.resize(boost::extents[lx_][ly_]);
  }
  if (ellipse_fluxy_.shape()[0] != lx_ || ellipse_fluxy_.shape()[1] != ly_) {
    ellipse_fluxy_.resize(boost::extents[lx_][ly_]);
  }
  for (unsigned int i = 0; i < lx_; i++) {
    for (unsigned int j = 0; j < ly_; j++) {
      ellipse_fluxx_[i][j] = 0.0;
      ellipse_fluxy_[i][j] = 0.0;
    }
  }
  for (auto gd : geo_divs_) {
    double rho_p = target_area_at(gd.id()) / gd.area();
    for (unsigned int pgon = 0; pgon < gd.n_polygons_with_holes(); ++pgon) {
      Ellipse ell = gd.min_ellipses()[pgon];
      Polygon_with_holes pwh = gd.polygons_with_holes()[pgon];
      const auto ext_ring = pwh.outer_boundary();
      double pwh_area = ext_ring.area();
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        pwh_area += h->area();
      }
      double cos_theta = cos(ell.theta);
      double sin_theta = sin(ell.theta);
      for (double x_grid = 0.5; x_grid < lx_; ++x_grid) {
        for (double y_grid = 0.5; y_grid < ly_; ++y_grid) {
          for (int i = -2; i <= 2; ++i) {
            double x = ((i + abs(i) % 2) * static_cast<int>(lx_)) +
                       (x_grid * (i % 2 == 0 ? 1 : -1));
            for (int j = -2; j <= 2; ++j) {
              double y = ((j + abs(j) % 2) * static_cast<int>(ly_)) +
                         (y_grid * (j % 2 == 0 ? 1 : -1));
              double x_tilde = ((x - ell.center.x()) * cos_theta +
                                (y - ell.center.y()) * sin_theta) /
                               ell.semimajor;
              double y_tilde = ((-(x - ell.center.x()) * sin_theta) +
                                (y - ell.center.y()) * cos_theta) /
                               ell.semiminor;
              double r_tilde_sq = (x_tilde * x_tilde) + (y_tilde * y_tilde);
              double prefac = ellipse_flux_prefactor(
                ell,
                r_tilde_sq,
                rho_p,
                rho_mean,
                pwh_area,
                nu);
              double flux_tilde_x = prefac * x_tilde;
              double flux_tilde_y = prefac * y_tilde;
              double flux_x = ell.semimajor * flux_tilde_x * cos_theta -
                              ell.semiminor * flux_tilde_y * sin_theta;
              double flux_y = ell.semimajor * flux_tilde_x * sin_theta +
                              ell.semiminor * flux_tilde_y * cos_theta;
              ellipse_fluxx_[x_grid][y_grid] += (i % 2 == 0 ? 1 : -1) * flux_x;
              ellipse_fluxy_[x_grid][y_grid] += (j % 2 == 0 ? 1 : -1) * flux_y;
            }
          }
        }
      }
    }
  }
  if (plot_flux) {
    std::string file_name = inset_name_ + "_ellipse_flux_" +
                            std::to_string(n_finished_integrations_) + ".eps";
    std::cerr << "Writing " << file_name << std::endl;
    write_flux_to_eps(file_name, ellipse_fluxx_, ellipse_fluxy_);
  }
}
