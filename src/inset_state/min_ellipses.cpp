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
      std::cout << "ellipse has semimajor axis " << ell.semimajor << ",\n";
      std::cout << "            semiminor axis " << ell.semiminor << ",\n";
      std::cout << "            center (" << ell.center.x() << ", "
                << ell.center.y() << "),\n";
      std::cout << "            theta = " << 180 * ell.theta / pi
                << " degrees." << std::endl;
      gd.push_back_ellipse(ell);
    }
  }
}

double InsetState::delta_rho(
  Ellipse ell,
  double r_tilde_sq,
  double rho_p,
  double mean_density,
  double pwh_area)
{
  if (r_tilde_sq >= 4 * XI * XI) {
    return 0.0;
  }
  double xi_to_6 = XI * XI * XI * XI * XI * XI;
  double prefac = ((rho_p - mean_density) * pwh_area) /
                  (16 * pi * ell.semimajor * ell.semiminor * xi_to_6);
  double polynomial = -r_tilde_sq * r_tilde_sq * r_tilde_sq +
                      9 * r_tilde_sq * r_tilde_sq * XI * XI -
                      24 * r_tilde_sq * XI * XI * XI * XI + 16 * xi_to_6;
  return prefac * polynomial;
}

void InsetState::fill_with_ellipse_density(bool plot_density)
{
  std::cout << "In fill_with_ellipse_density()" << std::endl;
  for (unsigned int i = 0; i < lx_; i++) {
    for (unsigned int j = 0; j < ly_; j++) {
      rho_init_(i, j) = 0.0;
    }
  }
  double mean_density = total_target_area() / total_inset_area();
  for (auto gd : geo_divs_) {
    std::cout << gd.id() << std::endl;
    std::cout << gd.n_rings() << std::endl;
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
              rho_init_(
                static_cast<unsigned int>(x_grid),
                static_cast<unsigned int>(y_grid)) +=
                delta_rho(ell, r_tilde_sq, rho_p, mean_density, pwh_area);
            }
          }
        }
      }
    }
  }

  // Scale rho_init_ so that
  // - the minimum is equal to depletion_factor * mean_density, where
  //   0 < depletion_factor < 1, and
  // - the mean is equal to mean_density.

  double rho_min = dbl_inf;
  double rho_mean = 0.0;
  double rho_max = -dbl_inf;
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      rho_min = std::min(rho_min, rho_init_(i, j));
      rho_mean += rho_init_(i, j);
      rho_max = std::max(rho_max, rho_init_(i, j));
    }
  }
  rho_mean /= (lx_ * ly_);
  std::cout << "rho_min = " << rho_min << std::endl;
  std::cout << "rho_mean = " << rho_mean << std::endl;
  std::cout << "rho_max = " << rho_max << std::endl;

  double depletion_factor = 0.95;
  double a = ((1.0 - depletion_factor) * mean_density) / (rho_mean - rho_min);
  double b = mean_density * (depletion_factor * rho_mean - rho_min) /
             (rho_mean - rho_min);
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      rho_init_(i, j) = a * rho_init_(i, j) + b;
    }
  }

  rho_min = dbl_inf;
  rho_mean = 0.0;
  rho_max = -dbl_inf;
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      rho_min = std::min(rho_min, rho_init_(i, j));
      rho_mean += rho_init_(i, j);
      rho_max = std::max(rho_max, rho_init_(i, j));
    }
  }
  rho_mean /= (lx_ * ly_);
  std::cout << "rho_min = " << rho_min << std::endl;
  std::cout << "rho_mean = " << rho_mean << std::endl;
  std::cout << "rho_max = " << rho_max << std::endl;

  if (plot_density) {
    std::string file_name = inset_name_ + "_ellipse_density_" +
                            std::to_string(n_finished_integrations_) + ".eps";
    std::cerr << "Writing " << file_name << std::endl;
    write_density_to_eps(file_name, rho_init_.as_1d_array());
  }
}
