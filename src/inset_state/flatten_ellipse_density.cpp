#include "constants.h"
#include "inset_state.h"
#include <cmath>

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

double ellipse_density_prefactor(
  Ellipse ell,
  double rho_p,
  double rho_mean,
  double pwh_area,
  double nu)
{
  return nu * pwh_area * (rho_p - rho_mean) /
         (pi * ell.semimajor * ell.semiminor);
}

double ellipse_density_polynomial(double r_tilde_sq)
{
  return((r_tilde_sq - xi_sq) * (r_tilde_sq - 4 * xi_sq) * (r_tilde_sq - 4 * xi_sq) /
    (16 * xi_sq * xi_sq * xi_sq));
}

double ellipse_flux_prefactor(
  Ellipse ell,
  double r_tilde_sq,
  double rho_p,
  double rho_mean,
  double pwh_area,
  double nu)
{
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

  // TODO: Does the next line adjust rho_mean if the total inset area has
  //   	changed in previous runs of flatten_density() or
  //   	flatten_ellipse_density()?
  double rho_mean = total_target_area() / total_inset_area();
  std::cout << "rho_mean = " << rho_mean << std::endl;

  // Determine attenuation factor nu that keeps density changes caused by
  // any ellipse within a fraction f of the mean density.
  double f = 0.1;
  std::vector<Ellipse> ells;
  std::vector<double> pwh_areas;
  std::vector<double> ell_density_prefactors;
  std::vector<double> rho_p_vec;
  for (auto gd : geo_divs_) {
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
        ellipse_density_prefactor(ell, rho_p, rho_mean, pwh_area, 1.0));
      rho_p_vec.push_back(rho_p);
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

  double nu = 1.0;
  double acceptable_min = -f * rho_mean; // f = 0.1
  double acceptable_max = f * rho_mean;
  if (rho_min < acceptable_min || rho_max > acceptable_max) {
    double nu_min = acceptable_min / rho_min;
    double nu_max = acceptable_max / rho_max;
    if (std::max(nu_min, nu_max) < 1.0) {
      nu = (nu_min < nu_max) ? nu_min : nu_max;
    }
  }
  std::cout << "nu = " << nu << std::endl;

  for (unsigned int pgn_index = 0; pgn_index < ell_density_prefactors.size();
        ++pgn_index) {
    ell_density_prefactors[pgn_index] *= nu;
  }
  rho_min = *std::min_element(
    ell_density_prefactors.begin(),
    ell_density_prefactors.end());
  rho_max = *std::max_element(
    ell_density_prefactors.begin(),
    ell_density_prefactors.end());
  std::cout << "rho_min_post = " << rho_min << ", rho_max_post = " << rho_max
            << std::endl;

  // Initial time and step size
  double t = 0.0;
  double delta_t = 1e-2;  // Initial time step.
  unsigned int iter = 0;  // Counter for the number of iterations
  bool accepted = false;


  // Integrate
  while (t < 1.0) {
    for (const auto &[start_pt, curr_pt] : proj_qd_.triangle_transformation) {
      double rho = rho_mean; 
      double flux_x = 0.0;
      double flux_y = 0.0;

      // Calculate density, flux and velocity at curr_pt and current time
      for (unsigned int pgn_index = 0; pgn_index < ells.size(); ++pgn_index) {
        auto ell = ells[pgn_index];
        auto pwh_area = pwh_areas[pgn_index];
        auto rho_p = rho_p_vec[pgn_index];
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
            rho += (ell_density_prefactors[pgn_index] *
                    ellipse_density_polynomial(r_tilde_sq)) *
                    (1.0 - t);
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
      //Velocity values
      double v_x = flux_x / rho;
      double v_y = flux_y / rho;
      accepted = false;
      while (!accepted) {
        //Euler value

        Point eul_val(
          curr_pt.x() + v_x * delta_t,
          curr_pt.y() + v_y * delta_t);

      	// Explicit midpoint method:
      	// x <- x + delta_t * v_x(x + 0.5*delta_t*v_x(x,y,t),
    	  //                    	y + 0.5*delta_t*v_y(x,y,t),
    	  //                    	t + 0.5*delta_t)
        
        // Calculating the velocity values needed to find the midpoint as per the midpoint method
        double v_x_t = curr_pt.x() + 0.5 * delta_t * v_x; // x value for new velocity
        double v_y_t = curr_pt.y() + 0.5 * delta_t * v_y; // y value for new velocity
        double t_mid = t + 0.5 * delta_t;

        //Carry out same steps as above to find velocity values

        double rho_mid = 0.0;
        double flux_x_mid = 0.0;
        double flux_y_mid = 0.0;
        for (unsigned int pgn_index_mid = 0; pgn_index_mid < ells.size();
              ++pgn_index_mid) {
          auto ell_mid = ells[pgn_index_mid];
          auto pwh_area_mid = pwh_areas[pgn_index_mid];
          auto rho_p_mid = rho_p_vec[pgn_index_mid];
          for (int i = -2; i <= 2; ++i) {
            double x = ((i + abs(i) % 2) * static_cast<int>(lx_)) +
                       (v_x_t * (i % 2 == 0 ? 1 : -1));
            for (int j = -2; j <= 2; ++j) {
              double y = ((j + abs(j) % 2) * static_cast<int>(ly_)) +
                         (v_y_t * (j % 2 == 0 ? 1 : -1));
              double x_tilde_mid =
                ((x - ell_mid.center.x()) * ell_mid.cos_theta +
                 (y - ell_mid.center.y()) * ell_mid.sin_theta) /
                ell_mid.semimajor;
              double y_tilde_mid =
                ((-(x - ell_mid.center.x()) * ell_mid.sin_theta) +
                 (y - ell_mid.center.y()) * ell_mid.cos_theta) /
                ell_mid.semiminor;
              double r_tilde_sq_mid =
                (x_tilde_mid * x_tilde_mid) + (y_tilde_mid * y_tilde_mid);
              rho_mid += (ell_density_prefactors[pgn_index_mid] *
                          ellipse_density_polynomial(r_tilde_sq_mid)) *
                          (1.0 - t_mid);
              double prefac_mid = ellipse_flux_prefactor(
                ell_mid,
                r_tilde_sq_mid,
                rho_p_mid,
                rho_mean,
                pwh_area_mid,
                nu);
              double flux_tilde_x_mid = prefac_mid * x_tilde_mid;
              double flux_tilde_y_mid = prefac_mid * y_tilde_mid;
              flux_x_mid +=
                ell_mid.semimajor * flux_tilde_x_mid * ell_mid.cos_theta -
                ell_mid.semiminor * flux_tilde_y_mid * ell_mid.sin_theta;
              flux_y_mid +=
                ell_mid.semimajor * flux_tilde_x_mid * ell_mid.sin_theta +
                ell_mid.semiminor * flux_tilde_y_mid * ell_mid.cos_theta;
            }
          }
        }
        double v_x_t_mid = flux_x_mid / rho_mid;
        double v_y_t_mid = flux_y_mid / rho_mid;

        //New displacement proposed by the midpoint method
        double mid_x = curr_pt.x() + delta_t * v_x_t_mid;
        double mid_y = curr_pt.y() + delta_t * v_y_t_mid;

        const double sq_dist = (mid_x - eul_val.x()) * (mid_x - eul_val.x()) +
                               (mid_y - eul_val.y()) * (mid_y - eul_val.y());

        //Checking if the sq_dist is within the tolerance and within the x and y limits
        if (sq_dist <= abs_tol && mid_x >= 0 && mid_x <= lx_ && mid_y >= 0 && mid_y <= ly_) {  
          accepted = true;
          proj_qd_.triangle_transformation[start_pt] = Point(mid_x, mid_y);
          delta_t *= inc_after_acc;
          std::cout << "delta: " << start_pt.x() - mid_x << ", " << start_pt.y() - mid_y << std::endl;
        }
        else{
          delta_t *= dec_after_not_acc;
        }
        

      
      }
      std::cout << "delta_t: " << delta_t << std::endl;
    }
    iter++;
    t += delta_t;
    std::cout << "Number of iterations: " << iter << std::endl;
    std::cout << "t = " << t << std::endl;
  }
}
