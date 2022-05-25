#include "constants.h"
#include "inset_state.h"

void InsetState::min_ellipses()
{
  for (auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      auto ext_ring = pwh.outer_boundary();
      std::cout << "Building minimum enclosing ellipse" << std::endl;
      Min_ellipse me2(
        ext_ring.vertices_begin(),
        ext_ring.vertices_end(),
        true);
      double a, b, c, d, e, f;

      // Following example at
      // https://doc.cgal.org/latest/Bounding_volumes/classCGAL_1_1Min__ellipse__2.html
      // The order "a, c, b" is deliberate so that it is easier to match the
      // coefficients with those at
      // https://en.wikipedia.org/wiki/Ellipse#General_ellipse
      me2.ellipse().double_coefficients(a, c, b, d, e, f);

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
      std::cout << "Ellipse has coefficients: " << a << ", " << b << ", " << c
                << ", " << d << ", " << e << ", " << f << std::endl;

      double denom = (b * b) - (4 * a * c);
      double fac1 =
        (a * e * e) + (c * d * d) - (b * d * e) + ((b * b - 4 * a * c) * f);
      double inner_sqrt = sqrt(((a - c) * (a - c)) + (b * b));
      Ellipse ell;
      ell.semimajor = -sqrt(2 * fac1 * (a + c + inner_sqrt)) / denom;
      ell.semiminor = -sqrt(2 * fac1 * (a + c - inner_sqrt)) / denom;
      ell.center = Point(
        ((2 * c * d) - (b * e)) / denom,
        ((2 * a * e) - (b * d)) / denom);
      ell.theta = (a < c) ? 0.0 : pi;
      if (b != 0.0) {
        ell.theta = atan((c - a - inner_sqrt) / b);
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
