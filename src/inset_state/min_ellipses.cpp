#include "inset_state.h"

void InsetState::min_ellipses()
{
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      auto ext_ring = pwh.outer_boundary();
      std::cout << "Building minimum enclosing ellipse" << std::endl;
      Min_ellipse me2(
        ext_ring.vertices_begin(),
        ext_ring.vertices_end(),
        true);
      double a, b, c, d, e, f;
      me2.ellipse().double_coefficients(a, c, b, d, e, f);
      std::cout << "ellipse has the equation " << a << " x^2 + " << c
                << " y^2 + " << b << " xy + " << d << " x + " << e << " y + "
                << f << " = 0." << std::endl;
      double denom = (b * b) - (4 * a * c);

      std::cout << "denom = " << denom << std::endl;

      if (me2.is_degenerate()) {
        std::cout << "Ellipse is degenerate" << std::endl;
      } else {
        std::cout << "Ellipse is not degenerate" << std::endl;
      }

      double fac1 =
        (a * e * e) + (c * d * d) - (b * d * e) + ((b * b - 4 * a * c) * f);
      double inner_sqrt = sqrt(((a - c) * (a - c)) + (b * b));
      double little_a = -sqrt(2 * fac1 * (a + c + inner_sqrt)) / denom;
      double little_b = -sqrt(2 * fac1 * (a + c - inner_sqrt)) / denom;
      std::cout << "ellipse has semimajor axis " << little_a << ",\n";
      std::cout << "            semiminor axis " << little_b << std::endl;
    }
  }
}
