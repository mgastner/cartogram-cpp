#define BOOST_TEST_MODULE test_triangulation_locate_uniform
#include "cgal_typedef.hpp"
#include "quadtree_corner.hpp"
#include "round_point.hpp"
#include "triangulation.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/number_utils.h>
#include <algorithm>
#include <array>
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <cstdint>
#include <random>
#include <type_traits>
#include <vector>

using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPoint = EPICK::Point_2;

static constexpr uint32_t GRID_W = 8;
static constexpr uint32_t GRID_H = 8;

// Uniform 8x8 grid quadtree with unit leaves
struct UniformQuadtreeLocator {
  struct Leaf {
    uint32_t x, y, size;
  };

  UniformQuadtreeLocator(uint32_t w = GRID_W, uint32_t h = GRID_H)
      : w_(w), h_(h)
  {
    leaves_.reserve(w_ * h_);
    for (uint32_t y = 0; y < h_; ++y)
      for (uint32_t x = 0; x < w_; ++x)
        leaves_.push_back({x, y, 1});
  }

  const std::vector<Leaf> &leaves() const
  {
    return leaves_;
  }

  Leaf locate(const Point &p) const
  {
    const double xd = CGAL::to_double(p.x());
    const double yd = CGAL::to_double(p.y());
    const int ix = std::max(
      0,
      std::min<int>(
        static_cast<int>(std::floor(xd)),
        static_cast<int>(w_) - 1));
    const int iy = std::max(
      0,
      std::min<int>(
        static_cast<int>(std::floor(yd)),
        static_cast<int>(h_) - 1));
    return {static_cast<uint32_t>(ix), static_cast<uint32_t>(iy), 1};
  }

  uint32_t width() const
  {
    return w_;
  }
  uint32_t height() const
  {
    return h_;
  }

private:
  uint32_t w_, h_;
  std::vector<Leaf> leaves_;
};

// Base projection with grid offset mapping
struct ProjectionBase {
  explicit ProjectionBase(uint32_t w = GRID_W, uint32_t h = GRID_H)
      : w_(w), h_(h)
  {
  }

  // Build only queries midpoints for s>1; we keep this permissive.
  bool is_valid_corner(uint32_t, uint32_t) const
  {
    return true;
  }

  // Map a leaf (x,y) -> bucket index; TRI::build sized buckets generously.
  uint32_t offset(uint32_t x, uint32_t y) const
  {
    return y * w_ + x;
  }

  uint32_t num_unique_corners() const
  {
    return (w_ + 1) * (h_ + 1);
  }

protected:
  uint32_t w_, h_;
};

// Shear in X: (x', y') = (x + k*y, y)
struct ShearXProjection : ProjectionBase {
  explicit ShearXProjection(double k, uint32_t w = GRID_W, uint32_t h = GRID_H)
      : ProjectionBase(w, h), k_(k)
  {
  }

  Point get(uint32_t x, uint32_t y) const
  {
    return Point(
      static_cast<double>(x) + k_ * static_cast<double>(y),
      static_cast<double>(y));
  }

private:
  double k_;
};

// Rotation about center C=(W/2,H/2)
struct RotateProjection : ProjectionBase {
  explicit RotateProjection(
    double theta,
    uint32_t w = GRID_W,
    uint32_t h = GRID_H)
      : ProjectionBase(w, h), theta_(theta), cx_(0.5 * w), cy_(0.5 * h),
        c_(std::cos(theta_)), s_(std::sin(theta_))
  {
  }

  Point get(uint32_t x, uint32_t y) const
  {
    const double X = static_cast<double>(x) - cx_;
    const double Y = static_cast<double>(y) - cy_;
    return Point(c_ * X - s_ * Y + cx_, s_ * X + c_ * Y + cy_);
  }

private:
  double theta_, cx_, cy_, c_, s_;
};

// Radially shrink toward center by factor (1 - alpha), alpha in (0,1)
struct RadialInwardProjection : ProjectionBase {
  explicit RadialInwardProjection(
    double alpha,
    uint32_t w = GRID_W,
    uint32_t h = GRID_H)
      : ProjectionBase(w, h), alpha_(alpha), cx_(0.5 * w), cy_(0.5 * h)
  {
  }

  Point get(uint32_t x, uint32_t y) const
  {
    const double X = static_cast<double>(x) - cx_;
    const double Y = static_cast<double>(y) - cy_;
    const double s = (1.0 - alpha_);
    return Point(s * X + cx_, s * Y + cy_);
  }

private:
  double alpha_, cx_, cy_;
};

// Smooth "wavy" warp (small, orientation-preserving)
struct WavyProjection : ProjectionBase {
  explicit WavyProjection(double eps, uint32_t w = GRID_W, uint32_t h = GRID_H)
      : ProjectionBase(w, h), eps_(eps)
  {
  }

  Point get(uint32_t x, uint32_t y) const
  {
    const double xd = static_cast<double>(x);
    const double yd = static_cast<double>(y);
    return Point(
      xd + eps_ * std::sin(0.7 * yd),
      yd + eps_ * std::sin(0.5 * xd));
  }

private:
  double eps_;
};

// Mirror across vertical axis: (x', y') = (W - x, y) -> orientation flips
struct MirrorXProjection : ProjectionBase {
  explicit MirrorXProjection(uint32_t w = GRID_W, uint32_t h = GRID_H)
      : ProjectionBase(w, h)
  {
  }

  Point get(uint32_t x, uint32_t y) const
  {
    return Point(
      static_cast<double>(w_) - static_cast<double>(x),
      static_cast<double>(y));
  }
};

// Contains test (same robust predicate as used in Triangulation::locate)
template <class TriangulationT>
static bool triangle_contains(
  const typename TriangulationT::Triangle &T,
  const Point &p)
{
  const EPoint A(T.vertices[0].x(), T.vertices[0].y());
  const EPoint B(T.vertices[1].x(), T.vertices[1].y());
  const EPoint C(T.vertices[2].x(), T.vertices[2].y());
  const EPoint P(p.x(), p.y());

  const auto orientation = CGAL::orientation(A, B, C);
  BOOST_REQUIRE(orientation != CGAL::COLLINEAR);

  if (orientation == CGAL::LEFT_TURN) {
    return CGAL::orientation(A, B, P) != CGAL::RIGHT_TURN &&
           CGAL::orientation(B, C, P) != CGAL::RIGHT_TURN &&
           CGAL::orientation(C, A, P) != CGAL::RIGHT_TURN;
  } else {
    return CGAL::orientation(A, B, P) != CGAL::LEFT_TURN &&
           CGAL::orientation(B, C, P) != CGAL::LEFT_TURN &&
           CGAL::orientation(C, A, P) != CGAL::LEFT_TURN;
  }
}

// Verify that a triangle's 3 integer vertices are the corners of a specific
// leaf.
template <class TriangulationT>
static bool triangle_in_leaf(
  const typename TriangulationT::Triangle &T,
  uint32_t leaf_x,
  uint32_t leaf_y)
{
  const uint32_t Xs[2] = {leaf_x, leaf_x + 1};
  const uint32_t Ys[2] = {leaf_y, leaf_y + 1};
  int hit = 0;
  for (size_t i = 0; i < 3; ++i) {
    const uint32_t vx = T.vertices[i].x();
    const uint32_t vy = T.vertices[i].y();
    const bool okx = (vx == Xs[0]) || (vx == Xs[1]);
    const bool oky = (vy == Ys[0]) || (vy == Ys[1]);
    if (okx && oky)
      ++hit;
    else
      return false;
  }
  return hit == 3;
}

// Deterministic sampler of N points in [0,8) x [0,8) (avoid exact border 8)
static std::vector<Point> sample_points(size_t N)
{
  std::vector<Point> out;
  out.reserve(N);
  std::mt19937_64 rng(0xC0FFEEULL);
  std::uniform_real_distribution<double> U(1e-7, 8.0 - 1e-7);
  for (size_t i = 0; i < N; ++i) {
    out.emplace_back(Point(U(rng), U(rng)));
  }
  return out;
}

template <class Proj>
static void run_projection_case(
  const char *name,
  const Proj &proj,
  bool expect_ok)
{
  UniformQuadtreeLocator qt(GRID_W, GRID_H);
  Triangulation<UniformQuadtreeLocator, Proj> tri;

  BOOST_TEST_CONTEXT(name)
  {
    const bool ok = tri.build(&qt, &proj);
    BOOST_CHECK_EQUAL(ok, expect_ok);
    if (!ok)
      return;

    // Triangulating a quad (no midpoints) -> 2 triangles per leaf
    const auto &Ts = tri.triangles();
    BOOST_REQUIRE_EQUAL(Ts.size(), static_cast<size_t>(GRID_W * GRID_H * 2));

    // Centers of all cells must locate to a triangle within the owning leaf.
    for (uint32_t y = 0; y < GRID_H; ++y) {
      for (uint32_t x = 0; x < GRID_W; ++x) {
        const Point p(
          static_cast<double>(x) + 0.5,
          static_cast<double>(y) + 0.5);
        const auto it = tri.locate(p);
        BOOST_REQUIRE(it != Ts.end());

        using TriT = Triangulation<UniformQuadtreeLocator, Proj>;
        BOOST_CHECK(triangle_contains<TriT>(*it, p));

        // Check leaf ownership
        const auto leaf = qt.locate(p);
        BOOST_CHECK(triangle_in_leaf<TriT>(*it, leaf.x, leaf.y));

        // Hit the fast-path cache by locating again; iterator index must
        // match.
        const auto idx1 = static_cast<uint32_t>(std::distance(Ts.begin(), it));
        const auto it2 = tri.locate(p);
        const auto idx2 =
          static_cast<uint32_t>(std::distance(Ts.begin(), it2));
        BOOST_CHECK_EQUAL(idx1, idx2);
      }
    }

    // Deterministic random points test
    const auto pts = sample_points(1000);
    using TriT = Triangulation<UniformQuadtreeLocator, Proj>;
    for (const auto &p : pts) {
      const auto it = tri.locate(p);
      BOOST_REQUIRE(it != Ts.end());
      BOOST_CHECK(triangle_contains<TriT>(*it, p));

      const auto leaf = qt.locate(p);
      BOOST_CHECK(triangle_in_leaf<TriT>(*it, leaf.x, leaf.y));
    }
  }
}

BOOST_AUTO_TEST_CASE(Build_And_Locate_Shear_Positive)
{
  run_projection_case(
    "Shear +X (k=+0.13)",
    ShearXProjection(+0.13),
    /*expect_ok=*/true);
}

BOOST_AUTO_TEST_CASE(Build_And_Locate_Shear_Negative)
{
  run_projection_case(
    "Shear +X (k=-0.21)",
    ShearXProjection(-0.21),
    /*expect_ok=*/true);
}

BOOST_AUTO_TEST_CASE(Build_And_Locate_Rotate)
{
  run_projection_case(
    "Rotate 25.7 deg",
    RotateProjection(M_PI / 7.0),
    /*expect_ok=*/true);
}

BOOST_AUTO_TEST_CASE(Build_And_Locate_Radial_Shrink)
{
  run_projection_case(
    "Radial inward alpha=0.30",
    RadialInwardProjection(0.30),
    /*expect_ok=*/true);
}

BOOST_AUTO_TEST_CASE(Build_And_Locate_Wavy)
{
  run_projection_case(
    "Wavy eps=0.05",
    WavyProjection(0.05),
    /*expect_ok=*/true);
}

BOOST_AUTO_TEST_CASE(Flip_Is_Detected)
{
  // Mirror in X should flip triangle orientation -> build() must return false
  run_projection_case(
    "Mirror X (flip expected)",
    MirrorXProjection(),
    /*expect_ok=*/false);
}
