// Adapter between the cartogram pipeline and the conformality core: seeds
// tracers of the cumulative map T0, moves them through every successful
// flattening pass with the same piecewise-affine map that projects the
// polygons, and reports diagnostics of the recorded map.

#include "conformality/objective.hpp"
#include "conformality/optimizer.hpp"
#include "conformality/parallel.hpp"
#include "inset_state.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

using conformality::det;
using conformality::frobenius_sq;
using conformality::Mat2;
using conformality::NodeSet;
using conformality::Vec2;

namespace
{

// R2 low-discrepancy sequence (Roberts 2018), based on the plastic number.
// Offsets of the form frac(0.5 + k alpha) avoid the lattice diagonals on
// which symmetric offsets would sit.
constexpr double plastic = 1.32471795724474602596;
constexpr double r2_alpha_x = 1.0 / plastic;
constexpr double r2_alpha_y = 1.0 / (plastic * plastic);

double fractional_part(const double v)
{
  return v - std::floor(v);
}

// q^2 nodes in each unit cell of [0, lx] x [0, ly], at R2 offsets that
// continue the sequence from index `first_index`. T is the identity.
NodeSet make_r2_nodes(
  const unsigned int lx,
  const unsigned int ly,
  const unsigned int q,
  const std::size_t first_index)
{
  const std::size_t per_cell = static_cast<std::size_t>(q) * q;
  const std::size_t n = static_cast<std::size_t>(lx) * ly * per_cell;
  NodeSet nodes;
  nodes.x.reserve(n);
  std::size_t k = first_index;
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      for (std::size_t l = 0; l < per_cell; ++l) {
        const double kd = static_cast<double>(k++);
        nodes.x.push_back(
          {static_cast<double>(i) + fractional_part(0.5 + kd * r2_alpha_x),
           static_cast<double>(j) + fractional_part(0.5 + kd * r2_alpha_y)});
      }
    }
  }
  nodes.w.assign(n, 1.0 / static_cast<double>(per_cell));
  nodes.g0.assign(n, 1.0);
  nodes.y = nodes.x;
  nodes.jac.assign(n, conformality::identity_mat2());
  return nodes;
}

// Even-odd point-in-polygon tests against all GeoDivs, with polygon edges
// bucketed by lattice row so that each query only scans nearby edges
class RegionLocator
{
public:
  RegionLocator(const std::vector<GeoDiv> &geo_divs, const unsigned int ly)
      : rows_(ly), parity_(geo_divs.size(), 0)
  {
    for (std::size_t g = 0; g < geo_divs.size(); ++g) {
      for (const auto &pwh : geo_divs[g].polygons_with_holes()) {
        add_ring(pwh.outer_boundary(), static_cast<int>(g));
        for (const auto &hole : pwh.holes()) {
          add_ring(hole, static_cast<int>(g));
        }
      }
    }
  }

  // Index of the GeoDiv containing (px, py), or -1. Points inside more than
  // one GeoDiv (overlapping input) are counted in n_conflicts.
  int region_of(const double px, const double py)
  {
    if (rows_.empty()) {
      return -1;
    }
    const double r_max = static_cast<double>(rows_.size() - 1);
    const auto r =
      static_cast<std::size_t>(std::clamp(std::floor(py), 0.0, r_max));
    for (const std::size_t e : rows_[r]) {
      const Edge &edge = edges_[e];

      // Half-open rule: the edge crosses the horizontal line through py if
      // exactly one endpoint lies strictly above it
      if ((edge.ay > py) != (edge.by > py)) {
        const double xi =
          edge.ax + (py - edge.ay) * (edge.bx - edge.ax) / (edge.by - edge.ay);
        if (px < xi) {
          const auto g = static_cast<std::size_t>(edge.g);
          if (parity_[g] == 0) {
            touched_.push_back(g);
          }
          parity_[g] = static_cast<unsigned char>(parity_[g] ^ 1u);
        }
      }
    }
    int result = -1;
    for (const std::size_t g : touched_) {
      if (parity_[g] != 0) {
        if (result < 0) {
          result = static_cast<int>(g);
        } else {
          ++n_conflicts;
        }
      }
      parity_[g] = 0;
    }
    touched_.clear();
    return result;
  }

  std::size_t n_conflicts{0};

private:
  struct Edge {
    double ax, ay, bx, by;
    int g;
  };

  void add_ring(const Polygon &ring, const int g)
  {
    const std::vector<Point> pts(ring.begin(), ring.end());
    const std::size_t n = pts.size();
    if (n < 3 || rows_.empty()) {
      return;
    }
    const double r_max = static_cast<double>(rows_.size() - 1);
    for (std::size_t i = 0; i < n; ++i) {
      const Point &a = pts[i];
      const Point &b = pts[(i + 1) % n];

      // Horizontal edges never satisfy the crossing test
      if (!(a.y() < b.y()) && !(b.y() < a.y())) {
        continue;
      }
      const std::size_t e = edges_.size();
      edges_.push_back({a.x(), a.y(), b.x(), b.y(), g});
      const auto r0 = static_cast<std::size_t>(
        std::clamp(std::floor(std::min(a.y(), b.y())), 0.0, r_max));
      const auto r1 = static_cast<std::size_t>(
        std::clamp(std::floor(std::max(a.y(), b.y())), 0.0, r_max));
      for (std::size_t r = r0; r <= r1; ++r) {
        rows_[r].push_back(e);
      }
    }
  }

  std::vector<Edge> edges_;
  std::vector<std::vector<std::size_t>> rows_;
  std::vector<unsigned char> parity_;
  std::vector<std::size_t> touched_;
};

// Summary statistics of one node cloud, used in the report
struct CloudSummary {
  std::size_t n{0};
  double min_det{0.0};
  std::size_t n_non_positive{0};
  double total_det_rel_error{0.0};
  double energy{0.0};
  double energy_land{0.0};
  double max_region_rel_error{0.0};
  double mean_region_rel_error{0.0};
  int worst_region{-1};
};

CloudSummary summarize_cloud(
  const NodeSet &nodes,
  const std::vector<int> &region,
  const std::vector<GeoDiv> &geo_divs,
  const double domain_area)
{
  CloudSummary s;
  s.n = nodes.size();
  s.min_det = std::numeric_limits<double>::infinity();
  double total = 0.0;
  std::vector<double> region_sum(geo_divs.size(), 0.0);
  NodeSet land;
  for (std::size_t j = 0; j < nodes.size(); ++j) {
    const double d = det(nodes.jac[j]);
    s.min_det = std::min(s.min_det, d);
    if (!(d > 0.0)) {
      ++s.n_non_positive;
    }
    total += nodes.w[j] * d;
    if (region[j] >= 0) {
      region_sum[static_cast<std::size_t>(region[j])] += nodes.w[j] * d;
      land.x.push_back(nodes.x[j]);
      land.w.push_back(nodes.w[j]);
      land.g0.push_back(nodes.g0[j]);
      land.y.push_back(nodes.y[j]);
      land.jac.push_back(nodes.jac[j]);
    }
  }
  s.total_det_rel_error = total / domain_area - 1.0;
  s.energy = conformality::objective(nodes);
  s.energy_land = conformality::objective(land);

  // Area-weighted mean and maximum of the relative error of the quadrature
  // of det J over each region, compared with the area of the output polygon
  double weighted = 0.0;
  double total_area = 0.0;
  for (std::size_t g = 0; g < geo_divs.size(); ++g) {
    const double area = geo_divs[g].area();
    if (!(area > 0.0)) {
      continue;
    }
    const double rel = region_sum[g] / area - 1.0;
    weighted += area * std::abs(rel);
    total_area += area;
    if (std::abs(rel) > s.max_region_rel_error) {
      s.max_region_rel_error = std::abs(rel);
      s.worst_region = static_cast<int>(g);
    }
  }
  s.mean_region_rel_error = total_area > 0.0 ? weighted / total_area : 0.0;
  return s;
}

nlohmann::json cloud_to_json(
  const CloudSummary &s,
  const std::vector<GeoDiv> &geo_divs,
  const unsigned int quad)
{
  return {
    {"nodes_per_cell_axis", quad},
    {"n_nodes", s.n},
    {"min_det_J", s.min_det},
    {"n_non_positive_det", s.n_non_positive},
    {"sum_w_det_J_over_area_minus_1", s.total_det_rel_error},
    {"E0", s.energy},
    {"E0_land", s.energy_land},
    {"region_area_rel_error_max", s.max_region_rel_error},
    {"region_area_rel_error_area_weighted_mean", s.mean_region_rel_error},
    {"region_area_worst",
     s.worst_region >= 0
       ? geo_divs[static_cast<std::size_t>(s.worst_region)].id()
       : std::string()}};
}

// Distance of boundary tracers from their side, how far they moved beyond
// the side's extent, and how far the corners moved
struct BoundaryCheck {
  double max_side_offset{0.0};
  double max_along_overshoot{0.0};
  double max_corner_motion{0.0};
};

BoundaryCheck check_boundary(
  const ConformalityTracers &conf,
  const double lx,
  const double ly)
{
  BoundaryCheck b;
  for (std::size_t k = 0; k < conf.boundary.size(); ++k) {
    const Vec2 p = conf.boundary[k];
    double offset = 0.0;
    double along = 0.0;
    double extent = 0.0;
    switch (conf.boundary_side[k]) {
    case 0:
      offset = std::abs(p.x);
      along = p.y;
      extent = ly;
      break;
    case 1:
      offset = std::abs(p.x - lx);
      along = p.y;
      extent = ly;
      break;
    case 2:
      offset = std::abs(p.y);
      along = p.x;
      extent = lx;
      break;
    default:
      offset = std::abs(p.y - ly);
      along = p.x;
      extent = lx;
      break;
    }
    b.max_side_offset = std::max(b.max_side_offset, offset);
    b.max_along_overshoot =
      std::max({b.max_along_overshoot, -along, along - extent});
  }
  const std::vector<Vec2> targets =
    {{0.0, 0.0}, {lx, 0.0}, {0.0, ly}, {lx, ly}};
  for (std::size_t k = 0; k < conf.corners.size(); ++k) {
    b.max_corner_motion = std::max(
      b.max_corner_motion,
      std::hypot(
        conf.corners[k].x - targets[k].x,
        conf.corners[k].y - targets[k].y));
  }
  return b;
}

nlohmann::json boundary_to_json(const BoundaryCheck &b, const std::size_t n)
{
  return {
    {"n_tracers", n},
    {"max_offset_from_side", b.max_side_offset},
    {"max_overshoot_along_side", b.max_along_overshoot},
    {"max_corner_motion", b.max_corner_motion}};
}

// Composed Jacobians of the optimization nodes against finite differences
// of their companion tracers.
//
// Forward differences with step delta and delta/2 are combined by
// Richardson extrapolation, J_R = 2 J(delta/2) - J(delta), which is exact
// for the piecewise-affine T0 and has an O(delta^2) error after the smooth
// correction A. A pair whose companions fell into a different affine piece
// of T0 in some pass keeps an O(1) discrepancy and is counted as crossing.
constexpr double fd_crossing_threshold = 1e-6;

struct FdCheck {
  std::size_t n_within{0};
  std::size_t n_crossing{0};

  // Among pairs within the threshold: relative error of J_R, and of the
  // plain forward difference with step delta
  double richardson_median{0.0};
  double richardson_max{0.0};
  double forward_median{0.0};
  double forward_max{0.0};
};

FdCheck check_finite_differences(const ConformalityTracers &conf)
{
  auto rel_error = [](const Mat2 &a, const Mat2 &b) {
    const Mat2 diff{
      a.a11 - b.a11,
      a.a12 - b.a12,
      a.a21 - b.a21,
      a.a22 - b.a22};
    return std::sqrt(frobenius_sq(diff) / frobenius_sq(b));
  };
  auto forward = [](const Vec2 y, const Vec2 e1, const Vec2 e2, double d) {
    return Mat2{
      (e1.x - y.x) / d,
      (e2.x - y.x) / d,
      (e1.y - y.y) / d,
      (e2.y - y.y) / d};
  };
  FdCheck c;
  std::vector<double> richardson, plain;
  for (std::size_t k = 0; k < conf.fd_parent.size(); ++k) {
    const std::size_t j = conf.fd_parent[k];
    const Vec2 y = conf.nodes.y[j];
    const Mat2 &jac = conf.nodes.jac[j];
    const Mat2 full = forward(y, conf.fd_e1[k], conf.fd_e2[k], conf.fd_delta);
    const Mat2 half =
      forward(y, conf.fd_e1h[k], conf.fd_e2h[k], 0.5 * conf.fd_delta);
    const Mat2 extrapolated{
      2.0 * half.a11 - full.a11,
      2.0 * half.a12 - full.a12,
      2.0 * half.a21 - full.a21,
      2.0 * half.a22 - full.a22};
    const double rel = rel_error(extrapolated, jac);
    if (rel > fd_crossing_threshold) {
      ++c.n_crossing;
    } else {
      richardson.push_back(rel);
      plain.push_back(rel_error(full, jac));
    }
  }
  std::sort(richardson.begin(), richardson.end());
  std::sort(plain.begin(), plain.end());
  c.n_within = richardson.size();
  if (!richardson.empty()) {
    c.richardson_median = richardson[richardson.size() / 2];
    c.richardson_max = richardson.back();
    c.forward_median = plain[plain.size() / 2];
    c.forward_max = plain.back();
  }
  return c;
}

nlohmann::json fd_to_json(const FdCheck &c, const ConformalityTracers &conf)
{
  return {
    {"delta", conf.fd_delta},
    {"n_pairs", conf.fd_parent.size()},
    {"n_within_threshold", c.n_within},
    {"threshold", fd_crossing_threshold},
    {"richardson_median_rel_error", c.richardson_median},
    {"richardson_max_rel_error", c.richardson_max},
    {"forward_median_rel_error", c.forward_median},
    {"forward_max_rel_error", c.forward_max},
    {"n_crossing", c.n_crossing}};
}

// Distortion measures of Gastner, Seguy and More (2018), Table 1, over the
// nodes (uniform in the original map, as in that table) or over land only
struct DistortionStats {
  std::size_t n{0};
  double e_mean{0.0};
  double e_max{0.0};
  double e_p50{0.0};
  double e_p90{0.0};
  double e_p99{0.0};
  double angular_mean{0.0};
  double angular_max{0.0};
};

DistortionStats distortion_stats(
  const NodeSet &nodes,
  const std::vector<int> &region,
  const bool land_only)
{
  DistortionStats d;
  std::vector<double> e;
  e.reserve(nodes.size());
  double angular_sum = 0.0;
  for (std::size_t j = 0; j < nodes.size(); ++j) {
    if (land_only && region[j] < 0) {
      continue;
    }
    const double ej = conformality::log_axis_ratio(nodes.jac[j]);
    const double aj = conformality::max_angular_distortion(nodes.jac[j]);
    e.push_back(ej);
    angular_sum += aj;
    d.angular_max = std::max(d.angular_max, aj);
  }
  d.n = e.size();
  if (e.empty()) {
    return d;
  }
  double sum = 0.0;
  for (const double v : e) {
    sum += v;
  }
  d.e_mean = sum / static_cast<double>(e.size());
  d.angular_mean = angular_sum / static_cast<double>(e.size());
  std::sort(e.begin(), e.end());
  auto quantile = [&](const double p) {
    const auto i =
      static_cast<std::size_t>(p * static_cast<double>(e.size() - 1) + 0.5);
    return e[std::min(i, e.size() - 1)];
  };
  d.e_p50 = quantile(0.5);
  d.e_p90 = quantile(0.9);
  d.e_p99 = quantile(0.99);
  d.e_max = e.back();
  return d;
}

nlohmann::json distortion_to_json(const DistortionStats &d)
{
  return {
    {"n_nodes", d.n},
    {"e_mean", d.e_mean},
    {"e_max", d.e_max},
    {"e_p50", d.e_p50},
    {"e_p90", d.e_p90},
    {"e_p99", d.e_p99},
    {"angular_mean", d.angular_mean},
    {"angular_max", d.angular_max}};
}

}  // namespace

void InsetState::seed_conformality_tracers()
{
  timer.start("Conformality tracers");
  conf_ = ConformalityTracers{};
  conf_.active = true;
  conf_.lx0 = lx_;
  conf_.ly0 = ly_;
  conf_.quad = args_.conf_quad;
  conf_.eval_quad = args_.conf_eval_quad;

  // The evaluation cloud continues the R2 sequence far beyond the indices
  // used by the optimization nodes, so the two clouds share no points
  conf_.nodes = make_r2_nodes(lx_, ly_, conf_.quad, 1);
  conf_.eval_nodes = make_r2_nodes(lx_, ly_, conf_.eval_quad, 1u << 30);

  RegionLocator locator(geo_divs_, ly_);
  conf_.region.resize(conf_.nodes.size());
  for (std::size_t j = 0; j < conf_.nodes.size(); ++j) {
    conf_.region[j] =
      locator.region_of(conf_.nodes.x[j].x, conf_.nodes.x[j].y);
  }
  conf_.eval_region.resize(conf_.eval_nodes.size());
  for (std::size_t j = 0; j < conf_.eval_nodes.size(); ++j) {
    conf_.eval_region[j] =
      locator.region_of(conf_.eval_nodes.x[j].x, conf_.eval_nodes.x[j].y);
  }
  if (locator.n_conflicts > 0) {
    std::cerr << "WARNING: " << locator.n_conflicts
              << " conformality nodes lie in more than one GeoDiv."
              << std::endl;
  }

  // Boundary tracers at the midpoints of the lattice edges on each side
  const double lx = static_cast<double>(lx_);
  const double ly = static_cast<double>(ly_);
  for (unsigned int j = 0; j < ly_; ++j) {
    const double t = static_cast<double>(j) + 0.5;
    conf_.boundary.push_back({0.0, t});
    conf_.boundary_side.push_back(0);
    conf_.boundary.push_back({lx, t});
    conf_.boundary_side.push_back(1);
  }
  for (unsigned int i = 0; i < lx_; ++i) {
    const double t = static_cast<double>(i) + 0.5;
    conf_.boundary.push_back({t, 0.0});
    conf_.boundary_side.push_back(2);
    conf_.boundary.push_back({t, ly});
    conf_.boundary_side.push_back(3);
  }
  conf_.corners = {{0.0, 0.0}, {lx, 0.0}, {0.0, ly}, {lx, ly}};

  // Finite-difference companions for about 1000 evenly spaced nodes
  const std::size_t stride =
    std::max<std::size_t>(1, conf_.nodes.size() / 1000);
  for (std::size_t j = 0; j < conf_.nodes.size(); j += stride) {
    const Vec2 x = conf_.nodes.x[j];
    conf_.fd_parent.push_back(j);
    conf_.fd_e1.push_back({x.x + conf_.fd_delta, x.y});
    conf_.fd_e2.push_back({x.x, x.y + conf_.fd_delta});
    conf_.fd_e1h.push_back({x.x + 0.5 * conf_.fd_delta, x.y});
    conf_.fd_e2h.push_back({x.x, x.y + 0.5 * conf_.fd_delta});
  }

  // Square graticule: lines of constant x and of constant y
  const double spacing = conf_.grid_spacing;
  const double step = 1.0 / static_cast<double>(conf_.grid_samples_per_unit);
  auto sampled_line = [&](const bool vertical, const double c) {
    std::vector<Vec2> line;
    const double extent = vertical ? ly : lx;
    const auto n = static_cast<std::size_t>(std::lround(extent / step));
    for (std::size_t i = 0; i <= n; ++i) {
      const double t = std::min(extent, static_cast<double>(i) * step);
      line.push_back(vertical ? Vec2{c, t} : Vec2{t, c});
    }
    return line;
  };
  for (double c = 0.0; c <= lx + 1e-9; c += spacing) {
    conf_.grid_lines.push_back(sampled_line(true, std::min(c, lx)));
  }
  for (double c = 0.0; c <= ly + 1e-9; c += spacing) {
    conf_.grid_lines.push_back(sampled_line(false, std::min(c, ly)));
  }

  // Input polygons, for figures of the original map
  for (std::size_t g = 0; g < geo_divs_.size(); ++g) {
    for (const auto &pwh : geo_divs_[g].polygons_with_holes()) {
      auto add_ring = [&](const Polygon &ring) {
        std::vector<Vec2> pts;
        for (const auto &p : ring) {
          pts.push_back({p.x(), p.y()});
        }
        conf_.input_rings.push_back(std::move(pts));
        conf_.input_ring_region.push_back(static_cast<int>(g));
      };
      add_ring(pwh.outer_boundary());
      for (const auto &hole : pwh.holes()) {
        add_ring(hole);
      }
    }
  }

  std::size_t n_land = 0;
  for (const int r : conf_.region) {
    n_land += (r >= 0) ? 1 : 0;
  }
  std::cerr << "Seeded conformality tracers: " << conf_.nodes.size()
            << " nodes (" << n_land << " on land), " << conf_.eval_nodes.size()
            << " evaluation nodes, " << conf_.boundary.size()
            << " boundary tracers, " << conf_.fd_parent.size()
            << " finite-difference pairs" << std::endl;
  timer.stop("Conformality tracers");
}

void InsetState::update_conformality_tracers()
{
  timer.start("Conformality tracers");
  const double x_max = std::nextafter(static_cast<double>(lx_), 0.0);
  const double y_max = std::nextafter(static_cast<double>(ly_), 0.0);

  // Applies the affine map of the triangle that contains p. Only the lookup
  // point is clamped into the lattice; p itself is never moved.
  auto map_point = [&](const Vec2 p, Mat2 *jac) -> Vec2 {
    const Point query(
      std::clamp(p.x, 0.0, x_max),
      std::clamp(p.y, 0.0, y_max));
    const auto it = triang_.locate(query);
    if (it == triang_.triangles().cend()) {
      throw std::runtime_error("conformality tracer outside triangulation");
    }
    const auto &v = it->vertices;
    const double p0x = v[0].x(), p0y = v[0].y();
    const Point q0 = proj_data_.get(v[0].x(), v[0].y());
    const Point q1 = proj_data_.get(v[1].x(), v[1].y());
    const Point q2 = proj_data_.get(v[2].x(), v[2].y());

    // J = [q1 - q0, q2 - q0] [p1 - p0, p2 - p0]^(-1)
    const Mat2 pm{
      static_cast<double>(v[1].x()) - p0x,
      static_cast<double>(v[2].x()) - p0x,
      static_cast<double>(v[1].y()) - p0y,
      static_cast<double>(v[2].y()) - p0y};
    const double det_p = det(pm);
    const Mat2 pm_inv{
      pm.a22 / det_p,
      -pm.a12 / det_p,
      -pm.a21 / det_p,
      pm.a11 / det_p};
    const Mat2 qm{
      q1.x() - q0.x(),
      q2.x() - q0.x(),
      q1.y() - q0.y(),
      q2.y() - q0.y()};
    const Mat2 jt = qm * pm_inv;
    if (jac != nullptr) {
      *jac = jt * (*jac);
    }
    const double dx = p.x - p0x, dy = p.y - p0y;
    return {
      q0.x() + jt.a11 * dx + jt.a12 * dy,
      q0.y() + jt.a21 * dx + jt.a22 * dy};
  };

  for (NodeSet *cloud : {&conf_.nodes, &conf_.eval_nodes}) {
    for (std::size_t j = 0; j < cloud->size(); ++j) {
      cloud->y[j] = map_point(cloud->y[j], &cloud->jac[j]);
    }
  }
  for (auto *pts :
       {&conf_.boundary,
        &conf_.corners,
        &conf_.fd_e1,
        &conf_.fd_e2,
        &conf_.fd_e1h,
        &conf_.fd_e2h}) {
    for (auto &p : *pts) {
      p = map_point(p, nullptr);
    }
  }
  for (auto &line : conf_.grid_lines) {
    for (auto &p : line) {
      p = map_point(p, nullptr);
    }
  }
  ++conf_.n_passes;
  timer.stop("Conformality tracers");
}

void InsetState::scale_conformality_tracers(const double factor)
{
  for (NodeSet *cloud : {&conf_.nodes, &conf_.eval_nodes}) {
    for (std::size_t j = 0; j < cloud->size(); ++j) {
      cloud->y[j] = {factor * cloud->y[j].x, factor * cloud->y[j].y};
      Mat2 &jac = cloud->jac[j];
      jac = {
        factor * jac.a11,
        factor * jac.a12,
        factor * jac.a21,
        factor * jac.a22};
    }
  }
  for (auto *pts :
       {&conf_.boundary,
        &conf_.corners,
        &conf_.fd_e1,
        &conf_.fd_e2,
        &conf_.fd_e1h,
        &conf_.fd_e2h}) {
    for (auto &p : *pts) {
      p = {factor * p.x, factor * p.y};
    }
  }
  for (auto &line : conf_.grid_lines) {
    for (auto &p : line) {
      p = {factor * p.x, factor * p.y};
    }
  }
  ++conf_.n_doublings;
}

void InsetState::finalize_conformality_t0()
{
  // The frozen Jacobian of the postprocessor is the realized det DT0
  for (NodeSet *cloud : {&conf_.nodes, &conf_.eval_nodes}) {
    for (std::size_t j = 0; j < cloud->size(); ++j) {
      cloud->g0[j] = det(cloud->jac[j]);
    }
  }
  conf_.grid_lines_t0 = conf_.grid_lines;
}

void InsetState::report_conformality_t0() const
{
  const double lx = static_cast<double>(lx_);
  const double ly = static_cast<double>(ly_);
  const double domain_area = lx * ly;

  const BoundaryCheck bc = check_boundary(conf_, lx, ly);
  const FdCheck fc = check_finite_differences(conf_);

  const CloudSummary opt =
    summarize_cloud(conf_.nodes, conf_.region, geo_divs_, domain_area);
  const CloudSummary eval = summarize_cloud(
    conf_.eval_nodes,
    conf_.eval_region,
    geo_divs_,
    domain_area);

  const bool polygons_exact = !conf_.simplified_after_seeding &&
                              !args_.disable_simplification_densification;

  std::cerr << std::setprecision(6) << "\nConformality tracers for inset "
            << pos_ << ": " << conf_.n_passes << " passes, "
            << conf_.n_doublings << " grid doublings, final lattice " << lx_
            << " x " << ly_ << std::endl;
  std::cerr << "  E0(T0) = " << opt.energy << " (evaluation cloud "
            << eval.energy << "), on land " << opt.energy_land << " ("
            << eval.energy_land << ")" << std::endl;
  std::cerr << "  min det J = " << opt.min_det
            << ", sum w det J / area - 1 = " << opt.total_det_rel_error << " ("
            << eval.total_det_rel_error << ")" << std::endl;
  std::cerr << "  region areas vs polygons: max rel. error "
            << opt.max_region_rel_error << " (" << eval.max_region_rel_error
            << "), area-weighted mean " << opt.mean_region_rel_error << " ("
            << eval.mean_region_rel_error << ")" << std::endl;
  std::cerr << "  boundary: max offset from side " << bc.max_side_offset
            << ", max overshoot along side " << bc.max_along_overshoot
            << ", max corner motion " << bc.max_corner_motion << std::endl;
  std::cerr << "  finite differences: " << fc.n_within << " pairs within "
            << fd_crossing_threshold << " (Richardson median "
            << fc.richardson_median << ", max " << fc.richardson_max << "), "
            << fc.n_crossing << " crossing an affine piece" << std::endl;
  if (!polygons_exact) {
    std::cerr << "  WARNING: polygons were simplified or not densified, so "
                 "they are not exact images under T0; region-area "
                 "comparisons include that discrepancy."
              << std::endl;
  }

  const nlohmann::json report = {
    {"inset", pos_},
    {"initial_lattice", {conf_.lx0, conf_.ly0}},
    {"final_lattice", {lx_, ly_}},
    {"n_passes", conf_.n_passes},
    {"n_doublings", conf_.n_doublings},
    {"polygons_exact_images", polygons_exact},
    {"optimization_cloud", cloud_to_json(opt, geo_divs_, conf_.quad)},
    {"evaluation_cloud", cloud_to_json(eval, geo_divs_, conf_.eval_quad)},
    {"boundary", boundary_to_json(bc, conf_.boundary.size())},
    {"finite_differences", fd_to_json(fc, conf_)},
    {"distortion",
     {{"rectangle",
       distortion_to_json(
         distortion_stats(conf_.eval_nodes, conf_.eval_region, false))},
      {"land",
       distortion_to_json(
         distortion_stats(conf_.eval_nodes, conf_.eval_region, true))}}}};
  std::ofstream out(inset_name_ + "_conformality_t0.json");
  out << std::setw(2) << report << std::endl;
}

namespace
{

// Cell areas of a regular mesh of the final lattice after replaying the
// correction A. Each cell boundary is sampled at s points per edge, and the
// image cell is the polygon through the replayed samples. Besides any loss
// of area preservation by A, the polygon misses the curvature of the
// images of the edges, an O((h / s)^2) error for cells of size h. Two
// refinements separate the two:
//   - finer meshes with s = 1: the error falls by about 4 per halving of h
//     once in the asymptotic regime;
//   - more samples per edge on a fixed mesh: the error tends to the area
//     error of A itself.
struct MeshAreaCheck {
  std::size_t cells_long{0};
  std::size_t samples_per_edge{1};
  std::size_t n_cells{0};
  double max_rel_error{0.0};
  double mean_rel_error{0.0};
};

// Replays A once on a lattice with base_cells * max_refinement intervals
// along the longer side (max_refinement a power of 2) and returns the
// checks for meshes of base_cells * r cells with s = 1, followed by the
// mesh of base_cells cells with s = r, for r = 1, 2, 4, ..., max_refinement
std::vector<MeshAreaCheck> check_mesh_areas(
  const conformality::SpectralBasis &basis,
  const std::vector<conformality::AcceptedStep> &steps,
  const double lx,
  const double ly,
  const std::size_t base_cells,
  const std::size_t max_refinement)
{
  const double longer = std::max(lx, ly);
  const auto ncx = std::max<std::size_t>(
    1,
    static_cast<std::size_t>(
      std::lround(static_cast<double>(base_cells) * lx / longer)));
  const auto ncy = std::max<std::size_t>(
    1,
    static_cast<std::size_t>(
      std::lround(static_cast<double>(base_cells) * ly / longer)));
  const std::size_t nx = ncx * max_refinement;
  const std::size_t ny = ncy * max_refinement;
  const double fx = lx / static_cast<double>(nx);
  const double fy = ly / static_cast<double>(ny);
  std::vector<Vec2> pts;
  pts.reserve((nx + 1) * (ny + 1));
  for (std::size_t i = 0; i <= nx; ++i) {
    for (std::size_t j = 0; j <= ny; ++j) {
      pts.push_back(
        {static_cast<double>(i) * fx, static_cast<double>(j) * fy});
    }
  }
  static_cast<void>(conformality::replay(pts, nullptr, basis, steps));
  auto at = [&](const std::size_t i, const std::size_t j) {
    return pts[i * (ny + 1) + j];
  };

  // Mesh of cells spanning `stride` lattice intervals, with boundary
  // samples every `step` intervals
  auto evaluate = [&](const std::size_t stride, const std::size_t step) {
    MeshAreaCheck m;
    m.cells_long = base_cells * (max_refinement / stride);
    m.samples_per_edge = stride / step;
    m.n_cells = (nx / stride) * (ny / stride);
    const double nominal =
      static_cast<double>(stride) * fx * static_cast<double>(stride) * fy;
    double sum = 0.0;
    for (std::size_t i0 = 0; i0 < nx; i0 += stride) {
      for (std::size_t j0 = 0; j0 < ny; j0 += stride) {

        // Shoelace formula relative to the first corner, counterclockwise
        const Vec2 origin = at(i0, j0);
        double twice_area = 0.0;
        Vec2 prev{0.0, 0.0};
        auto visit = [&](const std::size_t i, const std::size_t j) {
          const Vec2 q = at(i, j);
          const Vec2 p{q.x - origin.x, q.y - origin.y};
          twice_area += prev.x * p.y - p.x * prev.y;
          prev = p;
        };
        for (std::size_t k = step; k <= stride; k += step) {
          visit(i0 + k, j0);
        }
        for (std::size_t k = step; k <= stride; k += step) {
          visit(i0 + stride, j0 + k);
        }
        for (std::size_t k = step; k <= stride; k += step) {
          visit(i0 + stride - k, j0 + stride);
        }
        for (std::size_t k = step; k <= stride; k += step) {
          visit(i0, j0 + stride - k);
        }
        const double rel = std::abs(0.5 * twice_area / nominal - 1.0);
        m.max_rel_error = std::max(m.max_rel_error, rel);
        sum += rel;
      }
    }
    m.mean_rel_error = sum / static_cast<double>(m.n_cells);
    return m;
  };
  std::vector<MeshAreaCheck> checks;
  for (std::size_t r = 1; r <= max_refinement; r *= 2) {
    checks.push_back(evaluate(max_refinement / r, max_refinement / r));
  }
  for (std::size_t r = 2; r <= max_refinement; r *= 2) {
    checks.push_back(evaluate(max_refinement, max_refinement / r));
  }
  return checks;
}

// Replays the correction on several point sets at once, so that every
// substep is one parallel loop over all their points
conformality::FlowStats replay_sets(
  const std::vector<std::vector<Vec2> *> &sets,
  const conformality::SpectralBasis &basis,
  const std::vector<conformality::AcceptedStep> &steps)
{
  std::vector<Vec2> all;
  std::size_t total = 0;
  for (const auto *set : sets) {
    total += set->size();
  }
  all.reserve(total);
  for (const auto *set : sets) {
    all.insert(all.end(), set->begin(), set->end());
  }
  const auto stats = conformality::replay(all, nullptr, basis, steps);
  std::size_t offset = 0;
  for (auto *set : sets) {
    for (auto &p : *set) {
      p = all[offset++];
    }
  }
  return stats;
}

nlohmann::json mesh_to_json(const MeshAreaCheck &m)
{
  return {
    {"cells_long", m.cells_long},
    {"samples_per_edge", m.samples_per_edge},
    {"n_cells", m.n_cells},
    {"max_rel_area_error", m.max_rel_error},
    {"mean_rel_area_error", m.mean_rel_error}};
}

}  // namespace

void InsetState::run_conformality_optimization()
{
  if (!converge_ || intersections_found_) {
    std::cerr << "WARNING: skipping the conformality optimization because "
                 "the cartogram did not converge or has intersections."
              << std::endl;
    return;
  }
  const double lx = static_cast<double>(lx_);
  const double ly = static_cast<double>(ly_);
  const double domain_area = lx * ly;

  // Modes along each axis, proportional to the lattice side lengths
  const double longer = std::max(lx, ly);
  const auto mx = std::max<std::size_t>(
    1,
    static_cast<std::size_t>(std::lround(args_.conf_modes * lx / longer)));
  const auto my = std::max<std::size_t>(
    1,
    static_cast<std::size_t>(std::lround(args_.conf_modes * ly / longer)));
  const conformality::SpectralBasis basis(lx, ly, mx, my);

  conformality::OptimizerSettings settings;
  settings.max_iterations = args_.conf_max_iter;
  settings.on_iteration = [](const conformality::IterationRecord &r) {
    std::cerr << "  conformality iteration " << r.iteration << ": E "
              << r.energy << " -> " << r.energy_after << ", R " << r.r
              << ", h " << r.h << " (strain " << r.strain << "), "
              << r.n_substeps << " substeps, " << r.backtracks
              << " backtracks, " << r.seconds << " s" << std::endl;
  };

  std::cerr << "\nConformality optimization for inset " << pos_ << ": " << mx
            << " x " << my << " modes, " << conf_.nodes.size() << " nodes, "
            << conformality::max_threads() << " threads" << std::endl;
  const DistortionStats rect_before =
    distortion_stats(conf_.eval_nodes, conf_.eval_region, false);
  const DistortionStats land_before =
    distortion_stats(conf_.eval_nodes, conf_.eval_region, true);
  const CloudSummary opt_before =
    summarize_cloud(conf_.nodes, conf_.region, geo_divs_, domain_area);
  const CloudSummary eval_before = summarize_cloud(
    conf_.eval_nodes,
    conf_.eval_region,
    geo_divs_,
    domain_area);

  timer.start("Conformality optimization");
  const auto t0 = std::chrono::steady_clock::now();
  const conformality::OptimizationResult result =
    conformality::optimize(conf_.nodes, basis, settings);
  const auto t1 = std::chrono::steady_clock::now();

  // Apply the same correction to every other tracer set
  conformality::FlowStats replay_stats =
    conformality::replay(conf_.eval_nodes, basis, result.steps);

  // Output polygons of T0, densified so that the curved images of their
  // edges under A are resolved at a quarter lattice unit; they are
  // corrected below together with the other passive tracers
  conf_.t_rings.clear();
  conf_.t_ring_region.clear();
  for (std::size_t g = 0; g < geo_divs_.size(); ++g) {
    for (const auto &pwh : geo_divs_[g].polygons_with_holes()) {
      auto add_ring = [&](const Polygon &ring) {
        const std::vector<Point> pts(ring.begin(), ring.end());
        std::vector<Vec2> dense;
        for (std::size_t i = 0; i < pts.size(); ++i) {
          const Point &a = pts[i];
          const Point &b = pts[(i + 1) % pts.size()];
          const double len = std::hypot(b.x() - a.x(), b.y() - a.y());
          const auto n = std::max<std::size_t>(
            1,
            static_cast<std::size_t>(std::ceil(len / 0.25)));
          for (std::size_t k = 0; k < n; ++k) {
            const double t = static_cast<double>(k) / static_cast<double>(n);
            dense.push_back(
              {a.x() + t * (b.x() - a.x()), a.y() + t * (b.y() - a.y())});
          }
        }
        conf_.t_rings.push_back(std::move(dense));
        conf_.t_ring_region.push_back(static_cast<int>(g));
      };
      add_ring(pwh.outer_boundary());
      for (const auto &hole : pwh.holes()) {
        add_ring(hole);
      }
    }
  }
  std::vector<std::vector<Vec2> *> passive = {
    &conf_.boundary,
    &conf_.corners,
    &conf_.fd_e1,
    &conf_.fd_e2,
    &conf_.fd_e1h,
    &conf_.fd_e2h};
  for (auto &line : conf_.grid_lines) {
    passive.push_back(&line);
  }
  for (auto &ring : conf_.t_rings) {
    passive.push_back(&ring);
  }
  conformality::accumulate(
    replay_stats,
    replay_sets(passive, basis, result.steps));
  conf_.optimized = true;

  const std::vector<MeshAreaCheck> mesh_checks =
    check_mesh_areas(basis, result.steps, lx, ly, 64, 8);
  const auto t2 = std::chrono::steady_clock::now();
  timer.stop("Conformality optimization");

  const double optimize_seconds =
    std::chrono::duration<double>(t1 - t0).count();
  const double replay_seconds = std::chrono::duration<double>(t2 - t1).count();
  const DistortionStats rect_after =
    distortion_stats(conf_.eval_nodes, conf_.eval_region, false);
  const DistortionStats land_after =
    distortion_stats(conf_.eval_nodes, conf_.eval_region, true);
  const CloudSummary opt_after =
    summarize_cloud(conf_.nodes, conf_.region, geo_divs_, domain_area);
  const CloudSummary eval_after = summarize_cloud(
    conf_.eval_nodes,
    conf_.eval_region,
    geo_divs_,
    domain_area);
  const double drift_opt = conformality::max_det_drift(conf_.nodes);
  const double drift_eval = conformality::max_det_drift(conf_.eval_nodes);
  const BoundaryCheck bc = check_boundary(conf_, lx, ly);
  const FdCheck fc = check_finite_differences(conf_);

  std::cerr << std::setprecision(6) << "Conformality optimization: "
            << conformality::to_string(result.outcome) << " (" << result.reason
            << "), " << result.steps.size() << " accepted steps in "
            << optimize_seconds << " s (replay " << replay_seconds << " s)"
            << std::endl;
  std::cerr << "  E0: " << opt_before.energy << " -> " << opt_after.energy
            << " (evaluation cloud " << eval_before.energy << " -> "
            << eval_after.energy << "); on land " << eval_before.energy_land
            << " -> " << eval_after.energy_land << std::endl;
  std::cerr << "  mean ln(a/b): rectangle " << rect_before.e_mean << " -> "
            << rect_after.e_mean << ", land " << land_before.e_mean << " -> "
            << land_after.e_mean << std::endl;
  std::cerr << "  max |det J / g0 - 1|: " << drift_opt << " (evaluation cloud "
            << drift_eval << ")" << std::endl;
  std::cerr << "  check-mesh cell areas, max rel. error:";
  for (const auto &m : mesh_checks) {
    std::cerr << " " << m.max_rel_error << " (" << m.cells_long << " cells, "
              << m.samples_per_edge << "/edge)";
  }
  std::cerr << std::endl;
  std::cerr << "  boundary: max offset from side " << bc.max_side_offset
            << ", max corner motion " << bc.max_corner_motion
            << "; finite differences: " << fc.n_within << " pairs within "
            << fd_crossing_threshold << " (Richardson max "
            << fc.richardson_max << ", forward max " << fc.forward_max << "), "
            << fc.n_crossing << " crossing" << std::endl;
  std::cerr << "  NOTE: output polygons are those of T0; the correction is "
               "not yet applied to them."
            << std::endl;

  nlohmann::json mesh_json = nlohmann::json::array();
  for (const auto &m : mesh_checks) {
    mesh_json.push_back(mesh_to_json(m));
  }

  // Replay log: applying these steps in order with the same basis and
  // solver settings reproduces A on any point set
  nlohmann::json steps_json = nlohmann::json::array();
  for (const auto &step : result.steps) {
    steps_json.push_back(
      {{"h", step.h},
       {"n_substeps", step.n_substeps},
       {"coeffs", step.coeffs}});
  }
  const conformality::ImplicitMidpointSettings &solver = settings.solver;
  const nlohmann::json step_log = {
    {"inset", pos_},
    {"lattice", {lx_, ly_}},
    {"modes", {mx, my}},
    {"coefficient_index", "c[(m - 1) * my + (n - 1)] for mode (m, n)"},
    {"solver",
     {{"residual_tol", solver.residual_tol},
      {"max_iterations", solver.max_iterations}}},
    {"steps", steps_json}};
  std::ofstream(inset_name_ + "_conformality_steps.json") << step_log << '\n';

  nlohmann::json history = nlohmann::json::array();
  for (const auto &r : result.history) {
    history.push_back(
      {{"iteration", r.iteration},
       {"E", r.energy},
       {"R", r.r},
       {"slope", r.slope},
       {"h", r.h},
       {"strain", r.strain},
       {"substeps", r.n_substeps},
       {"backtracks", r.backtracks},
       {"E_after", r.energy_after},
       {"armijo_margin", r.armijo_margin},
       {"max_det_drift", r.max_det_drift},
       {"max_solver_iterations", r.max_solver_iterations},
       {"max_solver_residual", r.max_solver_residual},
       {"seconds", r.seconds}});
  }
  const nlohmann::json report = {
    {"inset", pos_},
    {"lattice", {lx_, ly_}},
    {"modes", {mx, my}},
    {"settings",
     {{"sigma", settings.sigma},
      {"max_iterations", settings.max_iterations},
      {"gradient_tol", settings.gradient_tol},
      {"energy_tol", settings.energy_tol},
      {"stall_window", settings.stall_window},
      {"stall_tol", settings.stall_tol},
      {"max_strain", settings.max_strain},
      {"substep_factor", settings.substep_factor},
      {"det_drift_tol", settings.det_drift_tol},
      {"openmp", conformality::openmp_enabled()},
      {"threads", conformality::max_threads()}}},
    {"outcome", conformality::to_string(result.outcome)},
    {"reason", result.reason},
    {"n_accepted_steps", result.steps.size()},
    {"sqrt_R_over_R0",
     result.initial_r > 0.0 ? std::sqrt(result.final_r / result.initial_r)
                            : 0.0},
    {"seconds", {{"optimize", optimize_seconds}, {"replay", replay_seconds}}},
    {"before",
     {{"optimization_cloud", cloud_to_json(opt_before, geo_divs_, conf_.quad)},
      {"evaluation_cloud",
       cloud_to_json(eval_before, geo_divs_, conf_.eval_quad)},
      {"distortion_rectangle", distortion_to_json(rect_before)},
      {"distortion_land", distortion_to_json(land_before)}}},
    {"after",
     {{"optimization_cloud", cloud_to_json(opt_after, geo_divs_, conf_.quad)},
      {"evaluation_cloud",
       cloud_to_json(eval_after, geo_divs_, conf_.eval_quad)},
      {"distortion_rectangle", distortion_to_json(rect_after)},
      {"distortion_land", distortion_to_json(land_after)}}},
    {"checks",
     {{"max_det_drift_optimization_cloud", drift_opt},
      {"max_det_drift_evaluation_cloud", drift_eval},
      {"replay_converged", replay_stats.converged},
      {"mesh", mesh_json},
      {"boundary", boundary_to_json(bc, conf_.boundary.size())},
      {"finite_differences", fd_to_json(fc, conf_)}}},
    {"history", history}};
  std::ofstream out(inset_name_ + "_conformality_opt.json");
  out << std::setw(2) << report << std::endl;
}

namespace
{

// Coordinates rounded to 1e-4 lattice units keep the plot data compact
nlohmann::json points_to_json(const std::vector<Vec2> &pts)
{
  nlohmann::json out = nlohmann::json::array();
  for (const Vec2 p : pts) {
    out.push_back({std::round(p.x * 1e4) / 1e4, std::round(p.y * 1e4) / 1e4});
  }
  return out;
}

nlohmann::json rings_to_json(
  const std::vector<std::vector<Vec2>> &rings,
  const std::vector<int> &region)
{
  nlohmann::json out = nlohmann::json::array();
  for (std::size_t i = 0; i < rings.size(); ++i) {
    out.push_back({{"region", region[i]}, {"ring", points_to_json(rings[i])}});
  }
  return out;
}

nlohmann::json lines_to_json(const std::vector<std::vector<Vec2>> &lines)
{
  nlohmann::json out = nlohmann::json::array();
  for (const auto &line : lines) {
    out.push_back(points_to_json(line));
  }
  return out;
}

}  // namespace

void InsetState::write_conformality_plot_data() const
{
  // Output polygons of T0 (exact images of the input polygons when
  // simplification is off)
  std::vector<std::vector<Vec2>> t0_rings;
  std::vector<int> t0_region;
  for (std::size_t g = 0; g < geo_divs_.size(); ++g) {
    for (const auto &pwh : geo_divs_[g].polygons_with_holes()) {
      auto add_ring = [&](const Polygon &ring) {
        std::vector<Vec2> pts;
        for (const auto &p : ring) {
          pts.push_back({p.x(), p.y()});
        }
        t0_rings.push_back(std::move(pts));
        t0_region.push_back(static_cast<int>(g));
      };
      add_ring(pwh.outer_boundary());
      for (const auto &hole : pwh.holes()) {
        add_ring(hole);
      }
    }
  }
  nlohmann::json names = nlohmann::json::array();
  for (const auto &gd : geo_divs_) {
    names.push_back(gd.id());
  }
  const nlohmann::json data = {
    {"inset", pos_},
    {"initial_lattice", {conf_.lx0, conf_.ly0}},
    {"final_lattice", {lx_, ly_}},
    {"grid_spacing", conf_.grid_spacing},
    {"region_names", names},
    {"grid",
     {{"t0", lines_to_json(conf_.grid_lines_t0)},
      {"t",
       conf_.optimized ? lines_to_json(conf_.grid_lines) : nlohmann::json()}}},
    {"polygons",
     {{"original", rings_to_json(conf_.input_rings, conf_.input_ring_region)},
      {"t0", rings_to_json(t0_rings, t0_region)},
      {"t",
       conf_.optimized ? rings_to_json(conf_.t_rings, conf_.t_ring_region)
                       : nlohmann::json()}}}};
  std::ofstream out(inset_name_ + "_conformality_plot.json");
  out << data << std::endl;
}
