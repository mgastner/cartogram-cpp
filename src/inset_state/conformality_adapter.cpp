// Adapter between the cartogram pipeline and the conformality core: seeds
// tracers of the cumulative map T0, moves them through every successful
// flattening pass with the same piecewise-affine map that projects the
// polygons, and reports diagnostics of the recorded map.

#include "conformality/objective.hpp"
#include "inset_state.hpp"
#include <algorithm>
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
       {&conf_.boundary, &conf_.corners, &conf_.fd_e1, &conf_.fd_e2}) {
    for (auto &p : *pts) {
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
       {&conf_.boundary, &conf_.corners, &conf_.fd_e1, &conf_.fd_e2}) {
    for (auto &p : *pts) {
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
}

void InsetState::report_conformality_t0() const
{
  const double lx = static_cast<double>(lx_);
  const double ly = static_cast<double>(ly_);
  const double domain_area = lx * ly;

  // Distance of boundary tracers from their side, and how far they moved
  // outside the side's extent
  double max_side_offset = 0.0;
  double max_along_overshoot = 0.0;
  for (std::size_t k = 0; k < conf_.boundary.size(); ++k) {
    const Vec2 p = conf_.boundary[k];
    double offset = 0.0;
    double along = 0.0;
    double extent = 0.0;
    switch (conf_.boundary_side[k]) {
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
    max_side_offset = std::max(max_side_offset, offset);
    max_along_overshoot =
      std::max({max_along_overshoot, -along, along - extent});
  }
  const std::vector<Vec2> corner_targets =
    {{0.0, 0.0}, {lx, 0.0}, {0.0, ly}, {lx, ly}};
  double max_corner_motion = 0.0;
  for (std::size_t k = 0; k < conf_.corners.size(); ++k) {
    max_corner_motion = std::max(
      max_corner_motion,
      std::hypot(
        conf_.corners[k].x - corner_targets[k].x,
        conf_.corners[k].y - corner_targets[k].y));
  }

  // Composed Jacobians against forward differences of companion tracers.
  // A pair whose companions fell into a different affine piece in some pass
  // shows an O(1) discrepancy; it is counted separately.
  constexpr double fd_crossing_threshold = 1e-6;
  std::vector<double> fd_errors;
  std::size_t n_crossing = 0;
  for (std::size_t k = 0; k < conf_.fd_parent.size(); ++k) {
    const std::size_t j = conf_.fd_parent[k];
    const Vec2 y = conf_.nodes.y[j];
    const double d = conf_.fd_delta;
    const Mat2 fd{
      (conf_.fd_e1[k].x - y.x) / d,
      (conf_.fd_e2[k].x - y.x) / d,
      (conf_.fd_e1[k].y - y.y) / d,
      (conf_.fd_e2[k].y - y.y) / d};
    const Mat2 &jac = conf_.nodes.jac[j];
    const Mat2 diff{
      fd.a11 - jac.a11,
      fd.a12 - jac.a12,
      fd.a21 - jac.a21,
      fd.a22 - jac.a22};
    const double rel = std::sqrt(frobenius_sq(diff) / frobenius_sq(jac));
    if (rel > fd_crossing_threshold) {
      ++n_crossing;
    } else {
      fd_errors.push_back(rel);
    }
  }
  std::sort(fd_errors.begin(), fd_errors.end());
  const double fd_median =
    fd_errors.empty() ? 0.0 : fd_errors[fd_errors.size() / 2];
  const double fd_max = fd_errors.empty() ? 0.0 : fd_errors.back();

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
  std::cerr << "  boundary: max offset from side " << max_side_offset
            << ", max overshoot along side " << max_along_overshoot
            << ", max corner motion " << max_corner_motion << std::endl;
  std::cerr << "  finite differences: " << fd_errors.size() << " pairs within "
            << fd_crossing_threshold << " (median " << fd_median << ", max "
            << fd_max << "), " << n_crossing << " crossing an affine piece"
            << std::endl;
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
    {"boundary",
     {{"n_tracers", conf_.boundary.size()},
      {"max_offset_from_side", max_side_offset},
      {"max_overshoot_along_side", max_along_overshoot},
      {"max_corner_motion", max_corner_motion}}},
    {"finite_differences",
     {{"delta", conf_.fd_delta},
      {"n_pairs", conf_.fd_parent.size()},
      {"n_within_threshold", fd_errors.size()},
      {"threshold", fd_crossing_threshold},
      {"median_rel_error", fd_median},
      {"max_rel_error", fd_max},
      {"n_crossing", n_crossing}}}};
  std::ofstream out(inset_name_ + "_conformality_t0.json");
  out << std::setw(2) << report << std::endl;
}
