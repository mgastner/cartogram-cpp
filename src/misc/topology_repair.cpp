#include "topology_repair.hpp"
#include <cmath>
#include <geos_c.h>
#include <stdexcept>
#include <utility>

namespace
{

// GEOS context with its error handler; GEOS errors become exceptions
class GeosContext
{
public:
  GeosContext() : handle_(GEOS_init_r())
  {
    GEOSContext_setErrorMessageHandler_r(handle_, &store_message, &error_);
  }
  ~GeosContext()
  {
    GEOS_finish_r(handle_);
  }
  GeosContext(const GeosContext &) = delete;
  GeosContext &operator=(const GeosContext &) = delete;

  [[nodiscard]] GEOSContextHandle_t get() const
  {
    return handle_;
  }

  // Throws if p is null, with the last GEOS error message
  template <typename T> T *check(T *p, const char *what) const
  {
    if (p == nullptr) {
      throw std::runtime_error(std::string(what) + " failed: " + error_);
    }
    return p;
  }

private:
  static void store_message(const char *message, void *userdata)
  {
    *static_cast<std::string *>(userdata) = message;
  }

  GEOSContextHandle_t handle_;
  std::string error_;
};

// Owning pointer to a GEOS geometry
class Geom
{
public:
  Geom(const GeosContext &ctx, GEOSGeometry *g) : ctx_(&ctx), g_(g) {}
  ~Geom()
  {
    if (g_ != nullptr) {
      GEOSGeom_destroy_r(ctx_->get(), g_);
    }
  }
  Geom(Geom &&other) noexcept
      : ctx_(other.ctx_), g_(std::exchange(other.g_, nullptr))
  {
  }
  Geom &operator=(Geom &&other) noexcept
  {
    std::swap(ctx_, other.ctx_);
    std::swap(g_, other.g_);
    return *this;
  }
  Geom(const Geom &) = delete;
  Geom &operator=(const Geom &) = delete;

  [[nodiscard]] GEOSGeometry *get() const
  {
    return g_;
  }

  // Gives up ownership, e.g. to a collection that takes it
  GEOSGeometry *release()
  {
    return std::exchange(g_, nullptr);
  }

private:
  const GeosContext *ctx_;
  GEOSGeometry *g_;
};

GEOSGeometry *make_ring(const GeosContext &ctx, const Polygon &ring)
{
  std::vector<double> xy;
  xy.reserve(2 * (ring.size() + 1));
  for (const auto &p : ring) {
    xy.push_back(p.x());
    xy.push_back(p.y());
  }
  xy.push_back(ring.vertex(0).x());
  xy.push_back(ring.vertex(0).y());
  GEOSCoordSequence *seq = ctx.check(
    GEOSCoordSeq_copyFromBuffer_r(
      ctx.get(),
      xy.data(),
      static_cast<unsigned int>(ring.size() + 1),
      0,
      0),
    "GEOSCoordSeq_copyFromBuffer");
  return ctx.check(
    GEOSGeom_createLinearRing_r(ctx.get(), seq),
    "GEOSGeom_createLinearRing");
}

Geom to_geos(const GeosContext &ctx, const GeoDiv &gd)
{
  std::vector<GEOSGeometry *> polygons;
  for (const auto &pwh : gd.polygons_with_holes()) {
    GEOSGeometry *shell = make_ring(ctx, pwh.outer_boundary());
    std::vector<GEOSGeometry *> holes;
    for (const auto &h : pwh.holes()) {
      holes.push_back(make_ring(ctx, h));
    }
    polygons.push_back(ctx.check(
      GEOSGeom_createPolygon_r(
        ctx.get(),
        shell,
        holes.data(),
        static_cast<unsigned int>(holes.size())),
      "GEOSGeom_createPolygon"));
  }
  return {
    ctx,
    ctx.check(
      GEOSGeom_createCollection_r(
        ctx.get(),
        GEOS_MULTIPOLYGON,
        polygons.data(),
        static_cast<unsigned int>(polygons.size())),
      "GEOSGeom_createCollection")};
}

Polygon ring_from_geos(const GeosContext &ctx, const GEOSGeometry *ring)
{
  const GEOSCoordSequence *seq =
    ctx.check(GEOSGeom_getCoordSeq_r(ctx.get(), ring), "GEOSGeom_getCoordSeq");
  unsigned int n = 0;
  GEOSCoordSeq_getSize_r(ctx.get(), seq, &n);
  Polygon out;

  // The last coordinate repeats the first one
  for (unsigned int i = 0; i + 1 < n; ++i) {
    double x = 0.0;
    double y = 0.0;
    GEOSCoordSeq_getXY_r(ctx.get(), seq, i, &x, &y);
    out.push_back(Point(x, y));
  }
  return out;
}

void append_polygon(
  const GeosContext &ctx,
  const GEOSGeometry *poly,
  GeoDiv &gd)
{
  Polygon ext = ring_from_geos(
    ctx,
    ctx.check(GEOSGetExteriorRing_r(ctx.get(), poly), "GEOSGetExteriorRing"));
  if (ext.is_clockwise_oriented()) {
    ext.reverse_orientation();
  }
  std::vector<Polygon> holes;
  const int n_holes = GEOSGetNumInteriorRings_r(ctx.get(), poly);
  for (int i = 0; i < n_holes; ++i) {
    Polygon h = ring_from_geos(
      ctx,
      ctx.check(
        GEOSGetInteriorRingN_r(ctx.get(), poly, i),
        "GEOSGetInteriorRingN"));
    if (h.is_counterclockwise_oriented()) {
      h.reverse_orientation();
    }
    holes.push_back(std::move(h));
  }
  gd.push_back(Polygon_with_holes(ext, holes.begin(), holes.end()));
}

// GeoDiv with the given ID and the polygons of a polygonal GEOS geometry
GeoDiv from_geos(const GeosContext &ctx, const GEOSGeometry *g, std::string id)
{
  GeoDiv gd(std::move(id));
  const int type = GEOSGeomTypeId_r(ctx.get(), g);
  if (type == GEOS_POLYGON) {
    append_polygon(ctx, g, gd);
  } else if (type == GEOS_MULTIPOLYGON) {
    const int n = GEOSGetNumGeometries_r(ctx.get(), g);
    for (int i = 0; i < n; ++i) {
      append_polygon(ctx, GEOSGetGeometryN_r(ctx.get(), g, i), gd);
    }
  } else {
    throw std::runtime_error(
      "repair produced a geometry that is not polygonal");
  }
  return gd;
}

// Valid geometry with the polygonal part of what MakeValid returns
Geom make_valid(const GeosContext &ctx, const GEOSGeometry *g)
{
  GEOSMakeValidParams *params = ctx.check(
    GEOSMakeValidParams_create_r(ctx.get()),
    "GEOSMakeValidParams_create");
  GEOSMakeValidParams_setMethod_r(
    ctx.get(),
    params,
    GEOS_MAKE_VALID_STRUCTURE);
  GEOSMakeValidParams_setKeepCollapsed_r(ctx.get(), params, 0);
  GEOSGeometry *valid = GEOSMakeValidWithParams_r(ctx.get(), g, params);
  GEOSMakeValidParams_destroy_r(ctx.get(), params);
  const Geom result(ctx, ctx.check(valid, "GEOSMakeValidWithParams"));

  // Collect the polygons (MakeValid may add collapsed lines or points)
  std::vector<GEOSGeometry *> polygons;
  auto collect = [&](const auto &self, const GEOSGeometry *part) -> void {
    const int type = GEOSGeomTypeId_r(ctx.get(), part);
    if (type == GEOS_POLYGON) {
      polygons.push_back(
        ctx.check(GEOSGeom_clone_r(ctx.get(), part), "GEOSGeom_clone"));
    } else if (type == GEOS_MULTIPOLYGON || type == GEOS_GEOMETRYCOLLECTION) {
      const int n = GEOSGetNumGeometries_r(ctx.get(), part);
      for (int i = 0; i < n; ++i) {
        self(self, GEOSGetGeometryN_r(ctx.get(), part, i));
      }
    }
  };
  collect(collect, result.get());
  return {
    ctx,
    ctx.check(
      GEOSGeom_createCollection_r(
        ctx.get(),
        GEOS_MULTIPOLYGON,
        polygons.data(),
        static_cast<unsigned int>(polygons.size())),
      "GEOSGeom_createCollection")};
}

// Area enclosed by the exterior rings minus the area of the holes, from the
// shoelace formula, which is defined also for rings that are not simple
double shoelace_area(const GeoDiv &gd)
{
  double a = 0.0;
  for (const auto &pwh : gd.polygons_with_holes()) {
    a += std::abs(pwh.outer_boundary().area());
    for (const auto &h : pwh.holes()) {
      a -= std::abs(h.area());
    }
  }
  return a;
}

// Normalized copy of a polygonal geometry as a MultiPolygon. The cleaner
// returns a Polygon where the input was a MultiPolygon with one part.
Geom normalized_multipolygon(const GeosContext &ctx, const GEOSGeometry *g)
{
  GEOSGeometry *copy =
    ctx.check(GEOSGeom_clone_r(ctx.get(), g), "GEOSGeom_clone");
  if (GEOSGeomTypeId_r(ctx.get(), copy) == GEOS_POLYGON) {
    copy = ctx.check(
      GEOSGeom_createCollection_r(ctx.get(), GEOS_MULTIPOLYGON, &copy, 1),
      "GEOSGeom_createCollection");
  }
  Geom result(ctx, copy);
  GEOSNormalize_r(ctx.get(), result.get());
  return result;
}

// True if both geometries have the same vertices in the same rings,
// regardless of the starting vertex and the order of rings and polygons
bool same_geometry(
  const GeosContext &ctx,
  const GEOSGeometry *a,
  const GEOSGeometry *b)
{
  const Geom na = normalized_multipolygon(ctx, a);
  const Geom nb = normalized_multipolygon(ctx, b);
  return GEOSEqualsExact_r(ctx.get(), na.get(), nb.get(), 0.0) == 1;
}

}  // namespace

bool all_rings_simple(const GeoDiv &gd)
{
  for (const auto &pwh : gd.polygons_with_holes()) {
    if (!pwh.outer_boundary().is_simple()) {
      return false;
    }
    for (const auto &h : pwh.holes()) {
      if (!h.is_simple()) {
        return false;
      }
    }
  }
  return true;
}

TopologyRepairResult repair_topology(
  std::vector<GeoDiv> &geo_divs,
  const TopologyRepairSettings &settings)
{
  const GeosContext ctx;
  const size_t n = geo_divs.size();

  // Step 1: rebuild regions with rings that are not simple
  std::vector<Geom> input;
  std::vector<bool> self_intersecting(n, false);
  input.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    Geom g = to_geos(ctx, geo_divs[i]);
    if (!all_rings_simple(geo_divs[i])) {
      self_intersecting[i] = true;
      g = make_valid(ctx, g.get());
    }
    input.push_back(std::move(g));
  }

  // Step 2: clean the coverage. The collection takes ownership of clones, so
  // that the inputs remain available for comparison.
  std::vector<GEOSGeometry *> clones;
  clones.reserve(n);
  for (const auto &g : input) {
    clones.push_back(
      ctx.check(GEOSGeom_clone_r(ctx.get(), g.get()), "GEOSGeom_clone"));
  }
  const Geom collection(
    ctx,
    ctx.check(
      GEOSGeom_createCollection_r(
        ctx.get(),
        GEOS_GEOMETRYCOLLECTION,
        clones.data(),
        static_cast<unsigned int>(n)),
      "GEOSGeom_createCollection"));
  GEOSCoverageCleanParams *params = ctx.check(
    GEOSCoverageCleanParams_create_r(ctx.get()),
    "GEOSCoverageCleanParams_create");
  if (settings.snapping_distance >= 0.0) {
    GEOSCoverageCleanParams_setSnappingDistance_r(
      ctx.get(),
      params,
      settings.snapping_distance);
  }
  GEOSGeometry *cleaned_raw =
    GEOSCoverageCleanWithParams_r(ctx.get(), collection.get(), params);
  GEOSCoverageCleanParams_destroy_r(ctx.get(), params);
  const Geom cleaned(
    ctx,
    ctx.check(cleaned_raw, "GEOSCoverageCleanWithParams"));
  if (
    GEOSGetNumGeometries_r(ctx.get(), cleaned.get()) != static_cast<int>(n)) {
    throw std::runtime_error(
      "GEOSCoverageCleanWithParams returned a different number of regions");
  }

  // Step 3: compare and convert back
  TopologyRepairResult result{{}, true};
  std::vector<std::pair<size_t, GeoDiv>> replacements;
  for (size_t i = 0; i < n; ++i) {
    const GEOSGeometry *c =
      GEOSGetGeometryN_r(ctx.get(), cleaned.get(), static_cast<int>(i));
    if (!self_intersecting[i] && same_geometry(ctx, c, input[i].get())) {
      continue;
    }
    const double before = shoelace_area(geo_divs[i]);
    double after = 0.0;
    if (GEOSisEmpty_r(ctx.get(), c) != 1) {
      GEOSArea_r(ctx.get(), c, &after);
    }
    const double change =
      before > 0.0 ? std::abs(after - before) / before : 1.0;
    result.changed.push_back({geo_divs[i].id(), self_intersecting[i], change});
    if (!(change <= settings.max_relative_area_change)) {
      result.within_tolerance = false;
    }
    if (GEOSisEmpty_r(ctx.get(), c) == 1) {
      result.within_tolerance = false;
      continue;
    }
    replacements.emplace_back(i, from_geos(ctx, c, geo_divs[i].id()));
  }
  if (result.within_tolerance) {
    for (auto &[i, gd] : replacements) {
      geo_divs[i] = std::move(gd);
    }
  }
  return result;
}
