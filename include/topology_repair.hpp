#ifndef TOPOLOGY_REPAIR_HPP_
#define TOPOLOGY_REPAIR_HPP_

#include "geo_div.hpp"
#include <string>
#include <vector>

// Repair of the input topology with GEOS. Two kinds of defects are repaired:
//
// 1. Rings that are not simple (self-intersections). The polygon is rebuilt
//    with GEOS MakeValid ("structure" method), which keeps the area enclosed
//    by the exterior rings and removes the area enclosed by the holes.
// 2. Neighbouring regions that do not share their borders exactly. GEOS
//    CoverageCleaner snaps vertices that lie on or within the snapping
//    distance of a neighbour's border into that border, so both regions get
//    the same vertices there, and it merges any overlap into one of the
//    overlapping regions (the one with the longest shared border).
//
// Without the second repair, CGAL's polyline simplification may simplify the
// two copies of a shared border differently, which opens gaps and overlaps
// between neighbours.

struct TopologyRepairSettings {

  // A region's area may change by at most this fraction of its original area
  double max_relative_area_change;

  // Vertices closer than this distance (in input coordinates) are snapped
  // together. A negative value lets GEOS choose: the diameter of the
  // bounding box of all regions divided by 1e8.
  double snapping_distance;
};

struct RegionRepair {
  std::string id;

  // True if one of the region's rings was not simple
  bool had_self_intersection;

  // |area after - area before| / area before. The area before is computed
  // with the shoelace formula, also for rings that are not simple.
  double relative_area_change;
};

struct TopologyRepairResult {

  // Regions whose geometry the repair changes, in input order
  std::vector<RegionRepair> changed;

  // False if a region's area would change by more than
  // max_relative_area_change; the GeoDivs are then left unmodified
  bool within_tolerance;
};

// Repairs the GeoDivs in place. Regions whose geometry the repair leaves
// unchanged keep their original vertices in their original order.
[[nodiscard]] TopologyRepairResult repair_topology(
  std::vector<GeoDiv> &geo_divs,
  const TopologyRepairSettings &settings);

// True if every exterior and interior ring of the GeoDiv is simple
[[nodiscard]] bool all_rings_simple(const GeoDiv &gd);

#endif  // TOPOLOGY_REPAIR_HPP_
