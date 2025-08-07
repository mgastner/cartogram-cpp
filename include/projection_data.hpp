#pragma once

#include "cgal_typedef.hpp"
#include <cstdint>
#include <memory>
#include <numeric>
#include <vector>

struct QuadtreeCorner {
public:
  QuadtreeCorner(uint32_t x, uint32_t y) : x_{x}, y_{y} {}

  uint32_t x() const noexcept
  {
    return x_;
  }
  uint32_t y() const noexcept
  {
    return y_;
  }

  operator Point() const noexcept
  {
    return Point(static_cast<double>(x_), static_cast<double>(y_));
  }

  auto operator<=>(const QuadtreeCorner &) const noexcept = default;
  bool operator==(const QuadtreeCorner &) const noexcept = default;

private:
  uint32_t x_;
  uint32_t y_;
};

class ProjectionData
{
private:
  Delaunay dt_;

  // Current grid size
  uint32_t lx_{0}, ly_{0};

  // Flat buffer of size lx_ * ly_, uninitialized PODs
  std::unique_ptr<uint32_t[]> corner_to_idx_{};
  std::vector<Point> projection_;

  uint32_t get_buf_index(uint32_t x, uint32_t y) const noexcept
  {
    return x * ly_ + y;
  }

  size_t buf_size() const noexcept
  {
    return static_cast<size_t>(lx_) * ly_;
  }

public:
  ProjectionData() = default;
  ProjectionData(const ProjectionData &o)
      : dt_(o.dt_), lx_(o.lx_), ly_(o.ly_),
        corner_to_idx_(o.corner_to_idx_ ? new uint32_t[buf_size()] : nullptr),
        projection_(o.projection_)
  {
    if (corner_to_idx_) {
      std::copy_n(o.corner_to_idx_.get(), buf_size(), corner_to_idx_.get());
    }
  }

  ProjectionData &operator=(const ProjectionData &o)
  {
    if (this == &o)
      return *this;
    dt_ = o.dt_;
    lx_ = o.lx_;
    ly_ = o.ly_;
    if (o.corner_to_idx_) {
      corner_to_idx_.reset(new uint32_t[o.buf_size()]);
      std::copy_n(o.corner_to_idx_.get(), o.buf_size(), corner_to_idx_.get());
    } else {
      corner_to_idx_.reset();
    }
    projection_ = o.projection_;
    return *this;
  }

  ProjectionData(ProjectionData &&) noexcept = default;
  ProjectionData &operator=(ProjectionData &&) noexcept = default;
  ~ProjectionData() = default;

  void reserve(uint32_t new_lx, uint32_t new_ly)
  {
    if (new_lx * new_ly <= lx_ * ly_)
      return;

    lx_ = new_lx;
    ly_ = new_ly;

    assert(
      static_cast<uint64_t>(lx_) * ly_ <=
      std::numeric_limits<uint32_t>::max());

    corner_to_idx_.reset(new uint32_t[lx_ * ly_]);
  }

  std::vector<Point> &get_projection() noexcept
  {
    return projection_;
  }

  void build_fast_indexing(const std::vector<QuadtreeCorner> &keys) noexcept
  {
    for (uint32_t i = 0; i < keys.size(); ++i) {
      uint32_t x = keys[i].x();
      uint32_t y = keys[i].y();
      corner_to_idx_[get_buf_index(x, y)] = i;
    }
  }

  Point get(uint32_t x, uint32_t y) const noexcept
  {
    uint32_t idx = corner_to_idx_[get_buf_index(x, y)];
    return projection_[idx];
  }

  const Delaunay &get_dt() const noexcept
  {
    return dt_;
  }

  Delaunay &get_dt() noexcept
  {
    return dt_;
  }

  void set_dt(Delaunay &&new_dt) noexcept
  {
    dt_ = std::move(new_dt);
  }
};
