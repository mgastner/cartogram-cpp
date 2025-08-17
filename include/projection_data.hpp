#pragma once

#include "quadtree_corner.hpp"
#include <algorithm>
#include <cstdint>
#include <memory>
#include <numeric>
#include <vector>

class ProjectionData
{
private:
  // Current grid size
  uint32_t lx_{0}, ly_{0};

  // Flat buffer of size lx_ * ly_, uninitialized PODs
  // One 64-bit slot per cell: [ high 32 = generation | low 32 = index ]
  // TODO: consider if top 8 bits for generation and rest 24 bits for index is
  // enough
  std::unique_ptr<uint64_t[]> corner_slot_{};

  // Bumps each time you call build_fast_indexing()
  uint32_t gen_{1};

  std::vector<Point> projection_;

  uint32_t get_buf_index(uint32_t x, uint32_t y) const noexcept
  {
    return x * ly_ + y;
  }

  size_t buf_size() const noexcept
  {
    return static_cast<size_t>(lx_) * static_cast<size_t>(ly_);
  }

  bool in_bounds(uint32_t x, uint32_t y) const noexcept
  {
    return x < lx_ && y < ly_;
  }

public:
  ProjectionData() = default;

  ProjectionData(const ProjectionData &o)
      : lx_(o.lx_), ly_(o.ly_),
        corner_slot_(
          o.corner_slot_
            ? std::unique_ptr<uint64_t[]>{new uint64_t[buf_size()]}
            : nullptr),
        gen_(o.gen_), projection_(o.projection_)
  {
    if (corner_slot_)
      std::copy_n(o.corner_slot_.get(), buf_size(), corner_slot_.get());
  }

  ProjectionData &operator=(const ProjectionData &o)
  {
    if (this == &o)
      return *this;
    lx_ = o.lx_;
    ly_ = o.ly_;
    gen_ = o.gen_;
    if (o.corner_slot_) {
      corner_slot_.reset(new uint64_t[o.buf_size()]);
      std::copy_n(o.corner_slot_.get(), o.buf_size(), corner_slot_.get());
    } else {
      corner_slot_.reset();
    }
    projection_ = o.projection_;
    return *this;
  }

  ProjectionData(ProjectionData &&) noexcept = default;
  ProjectionData &operator=(ProjectionData &&) noexcept = default;
  ~ProjectionData() = default;

  void reserve(uint32_t new_lx, uint32_t new_ly)
  {
    const uint64_t new_area = uint64_t(new_lx) * uint64_t(new_ly);
    const uint64_t cap_area = uint64_t(lx_) * uint64_t(ly_);

    // Update shape even if capacity is already enough
    lx_ = new_lx;
    ly_ = new_ly;

    // If we already have enough capacity and a buffer, keep it
    if (corner_slot_ && new_area <= cap_area)
      return;

    assert(new_area <= std::numeric_limits<uint32_t>::max());

    corner_slot_.reset(new uint64_t[buf_size()]());
  }

  [[nodiscard]] std::vector<Point> &get_projection() noexcept
  {
    return projection_;
  }

  // We use `gen` because we would like to keep the same underlying buffer
  // active across multiple intgrations and also be able to check if a corner
  // is valid for the current integration. So, we use a generation
  // number that increments each time we build the fast indexing.
  // This allows us to check if a corner is valid for the current generation
  // quickly by just checking the top 32 bits of the index
  void build_fast_indexing(const std::vector<QuadtreeCorner> &keys) noexcept
  {
    ++gen_;
    assert(corner_slot_ && "reserve() must be called before indexing");
    for (uint32_t i = 0; i < keys.size(); ++i) {
      const uint32_t x = keys[i].x();
      const uint32_t y = keys[i].y();

      assert(in_bounds(x, y) && "key outside reserved grid");

      const uint32_t idx = get_buf_index(x, y);
      corner_slot_[idx] = (uint64_t(gen_) << 32) | uint64_t(i);
    }
  }

  // Check if the corner (x, y) is valid in the current generation
  [[nodiscard]] bool is_valid_corner(uint32_t x, uint32_t y) const noexcept
  {
    assert(corner_slot_ && "reserve() must be called before indexing");
    if (!in_bounds(x, y))
      return false;
    const uint64_t idx = corner_slot_[get_buf_index(x, y)];
    return uint32_t(idx >> 32) == gen_;
  }

  /// Precondition: (x, y) is valid for the current generation
  [[nodiscard]] uint32_t offset(uint32_t x, uint32_t y) const noexcept
  {
    assert(in_bounds(x, y) && "offset() out of bounds");
    assert(is_valid_corner(x, y) && "invalid corner");
    return static_cast<uint32_t>(corner_slot_[get_buf_index(x, y)]);
  }

  [[nodiscard]] Point get(uint32_t x, uint32_t y) const noexcept
  {
    const uint32_t idx = offset(x, y);
    assert(idx < projection_.size() && "projection index out of range");
    return projection_[idx];
  }

  [[nodiscard]] uint32_t num_unique_corners() const noexcept
  {
    return static_cast<uint32_t>(projection_.size());
  }
};
