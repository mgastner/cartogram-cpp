#pragma once
#include <cmath>
#include <cstdint>
#include <vector>

class QuadtreeLeafLocator
{
public:
  struct Leaf {
    uint32_t x, y, size;
    Leaf() = default;
    Leaf(uint32_t x_, uint32_t y_, uint32_t s_) noexcept
        : x(x_), y(y_), size(s_)
    {
    }
    bool operator==(const Leaf &) const noexcept = default;
  };

  struct Node {
    uint32_t x, y, size;
    int32_t first_child;
    Node() = default;
    Node(uint32_t x_, uint32_t y_, uint32_t s_, int32_t first_child_)
        : x(x_), y(y_), size(s_), first_child(first_child_)
    {
    }
    [[nodiscard]] bool is_leaf() const noexcept
    {
      return first_child < 0;
    }
  };

  QuadtreeLeafLocator() = default;

  template <typename QuadtreeNodeCollection>
  void build(uint32_t lx, uint32_t ly, const QuadtreeNodeCollection &nodes)
  {
    root_size_ = std::max(lx, ly);
    lx_ = lx;
    ly_ = ly;
    nodes_.clear();
    nodes_.reserve(nodes.size());
    for (const auto &node : nodes) {
      nodes_.emplace_back(node.x, node.y, node.size, node.first_child);
    }
  }

  [[nodiscard]] bool empty() const noexcept
  {
    return nodes_.empty();
  }

  [[nodiscard]] uint32_t root_size() const noexcept
  {
    return root_size_;
  }

  [[nodiscard]] uint32_t num_nodes() const noexcept
  {
    return static_cast<uint32_t>(nodes_.size());
  }

  template <typename Point>
  [[nodiscard]] Leaf locate(const Point &pt) const noexcept
  {
    return locate(pt.x(), pt.y());
  }

  [[nodiscard]] Leaf locate(double px, double py) const noexcept
  {
    const double rx = (px <= 0.0) ? 0.0
                      : (px < double(root_size_))
                        ? px
                        : std::nextafter(double(root_size_), 0.0);
    const double ry = (py <= 0.0) ? 0.0
                      : (py < double(root_size_))
                        ? py
                        : std::nextafter(double(root_size_), 0.0);

    const uint32_t ix = static_cast<uint32_t>(rx);
    const uint32_t iy = static_cast<uint32_t>(ry);
    const uint32_t idx = locate_leaf(ix, iy);
    const Node &n = nodes_[idx];
    return {n.x, n.y, n.size};
  }

  [[nodiscard]] std::vector<Leaf> leaves() const
  {
    std::vector<Leaf> out;
    out.reserve(nodes_.size());
    for (const auto &node : nodes_) {
      if (
        node.is_leaf() && node.x + node.size <= lx_ &&
        node.y + node.size <= ly_) {
        out.emplace_back(node.x, node.y, node.size);
      }
    }
    return out;
  }

private:
  uint32_t locate_leaf(uint32_t px, uint32_t py) const
  {
    uint32_t idx = 0;
    while (!nodes_[idx].is_leaf()) {
      const Node &n = nodes_[idx];
      const uint32_t half = n.size >> 1;
      const bool right = px >= n.x + half;
      const bool bottom = py >= n.y + half;
      idx = static_cast<uint32_t>(n.first_child + (bottom << 1) + right);
    }
    return idx;
  }
  uint32_t lx_{}, ly_{};
  uint32_t root_size_ = 0;
  std::vector<Node> nodes_;
};
