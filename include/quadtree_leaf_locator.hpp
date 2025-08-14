#pragma once
#include <cmath>
#include <cstdint>
#include <vector>

class QuadtreeLeafLocator
{
public:
  struct Leaf {
    uint32_t x, y, size;
  };

  struct Node {
    uint32_t x, y, size;
    int32_t first_child;
    [[nodiscard]] bool is_leaf() const noexcept
    {
      return first_child < 0;
    }
  };

  QuadtreeLeafLocator() = default;

  template <typename QuadtreeNodeCollection>
  void build(uint32_t root_size, const QuadtreeNodeCollection &nodes)
  {
    root_size_ = root_size;
    nodes_.clear();
    nodes_.reserve(nodes.size());
    for (const auto &node : nodes) {
      nodes_.emplace_back(node.x, node.y, node.size, node.first_child);
    }
  }

  bool empty() const noexcept
  {
    return nodes_.empty();
  }
  uint32_t root_size() const noexcept
  {
    return root_size_;
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

  uint32_t root_size_ = 0;
  std::vector<Node> nodes_;
};
