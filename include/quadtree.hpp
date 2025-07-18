#pragma once
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>

template <class MetricFn> class Quadtree
{
public:
  struct Leaf {
    uint32_t x, y, size;
  };

  Quadtree(uint32_t root_size, std::size_t target_cells, MetricFn metric)
      : root_size_(root_size), target_(target_cells),
        metric_(std::move(metric))
  {
    assert(
      root_size_ && (root_size_ & (root_size_ - 1)) == 0 &&
      "root_size must be a power of two");

    // Reserve much more than target_cells to account for potential nodes
    // during grading
    const std::size_t cap = target_ * 5 + 8;
    nodes_.reserve(cap);
    leaves_.reserve(cap);
    heap_.reserve(cap);
    leaf_pos_.reserve(cap);

    nodes_.push_back({0, 0, root_size_, 0, -1, metric_(0, 0, root_size_)});
    add_leaf(0);
    heap_.push_back(0);
  }

  void build()
  {
    auto cmp = [this](uint32_t a, uint32_t b) {
      return nodes_[a].priority < nodes_[b].priority;  // max heap
    };
    std::make_heap(heap_.begin(), heap_.end(), cmp);

    // Split nodes until we reach the target number of leaves
    while (leaves_.size() < target_ && !heap_.empty()) {
      std::pop_heap(heap_.begin(), heap_.end(), cmp);
      const uint32_t idx = heap_.back();
      heap_.pop_back();
      Node &n = nodes_[idx];

      // Cannot split further because already at the lowest resolution
      if (!n.is_leaf() || n.size == 1)
        continue;
      split(idx, cmp);
    }
  }

  // Keep splitting the neighbors as long as depth difference is more than 1
  // Time complexity: O(n log n) where n is the target number of leaves
  // Number of new nodes that can be added by grading is bounded by 3 * target
  void grade()
  {
    auto cmp = [this](uint32_t a, uint32_t b) {
      return nodes_[a].priority < nodes_[b].priority;
    };

    bool changed;
    do {
      changed = false;
      const std::size_t stable = leaves_.size();  // only old leaves
      for (std::size_t i = 0; i < stable; ++i) {
        const uint32_t idx = leaves_[i];

        const uint32_t x = nodes_[idx].x;
        const uint32_t y = nodes_[idx].y;
        const uint32_t sz = nodes_[idx].size;

        check_balance(idx, x ? x - 1 : x, y, changed, cmp);
        check_balance(idx, x + sz, y, changed, cmp);
        check_balance(idx, x, y ? y - 1 : y, changed, cmp);
        check_balance(idx, x, y + sz, changed, cmp);
      }
    } while (changed);
  }

  [[nodiscard]] std::vector<Leaf> leaves() const
  {
    std::vector<Leaf> out;
    out.reserve(leaves_.size());
    for (uint32_t idx : leaves_)
      out.push_back({nodes_[idx].x, nodes_[idx].y, nodes_[idx].size});
    return out;
  }

  [[nodiscard]] std::size_t num_leaves() const noexcept
  {
    return leaves_.size();
  }

private:
  struct Node {
    uint32_t x, y, size;
    uint16_t depth;
    int32_t first_child;  // leaf if < 0
    double priority;  // metric value
    [[nodiscard]] bool is_leaf() const noexcept
    {
      return first_child < 0;
    }
  };

  void ensure_leafpos(std::size_t idx)
  {
    if (idx >= leaf_pos_.size())
      leaf_pos_.resize(idx + 1, uint32_t(-1));
  }

  void add_leaf(uint32_t idx)
  {
    ensure_leafpos(idx);
    leaf_pos_[idx] = static_cast<uint32_t>(leaves_.size());
    leaves_.push_back(idx);
  }

  void remove_leaf(uint32_t idx)
  {
    const uint32_t pos = leaf_pos_[idx];
    const uint32_t last = leaves_.back();
    leaves_[pos] = last;
    leaf_pos_[last] = pos;
    leaves_.pop_back();
  }

  template <class Cmp> void heap_push(uint32_t idx, Cmp cmp)
  {
    heap_.push_back(idx);
    std::push_heap(heap_.begin(), heap_.end(), cmp);
  }

  template <class Cmp> void split(uint32_t idx, Cmp cmp)
  {
    Node &p = nodes_[idx];
    remove_leaf(idx);

    p.first_child = static_cast<int32_t>(nodes_.size());
    const uint32_t child_size = p.size >> 1;
    const uint16_t next_depth = p.depth + 1;

    for (uint32_t dy = 0; dy < 2; ++dy)
      for (uint32_t dx = 0; dx < 2; ++dx) {
        const uint32_t cx = p.x + dx * child_size;
        const uint32_t cy = p.y + dy * child_size;
        const double pr = metric_(cx, cy, child_size);

        nodes_.push_back({cx, cy, child_size, next_depth, -1, pr});
        const uint32_t cid = static_cast<uint32_t>(nodes_.size() - 1);
        add_leaf(cid);
        heap_push(cid, cmp);
      }
  }

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

  template <class Cmp>
  void check_balance(
    uint32_t idx,
    uint32_t qx,
    uint32_t qy,
    bool &changed,
    Cmp cmp)
  {
    const uint32_t nb = locate_leaf(qx, qy);
    if (nb == idx)
      return;

    const Node &a = nodes_[idx];
    const Node &b = nodes_[nb];
    if (std::abs(int(a.depth) - int(b.depth)) > 1) {
      const uint32_t shallow = (a.depth < b.depth) ? idx : nb;
      if (nodes_[shallow].size > 1) {
        split(shallow, cmp);
        changed = true;
      }
    }
  }

  uint32_t root_size_;
  std::size_t target_;
  MetricFn metric_;

  std::vector<Node> nodes_;
  std::vector<uint32_t> leaves_;
  std::vector<uint32_t> heap_;
  std::vector<uint32_t> leaf_pos_;
};
