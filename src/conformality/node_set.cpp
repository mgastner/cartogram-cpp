#include "conformality/node_set.hpp"
#include <stdexcept>

namespace conformality
{

std::size_t NodeSet::size() const noexcept
{
  return x.size();
}

void NodeSet::check_consistent() const
{
  const std::size_t n = x.size();
  if (w.size() != n || g0.size() != n || y.size() != n || jac.size() != n) {
    throw std::invalid_argument("NodeSet: arrays differ in length");
  }
}

NodeSet make_cell_centre_nodes(
  const double lx,
  const double ly,
  const std::size_t nx,
  const std::size_t ny)
{
  if (!(lx > 0.0) || !(ly > 0.0) || nx == 0 || ny == 0) {
    throw std::invalid_argument(
      "make_cell_centre_nodes: sizes must be positive");
  }
  const double hx = lx / static_cast<double>(nx);
  const double hy = ly / static_cast<double>(ny);
  const std::size_t n = nx * ny;
  NodeSet nodes;
  nodes.x.reserve(n);
  for (std::size_t i = 0; i < nx; ++i) {
    for (std::size_t j = 0; j < ny; ++j) {
      nodes.x.push_back(
        {(static_cast<double>(i) + 0.5) * hx,
         (static_cast<double>(j) + 0.5) * hy});
    }
  }
  nodes.w.assign(n, hx * hy);
  nodes.g0.assign(n, 1.0);
  nodes.y = nodes.x;
  nodes.jac.assign(n, identity_mat2());
  return nodes;
}

}  // namespace conformality
