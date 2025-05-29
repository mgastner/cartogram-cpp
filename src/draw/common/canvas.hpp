#pragma once

#include "colors.hpp"
#include <boost/geometry.hpp>
#include <boost/geometry/io/svg/svg_mapper.hpp>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string_view>

namespace bg = boost::geometry;

using CanvasPoint = bg::model::d2::point_xy<double>;
using CanvasLine = bg::model::linestring<CanvasPoint>;
using CanvasBox = bg::model::box<CanvasPoint>;

class Canvas
{
public:
  Canvas(std::string_view file, double width, double height)
      : out_(std::string(file)), mapper_(out_, width, height)
  {
    mapper_.add(CanvasPoint{0, 0});
    mapper_.add(CanvasPoint{width, height});
  }

  void set_stroke(Color c, double px = 1.0)
  {
    stroke_ = c;
    stroke_w_ = px;
  }

  void set_fill(Color c)
  {
    fill_ = c;
    have_fill_ = true;
  }

  void clear_fill()
  {
    have_fill_ = false;
  }

  void move_to(double x, double y)
  {
    path_.clear();
    path_.push_back({x, y});
  }

  void line_to(double x, double y)
  {
    path_.push_back({x, y});
  }

  void stroke()
  {
    if (path_.size() > 1)
      mapper_.map(path_, stroke_style());
    path_.clear();
  }

  void fill()
  {
    if (have_fill_ && path_.size() > 2) {
      bg::model::polygon<CanvasPoint> poly;
      poly.outer().assign(path_.begin(), path_.end());
      mapper_.map(poly, fill_style());
    }
    path_.clear();
    have_fill_ = false;
  }

  void circle(double x, double y, double r, bool filled = false)
  {
    mapper_.map(CanvasPoint{x, y}, filled ? fill_style() : stroke_style(), r);
  }

  void rectangle(double x, double y, double w, double h, bool filled = false)
  {
    mapper_.map(
      CanvasBox(CanvasPoint{x, y}, CanvasPoint{x + w, y + h}),
      filled ? fill_style() : stroke_style());
  }

  void text(
    double x,
    double y,
    std::string_view t,
    double px,
    std::string_view family = "sans-serif",
    std::string_view style = "normal",  // normal | italic | oblique
    std::string_view weight = "normal")  // normal | bold | 100-900
  {
    std::ostringstream css;
    css << "fill:" << rgb_hex(stroke_) << ";font-family:" << family
        << ";font-style:" << style << ";font-weight:" << weight
        << ";font-size:" << px << "px;" << "text-anchor:start";
    mapper_.text(CanvasPoint{x, y}, std::string(t), css.str());
  }

  void paint_white(double w, double h)
  {
    Color s = stroke_;
    set_fill(Color{1, 1, 1});
    rectangle(0, 0, w, h, true);
    stroke_ = s;
  }

  void fill_polygon(const bg::model::polygon<CanvasPoint> &poly)
  {
    mapper_.map(poly, fill_style());
  }

  void stroke_polygon_outline(const bg::model::polygon<CanvasPoint> &poly)
  {
    // outer
    CanvasLine ls;
    ls.assign(poly.outer().begin(), poly.outer().end());
    mapper_.map(ls, stroke_style());

    // holes
    for (auto const &inner : poly.inners()) {
      CanvasLine hls;
      hls.assign(inner.begin(), inner.end());
      mapper_.map(hls, stroke_style());
    }
  }

private:
  static std::string rgb_hex(Color c)
  {
    std::ostringstream ss;
    ss << '#' << std::hex << std::setfill('0') << std::setw(2)
       << int(c.r * 255) << std::setw(2) << int(c.g * 255) << std::setw(2)
       << int(c.b * 255);
    return ss.str();
  }

  std::string stroke_style() const
  {
    std::ostringstream s;
    s << "stroke:" << rgb_hex(stroke_) << ";stroke-width:" << stroke_w_
      << ";stroke-linecap:round;fill:none";
    return s.str();
  }

  std::string fill_style() const
  {
    std::ostringstream s;
    s << "fill:" << rgb_hex(fill_);
    return s.str();
  }

  std::ofstream out_;
  bg::svg_mapper<CanvasPoint> mapper_;

  CanvasLine path_;
  Color stroke_{0, 0, 0};
  Color fill_{0, 0, 0};
  double stroke_w_{1.0};
  bool have_fill_{false};
};
