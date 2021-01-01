#include <string>
#include <vector>

#include "pll.h"

PLL::PLL(int pos_,
    Polyline pll_,
    std::vector<int>
    gd_v_,
    std::vector<int> pgnwh_v_,
    bool bool_hole_) {
  pos = pos_;
  pll = pll_;
  gd_v = gd_v_;
  pgnwh_v = pgnwh_v_;

  v1 = *pll_.begin();
  auto v2_it = pll_.begin();
  v2_it++;
  v2 = *v2_it;
  vl = *pll_.rbegin();
  auto v2l_it = pll_.rbegin();
  v2l_it++;
  v2l = *v2l_it;

  bool_hole = bool_hole_;
}

int PLL::get_pos() {
  return pos;
}

Polyline PLL::get_pll() {
  return pll;
}

Point PLL::get_v1() {
  return v1;
}

Point PLL::get_v2() {
  return v2;
}

Point PLL::get_vl() {
  return vl;
}

Point PLL::get_v2l() {
  return v2l;
}

std::vector<int> PLL::get_gd_v() {
  return gd_v;
}

std::vector<int> PLL::get_pgnwh_v() {
  return pgnwh_v;
}

bool PLL::get_bool_hole() {
  return bool_hole;
}

void PLL::set_pll(Polyline pll_) {
  pll = pll_;
  v1 = *pll_.begin();
  auto v2_it = pll_.begin();
  v2_it++;
  v2 = *v2_it;
  vl = *pll_.rbegin();
  auto v2l_it = pll_.rbegin();
  v2l_it++;
  v2l = *v2l_it;
}

void PLL::set_bool_hole(bool b) {
  bool_hole = b;
}
